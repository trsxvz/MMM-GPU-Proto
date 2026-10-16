/*!
 * \file   src/TrilinosBackend.cxx
 * \brief  This file implements the interface of the library of mfem-mgis
 * built against Trilinos (see `TrilinosBackend.hxx`), compiled by the
 * compiler of Kokkos
 * \date   15/10/2026
 */

#include <array>
#include <cmath>
#include <memory>
#include <string>
#include <vector>
#include <variant>
#include <utility>
#include <algorithm>
#include <stdexcept>
#include <type_traits>
#include "Kokkos_Core.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_DefaultMpiComm.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Operator.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosSolverFactory.hpp"
#include "BelosTpetraAdapter.hpp"
#include "Ifpack2_Factory.hpp"
#include "MueLu_CreateTpetraPreconditioner.hpp"
#include "Amesos2.hpp"
#include "Xpetra_CrsMatrixWrap.hpp"
#include "Xpetra_TpetraCrsMatrix.hpp"
#include "Xpetra_TpetraMultiVector.hpp"
#include "FROSch_TwoLevelPreconditioner_def.hpp"
#include "Ifpack2_Parameters.hpp"
#include "MueLu_MasterList.hpp"
#include "FROSch_TpetraPreconditioner_def.hpp"
#include "TrilinosBackend.hxx"

namespace mfem_mgis::trilinos {

  //! \brief scalar type
  using Scalar = double;
  //! \brief map of Tpetra (default ordinals and node)
  using Map = Tpetra::Map<>;
  //! \brief local ordinal
  using LO = Map::local_ordinal_type;
  //! \brief global ordinal
  using GO = Map::global_ordinal_type;
  //! \brief node of Tpetra
  using Node = Map::node_type;
  //! \brief sparse matrix of Tpetra
  using CrsMatrix = Tpetra::CrsMatrix<Scalar, LO, GO, Node>;
  //! \brief row matrix of Tpetra
  using RowMatrix = Tpetra::RowMatrix<Scalar, LO, GO, Node>;
  //! \brief multivector of Tpetra
  using MV = Tpetra::MultiVector<Scalar, LO, GO, Node>;
  //! \brief operator of Tpetra
  using OP = Tpetra::Operator<Scalar, LO, GO, Node>;
  //! \brief execution space of Tpetra
  using ExecutionSpace = Node::execution_space;
  //! \brief memory space of Tpetra
  using MemorySpace = Node::memory_space;
  //! \brief view of an array of the memory space of Tpetra
  template <typename T>
  using DeviceArray = Kokkos::View<T*, MemorySpace>;
  //! \brief view of an array of the memory space of Tpetra, not owned
  template <typename T>
  using UnmanagedDeviceArray =
      Kokkos::View<T*, MemorySpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  double Operator::getRelativeTolerance() const { return -1; }

  bool Operator::setRelativeTolerance(const double) { return false; }

  Operator::~Operator() = default;

  bool computesOnTheGPU() noexcept {
    return !Kokkos::SpaceAccessibility<Node::execution_space,
                                       Kokkos::HostSpace>::accessible;
  }  // end of computesOnTheGPU

  /*!
   * \brief check that Kokkos is initialized (by mfem-mgis, see
   * `initializeKokkosOnDemand`, or by the program)
   */
  static void checkKokkosInitialization() {
    if ((!Kokkos::is_initialized()) || (Kokkos::is_finalized())) {
      throw std::runtime_error("Trilinos: Kokkos is not initialized");
    }
  }  // end of checkKokkosInitialization

  /*!
   * \brief fill a parameter list of Teuchos with options
   * \param[out] l: list
   * \param[in] options: options
   */
  static void fill(Teuchos::ParameterList& l, const Options& options) {
    for (const auto& [n, o] : options) {
      std::visit(
          [&l, &n](const auto& v) {
            using Type = std::decay_t<decltype(v)>;
            if constexpr (std::is_same_v<Type, Options>) {
              fill(l.sublist(n), v);
            } else {
              l.set(n, v);
            }
          },
          o.value);
    }
  }  // end of fill

  /*!
   * \return a parameter list of Teuchos
   * \param[in] options: options
   * \param[in] name: name of the list (Amesos2 only reads a list named
   * `Amesos2`)
   */
  static Teuchos::RCP<Teuchos::ParameterList> makeParameterList(
      const Options& options, const std::string& name = "ANONYMOUS") {
    auto l = Teuchos::rcp(new Teuchos::ParameterList(name));
    fill(*l, options);
    return l;
  }  // end of makeParameterList

  /*!
   * \return the local index, in the ordering by vector dimension, of an
   * unknown of a process ordered by nodes
   * \param[in] i: local index of the unknown, ordered by nodes
   * \param[in] n: number of unknowns of the process
   * \param[in] c: number of components
   */
  static long long interleave(const long long i,
                              const long long n,
                              const long long c) {
    const auto nn = n / c;
    return (i % nn) * c + i / nn;
  }  // end of interleave

  /*!
   * \brief sparse matrix of Tpetra built from a matrix of hypre: the row map
   * is the contiguous map of the rows of the processes, the column map the
   * rows of the process followed by the columns of the off-diagonal block.
   * The unknowns ordered by nodes with several components are permuted in
   * each process (see `DistributedMatrixView`).
   */
  struct TpetraMatrix {
    /*!
     * \brief build the matrix
     * \param[in] m: matrix of hypre
     */
    void set(const DistributedMatrixView& m) {
      checkKokkosInitialization();
      this->communicator = m.communicator;
      const auto comm = Teuchos::rcp(new Teuchos::MpiComm<int>(m.communicator));
      const auto n = m.local_size;
      const auto no = m.number_of_offd_columns;
      const auto vdim = m.number_of_components;
      // first rows and numbers of rows of all the processes, checked on all
      // the processes, so that they all throw on error
      const auto nprocs = static_cast<std::size_t>(comm->getSize());
      auto rows = std::vector<long long>(2 * nprocs);
      const auto local = std::array<long long, 2>{m.first_row, n};
      MPI_Allgather(local.data(), 2, MPI_LONG_LONG, rows.data(), 2,
                    MPI_LONG_LONG, m.communicator);
      for (std::size_t p = 0; p != nprocs; ++p) {
        if ((vdim < 1) || (rows[2 * p + 1] % vdim != 0)) {
          throw std::runtime_error(
              "Trilinos: the number of rows of a process is not a multiple "
              "of the number of components of the unknowns");
        }
      }
      const auto permuted = (vdim > 1) && (m.ordered_by_nodes);
      // permutation of the rows of the process
      this->permutation.clear();
      this->device_permutation = DeviceArray<LO>{};
      if (permuted) {
        this->permutation.resize(static_cast<std::size_t>(n));
        for (int i = 0; i != n; ++i) {
          this->permutation[static_cast<std::size_t>(i)] =
              static_cast<LO>(interleave(i, n, vdim));
        }
        this->device_permutation = DeviceArray<LO>(
            Kokkos::view_alloc(Kokkos::WithoutInitializing, "permutation"),
            static_cast<std::size_t>(n));
        Kokkos::deep_copy(
            this->device_permutation,
            Kokkos::View<const LO*, Kokkos::HostSpace,
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
                this->permutation.data(), this->permutation.size()));
      }
      // columns: the rows of the process, then the off-diagonal columns,
      // permuted according to the rows of the processes owning them
      auto columns = std::vector<GO>(static_cast<std::size_t>(n + no));
      for (int i = 0; i != n; ++i) {
        columns[static_cast<std::size_t>(i)] = static_cast<GO>(m.first_row + i);
      }
      // first rows of the processes, in increasing order
      auto starts = std::vector<long long>(nprocs);
      for (std::size_t p = 0; p != nprocs; ++p) {
        starts[p] = rows[2 * p];
      }
      for (int k = 0; k != no; ++k) {
        const auto g = m.offd_global_columns[k];
        auto c = g;
        if (permuted) {
          // process owning the column: the last one whose first row is not
          // greater than the column (the processes without rows have the
          // first row of the next process)
          const auto p = static_cast<std::size_t>(
              std::upper_bound(starts.begin(), starts.end(), g) -
              starts.begin() - 1);
          c = rows[2 * p] + interleave(g - rows[2 * p], rows[2 * p + 1], vdim);
        }
        columns[static_cast<std::size_t>(n + k)] = static_cast<GO>(c);
      }
      this->row_map = Teuchos::rcp(
          new Map(static_cast<Tpetra::global_size_t>(m.global_size),
                  static_cast<std::size_t>(n), GO{0}, comm));
      const auto column_map = Teuchos::rcp(
          new Map(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),
                  Teuchos::ArrayView<const GO>(
                      columns.data(), static_cast<int>(columns.size())),
                  GO{0}, comm));
      // local matrix, the columns of each row being sorted
      using LocalMatrix = CrsMatrix::local_matrix_device_type;
      using Offsets = CrsMatrix::local_graph_device_type::row_map_type;
      using Entries = CrsMatrix::local_graph_device_type::entries_type;
      using Values = LocalMatrix::values_type;
      const auto has_offd = (no != 0) && (m.offd_offsets != nullptr);
      const auto nnz =
          n == 0 ? std::size_t{0}
                 : static_cast<std::size_t>(m.diag_offsets[n]) +
                       (has_offd ? static_cast<std::size_t>(m.offd_offsets[n])
                                 : std::size_t{0});
      auto offsets = typename Offsets::non_const_type(
          Kokkos::view_alloc(Kokkos::WithoutInitializing, "offsets"),
          static_cast<std::size_t>(n) + 1);
      auto entries = typename Entries::non_const_type(
          Kokkos::view_alloc(Kokkos::WithoutInitializing, "entries"), nnz);
      auto values = typename Values::non_const_type(
          Kokkos::view_alloc(Kokkos::WithoutInitializing, "values"), nnz);
      auto h_offsets = Kokkos::create_mirror_view(offsets);
      auto h_entries = Kokkos::create_mirror_view(entries);
      auto h_values = Kokkos::create_mirror_view(values);
      // rows of the process in the order of Tpetra
      auto source = std::vector<int>(static_cast<std::size_t>(n));
      for (int i = 0; i != n; ++i) {
        source[static_cast<std::size_t>(this->getIndex(i))] = i;
      }
      auto row = std::vector<std::pair<LO, Scalar>>{};
      auto k = std::size_t{0};
      for (int r = 0; r != n; ++r) {
        const auto i = source[static_cast<std::size_t>(r)];
        h_offsets(static_cast<std::size_t>(r)) = k;
        row.clear();
        for (auto j = m.diag_offsets[i]; j != m.diag_offsets[i + 1]; ++j) {
          row.emplace_back(this->getIndex(m.diag_columns[j]),
                           static_cast<Scalar>(m.diag_values[j]));
        }
        if (has_offd) {
          for (auto j = m.offd_offsets[i]; j != m.offd_offsets[i + 1]; ++j) {
            row.emplace_back(static_cast<LO>(n + m.offd_columns[j]),
                             static_cast<Scalar>(m.offd_values[j]));
          }
        }
        std::sort(row.begin(), row.end(), [](const auto& a, const auto& b) {
          return a.first < b.first;
        });
        for (const auto& [c, v] : row) {
          h_entries(k) = c;
          h_values(k) = v;
          ++k;
        }
      }
      h_offsets(static_cast<std::size_t>(n)) = k;
      Kokkos::deep_copy(offsets, h_offsets);
      Kokkos::deep_copy(entries, h_entries);
      Kokkos::deep_copy(values, h_values);
      this->matrix = Teuchos::rcp(
          new CrsMatrix(this->row_map, column_map, offsets, entries, values));
      this->matrix->fillComplete(this->row_map, this->row_map);
    }  // end of set
    /*!
     * \return the local index for Tpetra of an unknown of the process
     * \param[in] i: local index in the matrix of hypre
     */
    [[nodiscard]] LO getIndex(const int i) const noexcept {
      return this->permutation.empty()
                 ? static_cast<LO>(i)
                 : this->permutation[static_cast<std::size_t>(i)];
    }
    //! \brief communicator
    MPI_Comm communicator = MPI_COMM_NULL;
    //! \brief map of the rows
    Teuchos::RCP<const Map> row_map;
    //! \brief matrix
    Teuchos::RCP<CrsMatrix> matrix;
    /*!
     * \brief local indices for Tpetra of the rows of the process, empty if
     * the rows are not permuted
     */
    std::vector<LO> permutation;
    //! \brief copy of the permutation in the memory space of Tpetra
    DeviceArray<LO> device_permutation;
  };  // end of struct TpetraMatrix

  /*!
   * \brief checks that the columns of a multivector are contiguous in its
   * local views (the copies below index the columns of these views)
   * \param[in] v: multivector
   */
  static void checkConstantStride(const MV& v) {
    if (!v.isConstantStride()) {
      throw std::runtime_error(
          "Trilinos: multivectors without constant stride are not "
          "supported");
    }
  }  // end of checkConstantStride

  /*!
   * \brief copy a vector of the host in the column of a multivector
   * \param[out] d: multivector
   * \param[in] j: column
   * \param[in] s: vector
   * \param[in] m: matrix, giving the permutation of the unknowns
   */
  static void copy(MV& d,
                   const std::size_t j,
                   const double* s,
                   const TpetraMatrix& m) {
    checkConstantStride(d);
    const auto fill = [j, s, &m](const auto& v) {
      for (std::size_t i = 0; i != v.extent(0); ++i) {
        v(static_cast<std::size_t>(m.getIndex(static_cast<int>(i))), j) = s[i];
      }
    };
    // the values of a multivector of one column are all overwritten: they
    // are not copied from the device first
    if (d.getNumVectors() == 1) {
      fill(d.getLocalViewHost(Tpetra::Access::OverwriteAll));
    } else {
      fill(d.getLocalViewHost(Tpetra::Access::ReadWrite));
    }
  }  // end of copy

  /*!
   * \brief copy the column of a multivector in a vector of the host
   * \param[out] d: vector
   * \param[in] s: multivector
   * \param[in] j: column
   * \param[in] m: matrix, giving the permutation of the unknowns
   */
  static void copy(double* d,
                   const MV& s,
                   const std::size_t j,
                   const TpetraMatrix& m) {
    checkConstantStride(s);
    const auto v = s.getLocalViewHost(Tpetra::Access::ReadOnly);
    for (std::size_t i = 0; i != v.extent(0); ++i) {
      d[i] = v(static_cast<std::size_t>(m.getIndex(static_cast<int>(i))), j);
    }
  }  // end of copy

  /*!
   * \brief copy a vector of the memory space of Tpetra in the column of a
   * view of a multivector (a named function: nvcc does not accept the device
   * lambdas inside generic lambdas)
   * \param[out] v: view of the multivector
   * \param[in] j: column
   * \param[in] s: vector
   * \param[in] m: matrix, giving the permutation of the unknowns
   */
  template <typename View>
  static void copyToColumnOnDevice(const View& v,
                                   const std::size_t j,
                                   const double* s,
                                   const TpetraMatrix& m) {
    const auto n = v.extent(0);
    const auto src = UnmanagedDeviceArray<const double>(s, n);
    if (m.permutation.empty()) {
      Kokkos::deep_copy(Kokkos::subview(v, Kokkos::ALL, j), src);
      return;
    }
    const auto p = m.device_permutation;
    Kokkos::parallel_for(
        "mfem_mgis::trilinos::copyOnDevice",
        Kokkos::RangePolicy<ExecutionSpace>(0, n),
        KOKKOS_LAMBDA(const std::size_t i) {
          v(static_cast<std::size_t>(p(i)), j) = src(i);
        });
  }  // end of copyToColumnOnDevice

  /*!
   * \brief copy a vector of the memory space of Tpetra in the column of a
   * multivector
   * \param[out] d: multivector
   * \param[in] j: column
   * \param[in] s: vector
   * \param[in] m: matrix, giving the permutation of the unknowns
   */
  static void copyOnDevice(MV& d,
                           const std::size_t j,
                           const double* s,
                           const TpetraMatrix& m) {
    checkConstantStride(d);
    // the values of a multivector of one column are all overwritten
    if (d.getNumVectors() == 1) {
      copyToColumnOnDevice(d.getLocalViewDevice(Tpetra::Access::OverwriteAll),
                           j, s, m);
    } else {
      copyToColumnOnDevice(d.getLocalViewDevice(Tpetra::Access::ReadWrite), j,
                           s, m);
    }
  }  // end of copyOnDevice

  /*!
   * \brief copy the column of a multivector in a vector of the memory space
   * of Tpetra
   * \param[out] d: vector
   * \param[in] s: multivector
   * \param[in] j: column
   * \param[in] m: matrix, giving the permutation of the unknowns
   */
  static void copyOnDevice(double* d,
                           const MV& s,
                           const std::size_t j,
                           const TpetraMatrix& m) {
    checkConstantStride(s);
    const auto v = s.getLocalViewDevice(Tpetra::Access::ReadOnly);
    const auto n = v.extent(0);
    const auto dst = UnmanagedDeviceArray<double>(d, n);
    if (m.permutation.empty()) {
      Kokkos::deep_copy(dst, Kokkos::subview(v, Kokkos::ALL, j));
      return;
    }
    const auto p = m.device_permutation;
    Kokkos::parallel_for(
        "mfem_mgis::trilinos::copyOnDevice",
        Kokkos::RangePolicy<ExecutionSpace>(0, n),
        KOKKOS_LAMBDA(const std::size_t i) {
          dst(i) = v(static_cast<std::size_t>(p(i)), j);
        });
  }  // end of copyOnDevice

  /*!
   * \brief base class of the operators of Tpetra applied to vectors of the
   * host or of the memory space of Tpetra: the operator of Tpetra may be
   * given to Belos
   */
  struct TpetraOperatorBase : public Operator {
    void apply(const double* b, double* x) const override {
      this->allocateVectors();
      copy(*(this->b), 0, b, this->A);
      copy(*(this->x), 0, x, this->A);
      this->applyTpetraOperator(*(this->b), *(this->x));
      copy(x, *(this->x), 0, this->A);
    }
    void applyOnDevice(const double* b, double* x) const override {
      if (!computesOnTheGPU()) {
        this->apply(b, x);
        return;
      }
      this->allocateVectors();
      copyOnDevice(*(this->b), 0, b, this->A);
      copyOnDevice(*(this->x), 0, x, this->A);
      this->applyTpetraOperator(*(this->b), *(this->x));
      copyOnDevice(x, *(this->x), 0, this->A);
      Kokkos::fence("mfem_mgis::trilinos::applyOnDevice");
    }
    [[nodiscard]] Report getReport() const override {
      return {.converged = true,
              .iterations = 0,
              .achieved_tolerance = 0,
              .failure = {}};
    }
    //! \return the matrix
    [[nodiscard]] const TpetraMatrix& getMatrix() const noexcept {
      return this->A;
    }
    //! \return the operator of Tpetra, null if not set up
    [[nodiscard]] virtual Teuchos::RCP<OP> getTpetraOperator() const = 0;
    /*!
     * \brief apply the operator of Tpetra: `x` is overwritten, except by
     * the iterative solvers, which use it as initial guess
     * \param[in] b: right hand side
     * \param[in, out] x: result
     */
    virtual void applyTpetraOperator(const MV& b, MV& x) const {
      this->getTpetraOperator()->apply(b, x);
    }

   protected:
    //! \brief matrix
    TpetraMatrix A;

   private:
    /*!
     * \brief allocate the right hand side and the solution, on the maps of
     * the operator of Tpetra
     */
    void allocateVectors() const {
      const auto op = this->getTpetraOperator();
      if (op.is_null()) {
        throw std::runtime_error("Trilinos: the matrix is not set");
      }
      if ((this->b.is_null()) ||
          (this->b->getMap().get() != op->getRangeMap().get())) {
        this->b = Teuchos::rcp(new MV(op->getRangeMap(), 1));
        this->x = Teuchos::rcp(new MV(op->getDomainMap(), 1));
      }
    }
    //! \brief right hand side
    mutable Teuchos::RCP<MV> b;
    //! \brief solution
    mutable Teuchos::RCP<MV> x;
  };  // end of struct TpetraOperatorBase

  /*!
   * \brief checks that the names of the first level of a list of parameters
   * are known, so that a misspelt parameter is refused at the creation of
   * the solver (Ifpack2 and MueLu check them at the setup, or ignore them)
   * \param[in] library: name of the library
   * \param[in] p: parameters
   * \param[in] valid: known parameters
   * \param[in] prefixes: prefixes of the other names accepted
   */
  static void checkParameterNames(const std::string& library,
                                  const Teuchos::ParameterList& p,
                                  const Teuchos::ParameterList& valid,
                                  const std::vector<std::string>& prefixes) {
    for (auto i = p.begin(); i != p.end(); ++i) {
      const auto& n = p.name(i);
      const auto prefixed = std::ranges::any_of(
          prefixes,
          [&n](const std::string& prefix) { return n.starts_with(prefix); });
      if ((!valid.isParameter(n)) && (!prefixed)) {
        throw std::runtime_error(library + ": unknown parameter '" + n + "'");
      }
    }
  }  // end of checkParameterNames

  //! \brief preconditioner of Ifpack2
  struct Ifpack2Preconditioner final : public TpetraOperatorBase {
    /*!
     * \param[in] t: type
     * \param[in] o: options
     */
    Ifpack2Preconditioner(const std::string& t, const Options& o)
        : type(t), parameters(makeParameterList(o)) {
      auto valid = Teuchos::ParameterList{};
      Ifpack2::getValidParameters(valid);
      checkParameterNames("Ifpack2", *(this->parameters), valid, {});
    }
    void setMatrix(const DistributedMatrixView& m) override {
      this->A.set(m);
      auto p = Ifpack2::Factory::create<RowMatrix>(
          this->type, Teuchos::RCP<const RowMatrix>(this->A.matrix));
      if (p.is_null()) {
        throw std::runtime_error("invalid preconditioner '" + this->type + "'");
      }
      p->setParameters(*(this->parameters));
      p->initialize();
      p->compute();
      this->preconditioner = p;
    }
    [[nodiscard]] Teuchos::RCP<OP> getTpetraOperator() const override {
      return this->preconditioner;
    }

   private:
    //! \brief type
    const std::string type;
    //! \brief parameters
    const Teuchos::RCP<Teuchos::ParameterList> parameters;
    //! \brief preconditioner
    Teuchos::RCP<Ifpack2::Preconditioner<Scalar, LO, GO, Node>> preconditioner;
  };  // end of struct Ifpack2Preconditioner

  /*!
   * \return the number of nodes of the process, checked against the
   * coordinates on all the processes, which all throw on error
   * \param[in] comm: communicator
   * \param[in] map: map of the unknowns
   * \param[in] vdim: number of components of the unknowns
   * \param[in] c: coordinates of the nodes
   */
  static std::size_t getNumberOfNodes(const MPI_Comm comm,
                                      const Map& map,
                                      const int vdim,
                                      const NodesCoordinates& c) {
    const auto nn = map.getLocalNumElements() / static_cast<std::size_t>(vdim);
    auto mismatch = static_cast<int>(
        c.values.size() != static_cast<std::size_t>(c.dimension) * nn);
    MPI_Allreduce(MPI_IN_PLACE, &mismatch, 1, MPI_INT, MPI_LOR, comm);
    if (mismatch != 0) {
      throw std::runtime_error(
          "the number of nodes of the coordinates does not match the number "
          "of unknowns");
    }
    return nn;
  }  // end of getNumberOfNodes

  /*!
   * \return the coordinates of the nodes, on the contiguous map of the nodes
   * (the unknowns of a node being contiguous, see `TpetraMatrix`), or a null
   * pointer if the coordinates are not given (number of components of the
   * unknowns different from the space dimension)
   * \param[in] m: matrix
   * \param[in] vdim: number of components of the unknowns
   * \param[in] c: coordinates of the nodes
   */
  static Teuchos::RCP<MV> makeNodesCoordinates(const TpetraMatrix& m,
                                               const int vdim,
                                               const NodesCoordinates& c) {
    const auto d = c.dimension;
    if ((vdim != d) || ((d != 2) && (d != 3))) {
      return Teuchos::null;
    }
    const auto nn = getNumberOfNodes(m.communicator, *(m.row_map), vdim, c);
    const auto map = Teuchos::rcp(new Map(
        m.row_map->getGlobalNumElements() / static_cast<std::size_t>(vdim), nn,
        GO{0}, m.row_map->getComm()));
    auto x = Teuchos::rcp(new MV(map, static_cast<std::size_t>(d)));
    auto v = x->getLocalViewHost(Tpetra::Access::OverwriteAll);
    for (std::size_t k = 0; k != static_cast<std::size_t>(d); ++k) {
      for (std::size_t node = 0; node != nn; ++node) {
        v(node, k) = c.values[k * nn + node];
      }
    }
    return x;
  }  // end of makeNodesCoordinates

  /*!
   * \return the rigid body modes of the unknowns of the process, ordered by
   * vector dimension, or a null pointer if the number of components of the
   * unknowns is not the space dimension (2 or 3) or if the coordinates are
   * not given
   * \param[in] comm: communicator
   * \param[in] map: map of the unknowns
   * \param[in] vdim: number of components of the unknowns
   * \param[in] c: coordinates of the nodes
   */
  static Teuchos::RCP<MV> makeRigidBodyModes(const MPI_Comm comm,
                                             const Teuchos::RCP<const Map>& map,
                                             const int vdim,
                                             const NodesCoordinates& c) {
    const auto d = c.dimension;
    if ((vdim != d) || ((d != 2) && (d != 3))) {
      return Teuchos::null;
    }
    const auto nn = getNumberOfNodes(comm, *map, vdim, c);
    // coordinates relative to their mean over all the processes, for the
    // conditioning of the rotations (the same modes on all the processes)
    auto sums = std::array<double, 4>{0, 0, 0, static_cast<double>(nn)};
    for (std::size_t k = 0; k != static_cast<std::size_t>(d); ++k) {
      for (std::size_t node = 0; node != nn; ++node) {
        sums[k] += c.values[k * nn + node];
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, sums.data(), 4, MPI_DOUBLE, MPI_SUM, comm);
    auto mean = std::array<double, 3>{0, 0, 0};
    for (std::size_t k = 0; k != static_cast<std::size_t>(d); ++k) {
      mean[k] = sums[k] / std::max(sums[3], 1.);
    }
    const auto nmodes = std::size_t{d == 2 ? 3u : 6u};
    auto modes = Teuchos::rcp(new MV(map, nmodes));
    auto v = modes->getLocalViewHost(Tpetra::Access::OverwriteAll);
    for (std::size_t node = 0; node != nn; ++node) {
      auto x = std::array<double, 3>{0, 0, 0};
      for (std::size_t k = 0; k != static_cast<std::size_t>(d); ++k) {
        x[k] = c.values[k * nn + node] - mean[k];
      }
      // rotations around z (-y, x, 0), around x (0, -z, y) and around y
      // (z, 0, -x)
      const auto rz = std::array<double, 3>{-x[1], x[0], 0};
      const auto rx = std::array<double, 3>{0, -x[2], x[1]};
      const auto ry = std::array<double, 3>{x[2], 0, -x[0]};
      for (std::size_t k = 0; k != static_cast<std::size_t>(vdim); ++k) {
        const auto i = node * static_cast<std::size_t>(vdim) + k;
        for (std::size_t mode = 0; mode != nmodes; ++mode) {
          v(i, mode) = 0;
        }
        v(i, k) = 1;
        if (d == 2) {
          v(i, 2) = rz[k];
        } else {
          v(i, 3) = rz[k];
          v(i, 4) = rx[k];
          v(i, 5) = ry[k];
        }
      }
    }
    return modes;
  }  // end of makeRigidBodyModes

  //! \brief smoothed aggregation multigrid of MueLu
  struct MueLuPreconditioner final : public TpetraOperatorBase {
    /*!
     * \param[in] o: options
     * \param[in] c: coordinates of the nodes
     */
    MueLuPreconditioner(const Options& o, NodesCoordinates c)
        : parameters(makeParameterList(o)), coordinates(std::move(c)) {
      // the parameters of a level are given in a sublist `level N`
      checkParameterNames("MueLu", *(this->parameters),
                          *(MueLu::MasterList::List()), {"level "});
    }
    void setMatrix(const DistributedMatrixView& m) override {
      this->A.set(m);
      auto p = *(this->parameters);
      const auto vdim = m.number_of_components;
      const auto* const equations = "number of equations";
      if ((vdim > 1) && (!p.isParameter(equations))) {
        p.set(equations, vdim);
      }
      // the rigid body modes, unless the parameters give another number of
      // equations per node
      if (p.get<int>(equations, 1) == vdim) {
        const auto modes = makeRigidBodyModes(m.communicator, this->A.row_map,
                                              vdim, this->coordinates);
        if (!modes.is_null()) {
          p.sublist("user data").set("Nullspace", modes);
        }
      }
      // products of sparse matrices by the algorithm of Kokkos Kernels on
      // CUDA: the default one (cuSPARSE) ran out of memory on a GPU of 16 GB
      // from 421,443 unknowns of linear elasticity
      auto& products = p.sublist("matrixmatrix: kernel params");
      if (!products.isParameter("cuda: algorithm")) {
        products.set("cuda: algorithm", "SPGEMM_KK");
      }
      // reproducible aggregates (MueLu aggregates the nodes in a random
      // order by default), and silent setup
      if (!p.isParameter("aggregation: deterministic")) {
        p.set("aggregation: deterministic", true);
      }
      if (!p.isParameter("verbosity")) {
        p.set("verbosity", "none");
      }
      // coordinates of the nodes, used by some options of MueLu (dropping
      // based on the distances, repartitioning), when the nodes of MueLu are
      // the nodes of the mesh
      if (p.get<int>(equations, 1) == vdim) {
        const auto coordinates =
            makeNodesCoordinates(this->A, vdim, this->coordinates);
        if (!coordinates.is_null()) {
          p.sublist("user data").set("Coordinates", coordinates);
        }
      }
      this->preconditioner = MueLu::CreateTpetraPreconditioner(
          Teuchos::RCP<OP>(this->A.matrix), p);
    }
    [[nodiscard]] Teuchos::RCP<OP> getTpetraOperator() const override {
      return this->preconditioner;
    }

   private:
    //! \brief parameters
    const Teuchos::RCP<Teuchos::ParameterList> parameters;
    //! \brief coordinates of the nodes
    const NodesCoordinates coordinates;
    //! \brief preconditioner
    Teuchos::RCP<OP> preconditioner;
  };  // end of struct MueLuPreconditioner

  /*!
   * \brief makes the setup of FROSch silent (it prints tables at each setup
   * by default), unless a verbosity is given: the verbosity is read in the
   * parameters of the overlapping operator and of the first block of the
   * coarse operator
   * \param[in, out] p: parameters of FROSch
   */
  static void setDefaultFROSchVerbosity(Teuchos::ParameterList& p) {
    auto silent = [](Teuchos::ParameterList& l) {
      if (!l.isParameter("Verbosity")) {
        l.set("Verbosity", "None");
      }
    };
    silent(p.sublist("AlgebraicOverlappingOperator"));
    const auto coarse = p.isParameter("CoarseOperator Type")
                            ? p.get<std::string>("CoarseOperator Type")
                            : std::string{"IPOUHarmonicCoarseOperator"};
    silent(p.sublist(coarse).sublist("Blocks").sublist("1"));
  }  // end of setDefaultFROSchVerbosity

  //! \brief two level overlapping Schwarz preconditioner of FROSch
  struct FROSchPreconditioner final : public TpetraOperatorBase {
    /*!
     * \param[in] o: options
     * \param[in] c: coordinates of the nodes
     */
    FROSchPreconditioner(const Options& o, NodesCoordinates c)
        : parameters(makeParameterList(o)), coordinates(std::move(c)) {
      if (this->parameters->isParameter(
              "DirichletBoundary Numerical Tolerance")) {
        throw std::runtime_error(
            "FROSch: the parameter 'DirichletBoundary Numerical Tolerance' "
            "is not supported (the rows of the Dirichlet boundary "
            "conditions are not given to FROSch)");
      }
    }
    void setMatrix(const DistributedMatrixView& m) override {
      using XpetraCrsMatrix = Xpetra::TpetraCrsMatrix<Scalar, LO, GO, Node>;
      using XpetraMatrix = Xpetra::CrsMatrixWrap<Scalar, LO, GO, Node>;
      using XpetraMultiVector = Xpetra::TpetraMultiVector<Scalar, LO, GO, Node>;
      // FROSch 17.2.1 crashes when a process owns no unknown (constant
      // partition of unity): checked by all the processes, which fail
      // together
      auto empty = static_cast<int>(m.local_size == 0);
      MPI_Allreduce(MPI_IN_PLACE, &empty, 1, MPI_INT, MPI_LOR, m.communicator);
      if (empty != 0) {
        throw std::runtime_error(
            "FROSch: each process must own unknowns (FROSch does not support "
            "empty subdomains)");
      }
      this->A.set(m);
      const auto vdim = m.number_of_components;
      const auto d = this->coordinates.dimension;
      if ((d != 2) && (d != 3)) {
        throw std::runtime_error("FROSch: the space dimension must be 2 or 3");
      }
      auto p = Teuchos::rcp(new Teuchos::ParameterList(*(this->parameters)));
      // coordinates of the nodes, given when the number of components of the
      // unknowns is the space dimension
      const auto x = makeNodesCoordinates(this->A, vdim, this->coordinates);
      auto nodes =
          Teuchos::RCP<const Xpetra::MultiVector<Scalar, LO, GO, Node>>{};
      if (!x.is_null()) {
        nodes = Teuchos::rcp(new XpetraMultiVector(x));
      }
      if (!p->isParameter("Null Space Type")) {
        p->set("Null Space Type",
               x.is_null() ? "Laplace" : "Linear Elasticity");
      }
      // the values of the overlaps are added (additive Schwarz method), so
      // that the preconditioner is symmetric (restricted additive Schwarz
      // method by default in FROSch)
      auto& overlapping = p->sublist("AlgebraicOverlappingOperator");
      if (!overlapping.isParameter("Combine Values in Overlap")) {
        overlapping.set("Combine Values in Overlap", "Full");
      }
      setDefaultFROSchVerbosity(*p);
      // the conversions of the pointers are explicit, the constructors of
      // Xpetra::CrsMatrixWrap taking several kinds of pointers
      const auto xc = Teuchos::RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node>>(
          Teuchos::rcp(new XpetraCrsMatrix(this->A.matrix)));
      const auto xa = Teuchos::RCP<const Xpetra::Matrix<Scalar, LO, GO, Node>>(
          Teuchos::rcp(new XpetraMatrix(xc)));
      // the rows of the Dirichlet boundary conditions are not given to
      // FROSch: an empty list is given, FROSch building its own list
      // otherwise. With FROSch 17.2.1, the list of the rows whose entries are
      // null except the diagonal one (`FindOneEntryOnlyRowsGlobal` with a
      // tolerance, the matrices of MFEM keeping the eliminated entries) made
      // the initialization hang on several processes with three unknowns
      // per node, and it did not change the number of iterations on the
      // Laplace equation
      const auto repeated = FROSch::BuildRepeatedMap(xa->getCrsGraph());
      const auto no_dirichlet_rows =
          Teuchos::ArrayRCP<GO>(new GO[1], 0, 0, true);
      auto s = Teuchos::rcp(
          new FROSch::TwoLevelPreconditioner<Scalar, LO, GO, Node>(xa, p));
      // the return codes are local: the failure is decided on all the
      // processes before the collective calls
      const auto check = [&m](const int r) {
        auto failed = static_cast<int>(r < 0);
        MPI_Allreduce(MPI_IN_PLACE, &failed, 1, MPI_INT, MPI_LOR,
                      m.communicator);
        if (failed != 0) {
          throw std::runtime_error(
              "FROSch: the setup of the preconditioner failed");
        }
      };
      check(s->initialize(static_cast<unsigned>(d), static_cast<unsigned>(vdim),
                          p->get<int>("Overlap", 1), Teuchos::null, nodes,
                          FROSch::NodeWise, repeated, Teuchos::null,
                          no_dirichlet_rows));
      check(s->compute());
      this->preconditioner =
          Teuchos::rcp(new FROSch::TpetraPreconditioner<Scalar, LO, GO, Node>(
              Teuchos::RCP<Xpetra::Operator<Scalar, LO, GO, Node>>(s)));
    }
    [[nodiscard]] Teuchos::RCP<OP> getTpetraOperator() const override {
      return this->preconditioner;
    }

   private:
    //! \brief parameters
    const Teuchos::RCP<Teuchos::ParameterList> parameters;
    //! \brief coordinates of the nodes
    const NodesCoordinates coordinates;
    //! \brief preconditioner
    Teuchos::RCP<OP> preconditioner;
  };  // end of struct FROSchPreconditioner

  //! \brief direct solver of Amesos2, as an operator of Tpetra
  struct Amesos2Solver final : public TpetraOperatorBase {
    /*!
     * \param[in] s: solver
     * \param[in] o: options
     */
    Amesos2Solver(const std::string& s, const Options& o)
        : name(s), parameters(makeParameterList(o, "Amesos2")) {}
    void setMatrix(const DistributedMatrixView& m) override {
      this->A.set(m);
      this->solver = Amesos2::create<CrsMatrix, MV>(this->name, this->A.matrix);
      if (this->solver.is_null()) {
        throw std::runtime_error("Amesos2: the solver '" + this->name +
                                 "' is not available");
      }
      this->solver->setParameters(this->parameters);
      this->solver->symbolicFactorization();
      this->solver->numericFactorization();
    }
    [[nodiscard]] Teuchos::RCP<OP> getTpetraOperator() const override {
      return this->A.matrix;
    }
    void applyTpetraOperator(const MV& b, MV& x) const override {
      this->solver->solve(Teuchos::ptrFromRef(x), Teuchos::ptrFromRef(b));
    }

   private:
    //! \brief name of the solver
    const std::string name;
    //! \brief parameters
    const Teuchos::RCP<Teuchos::ParameterList> parameters;
    //! \brief solver
    Teuchos::RCP<Amesos2::Solver<CrsMatrix, MV>> solver;
  };  // end of struct Amesos2Solver

  /*!
   * \brief preconditioner given to Belos: a preconditioner of Trilinos,
   * applied to the vectors of Tpetra, or a preconditioner of mfem-mgis
   * applied to vectors of the host
   */
  struct PreconditionerOperator final : public OP {
    /*!
     * \param[in] m: matrix of the solver
     * \param[in] p: preconditioner of Trilinos (may be null)
     * \param[in] e: preconditioner of mfem-mgis (may be empty)
     */
    PreconditionerOperator(const TpetraMatrix& m,
                           const TpetraOperatorBase* const p,
                           ExternalPreconditioner e)
        : matrix(m), trilinos(p), external(std::move(e)) {}
    Teuchos::RCP<const Map> getDomainMap() const override {
      return this->matrix.row_map;
    }
    Teuchos::RCP<const Map> getRangeMap() const override {
      return this->matrix.row_map;
    }
    void apply(
        const MV& X,
        MV& Y,
        Teuchos::ETransp mode = Teuchos::NO_TRANS,
        Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
        Scalar beta = Teuchos::ScalarTraits<Scalar>::zero()) const override {
      if (mode != Teuchos::NO_TRANS) {
        throw std::runtime_error(
            "Belos: the transposed preconditioner is not available");
      }
      // exact comparisons: the result is the application of the
      // preconditioner only for these values
      const auto direct = (std::fpclassify(alpha - Scalar{1}) == FP_ZERO) &&
                          (std::fpclassify(beta) == FP_ZERO);
      if ((!direct) && ((this->r.is_null()) ||
                        (this->r->getNumVectors() != Y.getNumVectors()) ||
                        (this->r->getMap().get() != Y.getMap().get()))) {
        this->r = Teuchos::rcp(new MV(Y.getMap(), Y.getNumVectors()));
      }
      auto& result = direct ? Y : *(this->r);
      result.putScalar(Scalar{0});
      if (this->trilinos != nullptr) {
        this->trilinos->applyTpetraOperator(X, result);
      } else {
        // the vectors of the preconditioner of mfem-mgis are written and
        // read on the GPU when both compute there
        const auto device = this->external.device && computesOnTheGPU();
        for (std::size_t j = 0; j != X.getNumVectors(); ++j) {
          auto* const b = this->external.getRightHandSide();
          if (device) {
            copyOnDevice(b, X, j, this->matrix);
            Kokkos::fence("mfem_mgis::trilinos::PreconditionerOperator");
          } else {
            copy(b, X, j, this->matrix);
          }
          const auto* const x = this->external.apply();
          if (device) {
            // the result is read before mfem-mgis writes it again
            copyOnDevice(result, j, x, this->matrix);
            Kokkos::fence("mfem_mgis::trilinos::PreconditionerOperator");
          } else {
            copy(result, j, x, this->matrix);
          }
        }
      }
      if (!direct) {
        Y.update(alpha, result, beta);
      }
    }

   private:
    //! \brief matrix of the solver, giving the permutation of the unknowns
    const TpetraMatrix& matrix;
    //! \brief preconditioner of Trilinos
    const TpetraOperatorBase* const trilinos;
    //! \brief preconditioner of mfem-mgis
    const ExternalPreconditioner external;
    //! \brief result, when scaled and added to `Y`
    mutable Teuchos::RCP<MV> r;
  };  // end of struct PreconditionerOperator

  //! \brief Krylov solver of Belos
  struct BelosSolver final : public TpetraOperatorBase {
    /*!
     * \param[in] s: solver
     * \param[in] o: options
     * \param[in] p: preconditioner of Trilinos
     * \param[in] e: preconditioner of mfem-mgis
     */
    BelosSolver(const std::string& s,
                const Options& o,
                const TpetraOperatorBase* const p,
                ExternalPreconditioner e)
        : parameters(makeParameterList(o)),
          trilinos_preconditioner(p),
          external_preconditioner(std::move(e)) {
      auto factory = Belos::SolverFactory<Scalar, MV, OP>{};
      // Belos reports an unknown solver with its name in capitals
      if (!factory.isSupported(s)) {
        throw std::runtime_error("Belos: unknown solver '" + s + "'");
      }
      this->solver = factory.create(s, this->parameters);
    }
    void setMatrix(const DistributedMatrixView& m) override {
      // the matrix of a preconditioner of Trilinos, set before, is shared
      if (this->trilinos_preconditioner != nullptr) {
        this->A = this->trilinos_preconditioner->getMatrix();
        if (this->A.matrix.is_null()) {
          throw std::runtime_error(
              "Belos: the matrix of the preconditioner is not set");
        }
      } else {
        this->A.set(m);
      }
      this->problem = Teuchos::rcp(new Belos::LinearProblem<Scalar, MV, OP>());
      this->problem->setOperator(this->A.matrix);
      if ((this->trilinos_preconditioner != nullptr) ||
          (this->external_preconditioner)) {
        this->problem->setRightPrec(Teuchos::rcp(
            new PreconditionerOperator(this->A, this->trilinos_preconditioner,
                                       this->external_preconditioner)));
      }
    }
    [[nodiscard]] Teuchos::RCP<OP> getTpetraOperator() const override {
      return this->A.matrix;
    }
    void applyTpetraOperator(const MV& b, MV& x) const override {
      if (!this->problem->setProblem(Teuchos::rcpFromRef(x),
                                     Teuchos::rcpFromRef(b))) {
        throw std::runtime_error("Belos: invalid linear problem");
      }
      this->report = Report{};
      try {
        this->solver->setProblem(this->problem);
        const auto r = this->solver->solve();
        this->report.converged = r == Belos::Converged;
      } catch (std::exception& e) {
        this->report.failure = e.what();
      }
      // a failure of Belos (loss of positive definiteness, breakdown...) is
      // a resolution which did not converge, on all the processes
      auto failed = static_cast<int>(!this->report.failure.empty());
      MPI_Allreduce(MPI_IN_PLACE, &failed, 1, MPI_INT, MPI_LOR,
                    this->A.communicator);
      if (failed != 0) {
        this->report.converged = false;
        if (this->report.failure.empty()) {
          this->report.failure = "failure on another process";
        }
      }
      this->report.iterations = this->solver->getNumIters();
      try {
        // not implemented by all the solvers
        this->report.achieved_tolerance =
            static_cast<double>(this->solver->achievedTol());
      } catch (std::exception&) {
      }
    }
    [[nodiscard]] Report getReport() const override { return this->report; }
    [[nodiscard]] double getRelativeTolerance() const override {
      // the parameters of the solver, which may only contain the given
      // ones, then the default values
      const auto* const n = "Convergence Tolerance";
      for (const auto& p : {this->solver->getCurrentParameters(),
                            this->solver->getValidParameters()}) {
        if ((!p.is_null()) && (p->isType<Scalar>(n))) {
          return static_cast<double>(p->get<Scalar>(n));
        }
      }
      return -1;
    }
    [[nodiscard]] bool setRelativeTolerance(const double tol) override {
      this->parameters->set("Convergence Tolerance", static_cast<Scalar>(tol));
      this->solver->setParameters(this->parameters);
      return true;
    }

   private:
    //! \brief parameters
    const Teuchos::RCP<Teuchos::ParameterList> parameters;
    //! \brief preconditioner of Trilinos
    const TpetraOperatorBase* const trilinos_preconditioner;
    //! \brief preconditioner of mfem-mgis
    const ExternalPreconditioner external_preconditioner;
    //! \brief solver
    Teuchos::RCP<Belos::SolverManager<Scalar, MV, OP>> solver;
    //! \brief linear problem
    Teuchos::RCP<Belos::LinearProblem<Scalar, MV, OP>> problem;
    //! \brief report of the last resolution
    mutable Report report;
  };  // end of struct BelosSolver

  std::unique_ptr<Operator> makeIfpack2Preconditioner(const std::string& type,
                                                      const Options& options) {
    if (!Ifpack2::Factory::isSupported<RowMatrix>(type)) {
      auto msg = "invalid preconditioner '" + type + "' (supported:";
      for (const auto& n : Ifpack2::Factory::getSupportedNames<RowMatrix>()) {
        msg += " " + n;
      }
      throw std::runtime_error(msg + ")");
    }
    return std::make_unique<Ifpack2Preconditioner>(type, options);
  }  // end of makeIfpack2Preconditioner

  std::unique_ptr<Operator> makeMueLuPreconditioner(const Options& options,
                                                    NodesCoordinates c) {
    return std::make_unique<MueLuPreconditioner>(options, std::move(c));
  }  // end of makeMueLuPreconditioner

  std::unique_ptr<Operator> makeFROSchPreconditioner(const Options& options,
                                                     NodesCoordinates c) {
    return std::make_unique<FROSchPreconditioner>(options, std::move(c));
  }  // end of makeFROSchPreconditioner

  std::unique_ptr<Operator> makeAmesos2Solver(const std::string& solver,
                                              const Options& options) {
    return std::make_unique<Amesos2Solver>(solver, options);
  }  // end of makeAmesos2Solver

  bool isAmesos2SolverAvailable(const std::string& solver) {
    return Amesos2::query(solver);
  }  // end of isAmesos2SolverAvailable

  std::unique_ptr<Operator> makeBelosSolver(const std::string& solver,
                                            const Options& options,
                                            const Operator* const p,
                                            ExternalPreconditioner e) {
    const auto* const tp = dynamic_cast<const TpetraOperatorBase*>(p);
    if ((p != nullptr) && (tp == nullptr)) {
      throw std::runtime_error(
          "Belos: the preconditioner is not an operator of Tpetra");
    }
    return std::make_unique<BelosSolver>(solver, options, tp, std::move(e));
  }  // end of makeBelosSolver

}  // end of namespace mfem_mgis::trilinos
