/*!
 * \file   src/TrilinosTpetraBackend.hxx
 * \brief  This file implements the solvers and the preconditioners of
 * Trilinos for a node of Tpetra (`TpetraBackend`), explicitly instantiated
 * for the default node of Tpetra, on the GPU when Trilinos is built for
 * one, and for a node of the host (internal header, compiled by the
 * compiler of Kokkos, see `TrilinosBackend.hxx`)
 * \date   16/10/2026
 */

#ifndef LIB_MFEM_MGIS_SRC_TRILINOSTPETRABACKEND_HXX
#define LIB_MFEM_MGIS_SRC_TRILINOSTPETRABACKEND_HXX

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
#include "BelosSolverFactory_Generic.hpp"
#include "BelosLSQRSolMgr.hpp"
#include "BelosGCRODRSolMgr.hpp"
#include "Belos_Tpetra_SolverManagerBase.hpp"
#include "BelosTpetraAdapter.hpp"
#include "Ifpack2_Factory.hpp"
#include "Ifpack2_BlockRelaxation_decl.hpp"
#include "MueLu_CreateTpetraPreconditioner.hpp"
#include "Amesos2.hpp"
#include "Xpetra_CrsMatrixWrap.hpp"
#include "Xpetra_TpetraCrsMatrix.hpp"
#include "Xpetra_TpetraMultiVector.hpp"
#include "FROSch_TwoLevelPreconditioner_def.hpp"
#include "MueLu_MasterList.hpp"
#include "FROSch_TpetraPreconditioner_def.hpp"
#include "Tpetra_Details_DefaultTypes.hpp"
#include "KokkosSparse_SortCrs.hpp"
#include "TrilinosBackend.hxx"

// Tpetra computing on a GPU by default, and instantiated for the serial
// node: the solvers compute on the host when MFEM does. The OpenMP node is
// not used: with Kokkos built for a GPU, some products of sparse matrices of
// Tpetra 17.2.1 run in the default execution space of Kokkos on the data of
// this node (TpetraExt_MatrixMatrix_def.hpp)
#if (defined(HAVE_TPETRA_DEFAULTNODE_CUDAWRAPPERNODE) ||  \
     defined(HAVE_TPETRA_DEFAULTNODE_HIPWRAPPERNODE) ||   \
     defined(HAVE_TPETRA_DEFAULTNODE_SYCLWRAPPERNODE)) && \
    defined(HAVE_TPETRA_INST_SERIAL)
#define MFEM_MGIS_TRILINOS_HOST_NODE
#endif

namespace mfem_mgis::trilinos {

  //! \brief scalar type
  using Scalar = double;
  //! \brief local ordinal (default one of Tpetra)
  using LO = Tpetra::Map<>::local_ordinal_type;
  //! \brief global ordinal (default one of Tpetra)
  using GO = Tpetra::Map<>::global_ordinal_type;
  //! \brief default node of Tpetra, on the GPU when Trilinos is built for one
  using DefaultNode = Tpetra::Map<>::node_type;
#ifdef MFEM_MGIS_TRILINOS_HOST_NODE
  /*!
   * \brief node of Tpetra on the host: serial, each process of MPI
   * computing on one core
   */
  using HostNode = Tpetra::KokkosCompat::KokkosSerialWrapperNode;
#else  /* MFEM_MGIS_TRILINOS_HOST_NODE */
  /*!
   * \brief node of Tpetra on the host: the default one, the only one used
   * (on the host, or on the GPU when Tpetra is not instantiated for the
   * serial node)
   */
  using HostNode = DefaultNode;
#endif /* MFEM_MGIS_TRILINOS_HOST_NODE */

  /*!
   * \brief check that Kokkos is initialized (by mfem-mgis, see
   * `initializeKokkosOnDemand`, or by the program)
   */
  inline void checkKokkosInitialization() {
    if ((!Kokkos::is_initialized()) || (Kokkos::is_finalized())) {
      throw std::runtime_error("Trilinos: Kokkos is not initialized");
    }
  }  // end of checkKokkosInitialization

  /*!
   * \brief fill a parameter list of Teuchos with options
   * \param[out] l: list
   * \param[in] options: options
   */
  inline void fill(Teuchos::ParameterList& l, const Options& options) {
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
  inline Teuchos::RCP<Teuchos::ParameterList> makeParameterList(
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
  inline long long interleave(const long long i,
                              const long long n,
                              const long long c) {
    const auto nn = n / c;
    return (i % nn) * c + i / nn;
  }  // end of interleave

  /*!
   * \brief checks that the names of the first level of a list of parameters
   * are known, so that a misspelt parameter is refused at the creation of
   * the solver (Ifpack2 and MueLu check them at the setup, or ignore them)
   * \param[in] library: name of the library
   * \param[in] p: parameters
   * \param[in] valid: known parameters
   * \param[in] prefixes: prefixes of the other names accepted
   */
  inline void checkParameterNames(const std::string& library,
                                  const Teuchos::ParameterList& p,
                                  const Teuchos::ParameterList& valid,
                                  const std::vector<std::string>& prefixes) {
    for (auto i = p.begin(); i != p.end(); ++i) {
      const auto& n = p.name(i);
      const auto prefixed = std::ranges::any_of(
          prefixes, [&n](const std::string& prefix) noexcept {
            return n.starts_with(prefix);
          });
      if ((!valid.isParameter(n)) && (!prefixed)) {
        throw std::runtime_error(library + ": unknown parameter '" + n + "'");
      }
    }
  }  // end of checkParameterNames

  /*!
   * \brief makes the setup of FROSch less verbose (it prints tables at each
   * setup by default), unless a verbosity is given: the verbosity is read in
   * the parameters of the overlapping operator and of the first block of the
   * coarse operator. FROSch still prints some tables on the first process.
   * \param[in, out] p: parameters of FROSch
   */
  inline void setDefaultFROSchVerbosity(Teuchos::ParameterList& p) {
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

  /*!
   * \brief solvers and preconditioners of Trilinos computing on a node of
   * Tpetra
   * \tparam Node: node of Tpetra
   */
  template <typename Node>
  struct TpetraBackend {
    //! \brief map of Tpetra
    using Map = Tpetra::Map<LO, GO, Node>;
    //! \brief sparse matrix of Tpetra
    using CrsMatrix = Tpetra::CrsMatrix<Scalar, LO, GO, Node>;
    //! \brief row matrix of Tpetra
    using RowMatrix = Tpetra::RowMatrix<Scalar, LO, GO, Node>;
    //! \brief multivector of Tpetra
    using MV = Tpetra::MultiVector<Scalar, LO, GO, Node>;
    //! \brief operator of Tpetra
    using OP = Tpetra::Operator<Scalar, LO, GO, Node>;
    //! \brief execution space of Tpetra
    using ExecutionSpace = typename Node::execution_space;
    //! \brief memory space of Tpetra
    using MemorySpace = typename Node::memory_space;
    //! \brief view of an array of the memory space of Tpetra
    template <typename T>
    using DeviceArray = Kokkos::View<T*, MemorySpace>;
    //! \brief view of an array of the memory space of Tpetra, not owned
    template <typename T>
    using UnmanagedDeviceArray =
        Kokkos::View<T*, MemorySpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    /*!
     * \brief true if the node computes on a GPU (its execution space can't
     * access the memory of the host)
     */
    static constexpr bool gpu =
        !Kokkos::SpaceAccessibility<ExecutionSpace,
                                    Kokkos::HostSpace>::accessible;

    /*!
     * \return an array of the memory space of Tpetra: the given one, if it is
     * already there, or a copy
     * \param[in] v: array
     * \param[in] s: size
     * \param[in] device: if true, the array is in the memory of the GPU
     * which Tpetra computes on, and in the memory of the host otherwise
     */
    template <typename T>
    static DeviceArray<const T> toMemorySpace(const T* const v,
                                              const std::size_t s,
                                              const bool device) {
      if ((device) || (!gpu)) {
        return UnmanagedDeviceArray<const T>(v, s);
      }
      auto c = DeviceArray<T>(
          Kokkos::view_alloc(Kokkos::WithoutInitializing, "copy"), s);
      Kokkos::deep_copy(
          c, Kokkos::View<const T*, Kokkos::HostSpace,
                          Kokkos::MemoryTraits<Kokkos::Unmanaged>>(v, s));
      return c;
    }  // end of toMemorySpace
    /*!
     * \return the last value of an array of offsets, zero if it is empty
     * \param[in] o: offsets
     */
    static std::size_t getLastValue(const DeviceArray<const int>& o) {
      if (o.extent(0) == 0) {
        return 0;
      }
      auto v = int{};
      Kokkos::deep_copy(v, Kokkos::subview(o, o.extent(0) - 1));
      return static_cast<std::size_t>(v);
    }  // end of getLastValue
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
        const auto comm =
            Teuchos::rcp(new Teuchos::MpiComm<int>(m.communicator));
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
          columns[static_cast<std::size_t>(i)] =
              static_cast<GO>(m.first_row + i);
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
            c = rows[2 * p] +
                interleave(g - rows[2 * p], rows[2 * p + 1], vdim);
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
        // local matrix, built in the memory space of Tpetra, the columns of
        // each row being sorted
        using LocalMatrix = typename CrsMatrix::local_matrix_device_type;
        using Offsets =
            typename CrsMatrix::local_graph_device_type::row_map_type;
        using Entries =
            typename CrsMatrix::local_graph_device_type::entries_type;
        using Values = typename LocalMatrix::values_type;
        const auto un = static_cast<std::size_t>(n);
        const auto has_offd = (no != 0) && (m.offd_offsets != nullptr);
        // the offsets of a process without rows may not be allocated
        const auto no_offsets = std::size_t{n == 0 ? 0u : un + 1};
        const auto diag_offsets =
            toMemorySpace(m.diag_offsets, no_offsets, m.device);
        const auto diag_nnz = getLastValue(diag_offsets);
        const auto diag_columns =
            toMemorySpace(m.diag_columns, diag_nnz, m.device);
        const auto diag_values =
            toMemorySpace(m.diag_values, diag_nnz, m.device);
        auto offd_offsets = DeviceArray<const int>{};
        auto offd_columns = DeviceArray<const int>{};
        auto offd_values = DeviceArray<const double>{};
        auto offd_nnz = std::size_t{};
        if (has_offd) {
          offd_offsets = toMemorySpace(m.offd_offsets, no_offsets, m.device);
          offd_nnz = getLastValue(offd_offsets);
          offd_columns = toMemorySpace(m.offd_columns, offd_nnz, m.device);
          offd_values = toMemorySpace(m.offd_values, offd_nnz, m.device);
        }
        const auto nnz = diag_nnz + offd_nnz;
        auto offsets = typename Offsets::non_const_type(
            Kokkos::view_alloc(Kokkos::WithoutInitializing, "offsets"), un + 1);
        auto entries = typename Entries::non_const_type(
            Kokkos::view_alloc(Kokkos::WithoutInitializing, "entries"), nnz);
        auto values = typename Values::non_const_type(
            Kokkos::view_alloc(Kokkos::WithoutInitializing, "values"), nnz);
        // rows of the process in the order of Tpetra
        const auto p = this->device_permutation;
        auto source = DeviceArray<int>(
            Kokkos::view_alloc(Kokkos::WithoutInitializing, "source"), un);
        Kokkos::parallel_for(
            "mfem_mgis::trilinos::TpetraMatrix::source",
            Kokkos::RangePolicy<ExecutionSpace>(0, n),
            KOKKOS_LAMBDA(const int i) { source(permuted ? p(i) : i) = i; });
        // offsets of the rows
        Kokkos::parallel_scan(
            "mfem_mgis::trilinos::TpetraMatrix::offsets",
            Kokkos::RangePolicy<ExecutionSpace>(0, n + 1),
            KOKKOS_LAMBDA(const int r, std::size_t& o, const bool last) {
              if (last) {
                offsets(r) = o;
              }
              if (r != n) {
                const auto i = source(r);
                o += static_cast<std::size_t>(diag_offsets(i + 1) -
                                              diag_offsets(i));
                if (has_offd) {
                  o += static_cast<std::size_t>(offd_offsets(i + 1) -
                                                offd_offsets(i));
                }
              }
            });
        // entries of the rows, then sorted by Kokkos Kernels
        Kokkos::parallel_for(
            "mfem_mgis::trilinos::TpetraMatrix::entries",
            Kokkos::RangePolicy<ExecutionSpace>(0, n),
            KOKKOS_LAMBDA(const int r) {
              const auto i = source(r);
              auto k = offsets(r);
              for (auto j = diag_offsets(i); j != diag_offsets(i + 1); ++j) {
                const auto c = diag_columns(j);
                entries(k) = static_cast<LO>(permuted ? p(c) : c);
                values(k) = static_cast<Scalar>(diag_values(j));
                ++k;
              }
              if (has_offd) {
                for (auto j = offd_offsets(i); j != offd_offsets(i + 1); ++j) {
                  entries(k) = static_cast<LO>(n + offd_columns(j));
                  values(k) = static_cast<Scalar>(offd_values(j));
                  ++k;
                }
              }
            });
        // on the GPU, each row is sorted by one thread (and its vector
        // lanes), without the temporary arrays of the sort of all the
        // entries, which Kokkos Kernels selects for long or imbalanced rows
        // (about 32 bytes per entry)
        auto algorithm = KokkosSparse::SortAlgorithm::DEFAULT;
        if constexpr (gpu) {
          algorithm = KokkosSparse::SortAlgorithm::PARALLEL_THREAD_LEVEL;
        }
        KokkosSparse::sort_crs_matrix(ExecutionSpace{}, offsets, entries,
                                      values, static_cast<LO>(n + no),
                                      algorithm);
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
          v(static_cast<std::size_t>(m.getIndex(static_cast<int>(i))), j) =
              s[i];
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
      [[nodiscard]] bool computesOnTheGPU() const noexcept override {
        return gpu;
      }
      void apply(const double* b, double* x) const override {
        this->allocateVectors();
        copy(*(this->b), 0, b, this->A);
        copy(*(this->x), 0, x, this->A);
        this->applyTpetraOperator(*(this->b), *(this->x));
        copy(x, *(this->x), 0, this->A);
      }
      void applyOnDevice(const double* b, double* x) const override {
        if (!gpu) {
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
        const auto op = this->getTpetraOperator();
        if (op.is_null()) {
          // for instance, a preconditioner whose setup failed
          throw std::runtime_error("Trilinos: the operator is not set up");
        }
        op->apply(b, x);
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

    //! \brief preconditioner of Ifpack2
    struct Ifpack2Preconditioner final : public TpetraOperatorBase {
      /*!
       * \param[in] t: type
       * \param[in] o: options
       */
      Ifpack2Preconditioner(const std::string& t, const Options& o)
          : type(t), parameters(makeParameterList(o)) {
        // the parameters are given to a preconditioner of this type built
        // without matrix, so that the ones checked by Ifpack2 (the names of
        // the relaxation) are refused at the creation of the solver.
        // `Ifpack2::getValidParameters` does not list all the parameters of
        // the preconditioners. The block relaxations, which check their
        // parameters at the setup, read the matrix in `setParameters` (a
        // negative number of blocks gives the number of rows per block).
        auto p = Teuchos::RCP<Ifpack2::Preconditioner<Scalar, LO, GO, Node>>{};
        try {
          p = Ifpack2::Factory::create<RowMatrix>(
              this->type, Teuchos::RCP<const RowMatrix>{});
        } catch (std::exception&) {
          // a preconditioner which needs its matrix: the parameters are
          // checked at the setup
        }
        using BlockRelaxation = Ifpack2::BlockRelaxation<RowMatrix>;
        if ((!p.is_null()) &&
            (dynamic_cast<const BlockRelaxation*>(p.get()) == nullptr)) {
          auto copy = *(this->parameters);
          try {
            p->setParameters(copy);
          } catch (std::exception& e) {
            throw std::runtime_error("Ifpack2: " + std::string{e.what()});
          }
        }
      }
      void setMatrix(const DistributedMatrixView& m) override {
        this->A.set(m);
        auto p = Ifpack2::Factory::create<RowMatrix>(
            this->type, Teuchos::RCP<const RowMatrix>(this->A.matrix));
        if (p.is_null()) {
          throw std::runtime_error("invalid preconditioner '" + this->type +
                                   "'");
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
      Teuchos::RCP<Ifpack2::Preconditioner<Scalar, LO, GO, Node>>
          preconditioner;
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
      const auto nn =
          map.getLocalNumElements() / static_cast<std::size_t>(vdim);
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
          m.row_map->getGlobalNumElements() / static_cast<std::size_t>(vdim),
          nn, GO{0}, m.row_map->getComm()));
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
    static Teuchos::RCP<MV> makeRigidBodyModes(
        const MPI_Comm comm,
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
          : parameters(makeParameterList(o)),
            simplified(isSimplifiedSyntax(*(this->parameters))),
            coordinates(std::move(c)),
            reuse(this->parameters->isType<std::string>("reuse: type")
                      ? this->parameters->get<std::string>("reuse: type")
                      : std::string{"none"}) {
        // the names of the simplified syntax of MueLu (the default one), which
        // MueLu checks at the setup, are checked at the creation of the
        // solver. The parameters of a level are given in a sublist `level N`.
        // The other syntaxes are checked by MueLu only.
        if (this->simplified) {
          checkParameterNames("MueLu", *(this->parameters),
                              *(MueLu::MasterList::List()), {"level "});
        }
      }
      void setMatrix(const DistributedMatrixView& m) override {
        // the hierarchy of the previous operator is updated when the
        // parameter `reuse: type` of MueLu asks for it (`RP`: the
        // prolongators and the restrictions are kept, the coarse operators
        // computed again; `tP`: the tentative prolongators are kept...) and
        // when the distribution of the rows is unchanged
        const auto previous_rows = this->A.row_map;
        this->A.set(m);
        if ((!this->preconditioner.is_null()) && (this->reuse != "none") &&
            (!previous_rows.is_null()) &&
            (previous_rows->isSameAs(*(this->A.row_map)))) {
          try {
            MueLu::ReuseTpetraPreconditioner(this->A.matrix,
                                             *(this->preconditioner));
          } catch (...) {
            // the hierarchy, partially updated, is built again for the next
            // operator
            this->preconditioner = Teuchos::null;
            throw;
          }
          return;
        }
        auto p = *(this->parameters);
        const auto vdim = m.number_of_components;
        // number of equations per node, given in the simplified syntax, the
        // number of components otherwise (the other syntaxes ignore this
        // parameter, but the user data, given in all the syntaxes)
        const auto* const equations = "number of equations";
        if ((this->simplified) && (vdim > 1) && (!p.isParameter(equations))) {
          p.set(equations, vdim);
        }
        const auto neq =
            p.isType<int>(equations) ? p.get<int>(equations) : vdim;
        // the rigid body modes, unless the parameters give another number of
        // equations per node
        if (neq == vdim) {
          const auto modes = makeRigidBodyModes(m.communicator, this->A.row_map,
                                                vdim, this->coordinates);
          if (!modes.is_null()) {
            p.sublist("user data").set("Nullspace", modes);
          }
        }
        if (this->simplified) {
          setDefaultParameters(p, m.communicator);
        }
        // coordinates of the nodes, used by some options of MueLu (dropping
        // based on the distances, repartitioning), when the nodes of MueLu are
        // the nodes of the mesh
        if (neq == vdim) {
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
      /*!
       * \return if the parameters use the simplified syntax of MueLu (the
       * default one), and not the factories (sublist `Hierarchy`), the
       * facade (`MueLu preconditioner`) or the parameters of ML
       * \param[in] p: parameters
       */
      static bool isSimplifiedSyntax(const Teuchos::ParameterList& p) {
        const auto* const syntax = "parameterlist: syntax";
        return (!p.isSublist("Hierarchy")) &&
               (!p.isParameter("MueLu preconditioner")) &&
               (!((p.isType<std::string>(syntax)) &&
                  (p.get<std::string>(syntax) == "ml")));
      }  // end of isSimplifiedSyntax
#ifdef KOKKOS_ENABLE_CUDA
      /*!
       * \return true if a GPU is used by several processes of a node, on
       * all the processes: the UUIDs of the GPUs used by the processes of
       * each node are compared (the GPUs visible to a process may be
       * restricted, for instance by `CUDA_VISIBLE_DEVICES`; the instances of
       * a GPU partitioned by MIG have their own memories, and their own
       * UUIDs with the recent drivers, otherwise they are seen as one shared
       * GPU, which only keeps the algorithms using less memory)
       * \param[in] comm: communicator
       */
      static bool isGPUShared(const MPI_Comm comm) {
        auto node = MPI_Comm{};
        MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,
                            &node);
        auto processes = 0;
        MPI_Comm_size(node, &processes);
        const auto cuda = Kokkos::Cuda{};
        const auto& properties = cuda.cuda_device_prop();
        constexpr auto n = sizeof(properties.uuid.bytes);
        auto local = std::array<char, n>{};
        std::ranges::copy(properties.uuid.bytes, local.begin());
        auto uuids = std::vector<std::array<char, n>>(
            static_cast<std::size_t>(processes));
        MPI_Allgather(local.data(), static_cast<int>(n), MPI_BYTE, uuids.data(),
                      static_cast<int>(n), MPI_BYTE, node);
        MPI_Comm_free(&node);
        std::ranges::sort(uuids);
        auto shared =
            static_cast<int>(std::ranges::adjacent_find(uuids) != uuids.end());
        MPI_Allreduce(MPI_IN_PLACE, &shared, 1, MPI_INT, MPI_LOR, comm);
        return shared != 0;
      }  // end of isGPUShared
#endif /* KOKKOS_ENABLE_CUDA */
      /*!
       * \brief set the default parameters of mfem-mgis, unless given, in the
       * simplified syntax (the other syntaxes ignore them)
       * \param[in, out] p: parameters
       */
      static void setDefaultParameters(Teuchos::ParameterList& p,
                                       const MPI_Comm comm) {
#ifdef KOKKOS_ENABLE_CUDA
        // products of sparse matrices by the algorithm of Kokkos Kernels on
        // CUDA: the default one (cuSPARSE) ran out of memory on a GPU of 16 GB
        // from 421,443 unknowns of linear elasticity. Kokkos Kernels reserves
        // up to half the free memory of the GPU for each product
        // (KokkosSparse_spgemm_impl.hpp): several processes sharing a GPU ran
        // out of memory, and keep the default algorithm. Tpetra also reads
        // this parameter on the host.
        if constexpr (std::is_same_v<typename Node::execution_space,
                                     Kokkos::Cuda>) {
          auto& products = p.sublist("matrixmatrix: kernel params");
          if ((!products.isParameter("cuda: algorithm")) &&
              (!isGPUShared(comm))) {
            products.set("cuda: algorithm", "SPGEMM_KK");
          }
        }
#else  /* KOKKOS_ENABLE_CUDA */
        static_cast<void>(comm);
#endif /* KOKKOS_ENABLE_CUDA */
        // the factories of MueLu which are not written with Kokkos (dropping,
        // aggregation and tentative prolongator computed on the host), MueLu
        // using its factories written with Kokkos on all the nodes but the
        // serial one by default: up to Trilinos 17.2.1, the local QR
        // factorization of the tentative prolongator written with Kokkos
        // does not compute the Q factor of a nullspace of several vectors
        // (replaced in the development version of Trilinos after 17.2.1): on
        // the linear elasticity with the rigid body modes, twice the
        // iterations (22 against 12 at 421,443 unknowns on the GPU), and a
        // longer setup (0.81 against 0.51 seconds)
        if (!p.isParameter("use kokkos refactor")) {
          p.set("use kokkos refactor", false);
        }
        // reproducible aggregates with the factories written with Kokkos,
        // which aggregate the nodes in a random order by default (the other
        // ones follow `aggregation: ordering`, natural by default), and
        // silent setup
        if (!p.isParameter("aggregation: deterministic")) {
          p.set("aggregation: deterministic", true);
        }
        if (!p.isParameter("verbosity")) {
          p.set("verbosity", "none");
        }
      }  // end of setDefaultParameters
      //! \brief parameters
      const Teuchos::RCP<Teuchos::ParameterList> parameters;
      //! \brief if true, the parameters use the simplified syntax of MueLu
      const bool simplified;
      //! \brief coordinates of the nodes
      const NodesCoordinates coordinates;
      //! \brief kind of reuse of the hierarchy (parameter `reuse: type`)
      const std::string reuse;
      //! \brief preconditioner
      Teuchos::RCP<MueLu::TpetraOperator<Scalar, LO, GO, Node>> preconditioner;
    };  // end of struct MueLuPreconditioner

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
        using XpetraMultiVector =
            Xpetra::TpetraMultiVector<Scalar, LO, GO, Node>;
        // FROSch 17.2.1 crashes when a process owns no unknown (constant
        // partition of unity): checked by all the processes, which fail
        // together
        auto empty = static_cast<int>(m.local_size == 0);
        MPI_Allreduce(MPI_IN_PLACE, &empty, 1, MPI_INT, MPI_LOR,
                      m.communicator);
        if (empty != 0) {
          throw std::runtime_error(
              "FROSch: each process must own unknowns (FROSch does not support "
              "empty subdomains)");
        }
        this->A.set(m);
        const auto vdim = m.number_of_components;
        const auto d = this->coordinates.dimension;
        if ((d != 2) && (d != 3)) {
          throw std::runtime_error(
              "FROSch: the space dimension must be 2 or 3");
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
        const auto xa =
            Teuchos::RCP<const Xpetra::Matrix<Scalar, LO, GO, Node>>(
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
        check(s->initialize(
            static_cast<unsigned>(d), static_cast<unsigned>(vdim),
            p->get<int>("Overlap", 1), Teuchos::null, nodes, FROSch::NodeWise,
            repeated, Teuchos::null, no_dirichlet_rows));
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
        this->solver =
            Amesos2::create<CrsMatrix, MV>(this->name, this->A.matrix);
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
          const auto device = this->external.device && gpu;
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
       * \param[in] settings: settings of mfem-mgis
       * \param[in] p: preconditioner of Trilinos
       * \param[in] e: preconditioner of mfem-mgis
       */
      BelosSolver(const std::string& s,
                  const Options& o,
                  const BelosSettings& settings,
                  const TpetraOperatorBase* const p,
                  ExternalPreconditioner e)
          : name(s),
            parameters(makeParameterList(o)),
            trilinos_preconditioner(p),
            external_preconditioner(std::move(e)) {
        auto factory = Belos::SolverFactory<Scalar, MV, OP>{};
        // Belos reports an unknown solver with its name in capitals
        if (!factory.isSupported(s)) {
          throw std::runtime_error("Belos: unknown solver '" + s + "'");
        }
        try {
          // a solver built without parameters gives the parameters it
          // accepts (some solvers, such as MINRES, write their default
          // values in the parameters they receive)
          const auto probe = factory.create(s, Teuchos::parameterList());
          const auto valid = probe->getValidParameters();
          // the hybrid GMRES of Belos (GMRES polynomial used as
          // preconditioner) iterates with an outer solver, which receives
          // the convergence parameters
          this->hybrid = (!valid->isParameter("Convergence Tolerance")) &&
                         (valid->isParameter("Outer Solver Params"));
          if (this->hybrid) {
            this->setOuterSolver(settings);
            // the verbosity of the polynomial
            auto verbosity = BelosSettings{};
            verbosity.verbosity_level = settings.verbosity_level;
            addSettings(*(this->parameters), *valid, verbosity);
          } else {
            addSettings(*(this->parameters), *valid, settings);
          }
          // the solvers of Belos only check the parameters from the second
          // call to `setParameters` (a change of the tolerance): they are
          // checked at the creation, except for the solvers which ignore
          // the unknown parameters (LSQR, the native solvers of Tpetra)
          const auto* const lsqr =
              dynamic_cast<const Belos::LSQRSolMgr<Scalar, MV, OP>*>(
                  probe.get());
          const auto* const native = dynamic_cast<
              const BelosTpetra::Impl::SolverManagerBase<Scalar, MV, OP>*>(
              probe.get());
          if ((lsqr == nullptr) && (native == nullptr)) {
            this->parameters->validateParameters(*valid, 0);
          }
          this->solver = factory.create(s, this->parameters);
        } catch (std::exception& ex) {
          throw std::runtime_error("Belos: " + std::string{ex.what()});
        }
      }
      void setMatrix(const DistributedMatrixView& m) override {
        const auto previous_rows = this->A.row_map;
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
        this->problem =
            Teuchos::rcp(new Belos::LinearProblem<Scalar, MV, OP>());
        this->problem->setOperator(this->A.matrix);
        if ((this->trilinos_preconditioner != nullptr) ||
            (this->external_preconditioner)) {
          this->problem->setRightPrec(Teuchos::rcp(
              new PreconditionerOperator(this->A, this->trilinos_preconditioner,
                                         this->external_preconditioner)));
        }
        // the solvers of Belos keep data of the previous resolutions: the
        // recycled subspace of GCRODR is kept for a matrix of the same rows
        // (the purpose of GCRODR), the other solvers are built again (the
        // seed space of PCPG is only valid for the same matrix, and the
        // polynomial of the hybrid GMRES is built for the matrix; GCRODR and
        // the native solvers of Tpetra throw on a reset before the vectors
        // are given)
        using GCRODR = Belos::GCRODRSolMgr<Scalar, MV, OP>;
        const auto recycling =
            (dynamic_cast<const GCRODR*>(this->solver.get()) != nullptr) &&
            (!previous_rows.is_null()) &&
            (previous_rows->isSameAs(*(this->A.row_map)));
        if (!recycling) {
          auto factory = Belos::SolverFactory<Scalar, MV, OP>{};
          this->solver = factory.create(this->name, this->parameters);
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
        const auto* const n = "Convergence Tolerance";
        if (this->hybrid) {
          // the tolerance of the outer solver, given or by default
          if (this->outer_valid_parameters.is_null()) {
            return -1;
          }
          const Teuchos::ParameterList& p =
              this->parameters->sublist("Outer Solver Params");
          const Teuchos::ParameterList& valid = *(this->outer_valid_parameters);
          for (const auto* const l : {&p, &valid}) {
            if (l->template isType<Scalar>(n)) {
              return static_cast<double>(l->template get<Scalar>(n));
            }
          }
          return -1;
        }
        // the parameters of the solver, which may only contain the given
        // ones, then the default values
        for (const auto& p : {this->solver->getCurrentParameters(),
                              this->solver->getValidParameters()}) {
          if ((!p.is_null()) && (p->template isType<Scalar>(n))) {
            return static_cast<double>(p->template get<Scalar>(n));
          }
        }
        return -1;
      }
      [[nodiscard]] bool setRelativeTolerance(const double tol) override {
        const auto* const n = "Convergence Tolerance";
        if (this->hybrid) {
          // the tolerance of the outer solver of the hybrid GMRES, if any
          if (this->outer_valid_parameters.is_null()) {
            return false;
          }
          this->parameters->sublist("Outer Solver Params")
              .set(n, static_cast<Scalar>(tol));
        } else {
          if (!this->solver->getValidParameters()->isParameter(n)) {
            // LSQR has its own tolerances
            return false;
          }
          this->parameters->set(n, static_cast<Scalar>(tol));
        }
        this->solver->setParameters(this->parameters);
        return true;
      }

     private:
      /*!
       * \brief add the settings of mfem-mgis to the parameters of a solver
       * of Belos, when the solver accepts them and when they are not given
       * \param[in, out] p: parameters
       * \param[in] valid: parameters accepted by the solver
       * \param[in] s: settings
       */
      static void addSettings(Teuchos::ParameterList& p,
                              const Teuchos::ParameterList& valid,
                              const BelosSettings& s) {
        const auto add = [&p, &valid](const char* const n, const auto v) {
          if ((valid.isParameter(n)) && (!p.isParameter(n))) {
            p.set(n, v);
          }
        };
        if (s.relative_tolerance.has_value()) {
          if (!valid.isParameter("Convergence Tolerance")) {
            throw std::runtime_error(
                "the solver has no relative tolerance (parameter "
                "'Convergence Tolerance')");
          }
          add("Convergence Tolerance",
              static_cast<Scalar>(*(s.relative_tolerance)));
        }
        if (s.maximum_number_of_iterations.has_value()) {
          add("Maximum Iterations", *(s.maximum_number_of_iterations));
        }
        // restarted GMRES: 50 Krylov vectors (300 by default in Belos), and
        // the number of restarts allowing the maximum number of iterations
        if ((valid.isParameter("Num Blocks")) &&
            (valid.isParameter("Maximum Restarts"))) {
          add("Num Blocks", 50);
          const auto* const mi = "Maximum Iterations";
          const auto nb = std::max(p.get<int>("Num Blocks"), 1);
          const auto maxit = p.isParameter(mi)       ? p.get<int>(mi)
                             : valid.isParameter(mi) ? valid.get<int>(mi)
                                                     : 1000;
          add("Maximum Restarts", (maxit + nb - 1) / nb);
        }
        if (s.verbosity_level.has_value()) {
          // Belos::Warnings and Belos::FinalSummary, then
          // Belos::IterationDetails and Belos::StatusTestDetails at each
          // iteration
          const auto l = *(s.verbosity_level);
          add("Verbosity", l > 1 ? 1 + 2 + 8 + 32 : (l > 0 ? 1 + 8 : 0));
          if (l > 1) {
            add("Output Frequency", 1);
          }
        }
      }  // end of addSettings
      /*!
       * \brief set the outer solver of the hybrid GMRES (GMRES by default,
       * the polynomial being applied alone if an empty name is given), and
       * add the settings of mfem-mgis to its parameters
       * \param[in] settings: settings of mfem-mgis
       */
      void setOuterSolver(const BelosSettings& settings) {
        const auto* const n = "Outer Solver";
        if (!this->parameters->isParameter(n)) {
          this->parameters->set(n, std::string{"Block GMRES"});
        }
        const auto outer = this->parameters->get<std::string>(n);
        if (outer.empty()) {
          if ((settings.relative_tolerance.has_value()) ||
              (settings.maximum_number_of_iterations.has_value())) {
            throw std::runtime_error(
                "the hybrid GMRES without outer solver has no relative "
                "tolerance nor maximum number of iterations");
          }
          return;
        }
        // the outer solver is built by the hybrid GMRES with the generic
        // factory of Belos
        auto factory =
            Belos::GenericSolverFactory<Scalar, Belos::MultiVec<Scalar>,
                                        Belos::Operator<Scalar>>{};
        if (!factory.isSupported(outer)) {
          throw std::runtime_error("unknown outer solver '" + outer + "'");
        }
        this->outer_valid_parameters =
            factory.create(outer, Teuchos::parameterList())
                ->getValidParameters();
        auto& p = this->parameters->sublist("Outer Solver Params");
        addSettings(p, *(this->outer_valid_parameters), settings);
        p.validateParameters(*(this->outer_valid_parameters), 0);
      }  // end of setOuterSolver
      //! \brief name of the solver
      const std::string name;
      //! \brief if true, the solver is the hybrid GMRES of Belos
      bool hybrid = false;
      //! \brief parameters accepted by the outer solver of the hybrid GMRES
      Teuchos::RCP<const Teuchos::ParameterList> outer_valid_parameters;
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

    //! \return a preconditioner of Ifpack2 (see `makeIfpack2Preconditioner`)
    static std::unique_ptr<Operator> makeIfpack2Preconditioner(
        const std::string&, const Options&);
    //! \return the multigrid of MueLu (see `makeMueLuPreconditioner`)
    static std::unique_ptr<Operator> makeMueLuPreconditioner(const Options&,
                                                             NodesCoordinates);
    //! \return the preconditioner of FROSch (see `makeFROSchPreconditioner`)
    static std::unique_ptr<Operator> makeFROSchPreconditioner(const Options&,
                                                              NodesCoordinates);
    //! \return a direct solver of Amesos2 (see `makeAmesos2Solver`)
    static std::unique_ptr<Operator> makeAmesos2Solver(const std::string&,
                                                       const Options&);
    //! \return a Krylov solver of Belos (see `makeBelosSolver`)
    static std::unique_ptr<Operator> makeBelosSolver(const std::string&,
                                                     const Options&,
                                                     const BelosSettings&,
                                                     const Operator*,
                                                     ExternalPreconditioner);
  };  // end of struct TpetraBackend

  template <typename Node>
  std::unique_ptr<Operator> TpetraBackend<Node>::makeIfpack2Preconditioner(
      const std::string& type, const Options& options) {
    if (!Ifpack2::Factory::isSupported<RowMatrix>(type)) {
      // some names contain spaces
      auto msg = "invalid preconditioner '" + type + "' (supported: ";
      auto first = true;
      for (const auto& n : Ifpack2::Factory::getSupportedNames<RowMatrix>()) {
        msg += (first ? "'" : ", '") + n + "'";
        first = false;
      }
      throw std::runtime_error(msg + ")");
    }
    return std::make_unique<Ifpack2Preconditioner>(type, options);
  }  // end of makeIfpack2Preconditioner

  template <typename Node>
  std::unique_ptr<Operator> TpetraBackend<Node>::makeMueLuPreconditioner(
      const Options& options, NodesCoordinates c) {
    return std::make_unique<MueLuPreconditioner>(options, std::move(c));
  }  // end of makeMueLuPreconditioner

  template <typename Node>
  std::unique_ptr<Operator> TpetraBackend<Node>::makeFROSchPreconditioner(
      const Options& options, NodesCoordinates c) {
    return std::make_unique<FROSchPreconditioner>(options, std::move(c));
  }  // end of makeFROSchPreconditioner

  template <typename Node>
  std::unique_ptr<Operator> TpetraBackend<Node>::makeAmesos2Solver(
      const std::string& solver, const Options& options) {
    return std::make_unique<Amesos2Solver>(solver, options);
  }  // end of makeAmesos2Solver

  template <typename Node>
  std::unique_ptr<Operator> TpetraBackend<Node>::makeBelosSolver(
      const std::string& solver,
      const Options& options,
      const BelosSettings& settings,
      const Operator* const p,
      ExternalPreconditioner e) {
    const auto* const tp = dynamic_cast<const TpetraOperatorBase*>(p);
    if ((p != nullptr) && (tp == nullptr)) {
      throw std::runtime_error(
          "Belos: the preconditioner is not an operator of Tpetra computing "
          "on the node of the solver");
    }
    return std::make_unique<BelosSolver>(solver, options, settings, tp,
                                         std::move(e));
  }  // end of makeBelosSolver

  extern template struct TpetraBackend<DefaultNode>;
#ifdef MFEM_MGIS_TRILINOS_HOST_NODE
  extern template struct TpetraBackend<HostNode>;
#endif /* MFEM_MGIS_TRILINOS_HOST_NODE */

}  // end of namespace mfem_mgis::trilinos

#endif /* LIB_MFEM_MGIS_SRC_TRILINOSTPETRABACKEND_HXX */
