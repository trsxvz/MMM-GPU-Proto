/*!
 * \file   src/TrilinosTekoBackend.cxx
 * \brief  This file implements the block preconditioners of Teko (see
 * `makeTekoPreconditioner`), compiled by the compiler of Kokkos when
 * mfem-mgis uses Teko (option enable-teko)
 * \date   16/10/2026
 */

#include <string>
#include <vector>
#include <numeric>
#include <cctype>
#include <algorithm>
#include "MueLu_config.hpp"
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#if defined(HAVE_MUELU_STRATIMIKOS) && defined(HAVE_MUELU_THYRA)
#include "Stratimikos_MueLuHelpers.hpp"
#endif
#include "Teko_ConfigDefs.hpp"
#include "Teko_InverseLibrary.hpp"
#include "Teko_PreconditionerFactory.hpp"
#include "Teko_StridedTpetraOperator.hpp"
#include "Teko_TpetraInverseFactoryOperator.hpp"
#include "TrilinosTpetraBackend.hxx"

namespace mfem_mgis::trilinos {

  // Teko 17.2.1 is only instantiated for double, and for the ordinals and
  // the node of the default map of Tpetra (Teko_ConfigDefs.hpp)
  static_assert(std::is_same_v<Teko::ST, Scalar>,
                "Teko requires double precision");
  static_assert(std::is_same_v<Teko::LO, LO>);
  static_assert(std::is_same_v<Teko::GO, GO>);
  static_assert(std::is_same_v<Teko::NT, DefaultNode>);

  namespace {

    //! \brief backend of the default node of Tpetra, the node of Teko
    using Backend = TpetraBackend<DefaultNode>;
    //! \brief operator of Tpetra
    using OP = Backend::OP;

    //! \brief name of the inverse of the default library
    constexpr const char* default_inverse = "Block Jacobi";

    /*!
     * \return the default library of inverses: a block Jacobi
     * preconditioner whose diagonal blocks are approximated by one V-cycle
     * of the smoothed aggregation multigrid of MueLu, silent, with
     * reproducible aggregates
     */
    Teuchos::RCP<Teuchos::ParameterList> makeDefaultInverses() {
#if defined(HAVE_MUELU_STRATIMIKOS) && defined(HAVE_MUELU_THYRA)
      auto l = Teuchos::rcp(new Teuchos::ParameterList());
      auto& jacobi = l->sublist(default_inverse);
      jacobi.set("Type", "Block Jacobi");
      jacobi.set("Inverse Type", "MueLu of a block");
      auto& muelu = l->sublist("MueLu of a block");
      muelu.set("Type", "MueLu");
      muelu.set("verbosity", "none");
      muelu.set("aggregation: deterministic", true);
      return l;
#else  /* defined(HAVE_MUELU_STRATIMIKOS) && defined(HAVE_MUELU_THYRA) */
      throw std::runtime_error(
          "the parameter 'Parameters' is required (Trilinos is built without "
          "the adapter of MueLu for Stratimikos, used by the default "
          "inverse)");
#endif /* defined(HAVE_MUELU_STRATIMIKOS) && defined(HAVE_MUELU_THYRA) */
    }  // end of makeDefaultInverses

    //! \return the builder of Stratimikos used by Teko
    Teuchos::RCP<Stratimikos::DefaultLinearSolverBuilder>
    makeStratimikosBuilder() {
      auto builder =
          Teuchos::rcp(new Stratimikos::DefaultLinearSolverBuilder());
#if defined(HAVE_MUELU_STRATIMIKOS) && defined(HAVE_MUELU_THYRA)
      Stratimikos::enableMueLu<Teko::ST, Teko::LO, Teko::GO, Teko::NT>(
          *builder);
#endif /* defined(HAVE_MUELU_STRATIMIKOS) && defined(HAVE_MUELU_THYRA) */
      Teko::addToStratimikosBuilder(builder);
      return builder;
    }  // end of makeStratimikosBuilder

    /*!
     * \return the names of the entries of parameter lists
     * \param[in] lists: parameter lists
     */
    std::vector<std::string> getNames(
        const std::vector<const Teuchos::ParameterList*>& lists) {
      auto names = std::vector<std::string>{};
      for (const auto* const l : lists) {
        for (auto i = l->begin(); i != l->end(); ++i) {
          names.push_back(l->name(i));
        }
      }
      return names;
    }  // end of getNames

    /*!
     * \return the names, quoted and separated by commas
     * \param[in] names: names
     */
    std::string join(const std::vector<std::string>& names) {
      auto r = std::string{};
      for (const auto& n : names) {
        r += (r.empty() ? "'" : ", '") + n + "'";
      }
      return r;
    }  // end of join

    //! \brief names of the types of the inverses known to Teko
    struct InverseTypes {
      /*!
       * \param[in] builder: builder of Stratimikos
       */
      explicit InverseTypes(
          const Stratimikos::DefaultLinearSolverBuilder& builder) {
        Teko::PreconditionerFactory::getPreconditionerFactoryNames(
            this->blocks);
        const auto valid = builder.getValidParameters();
        this->preconditioners =
            getNames({&(valid->sublist("Preconditioner Types"))});
        this->solvers = getNames({&(valid->sublist("Linear Solver Types"))});
      }
      //! \return true if the type is a block preconditioner of Teko
      [[nodiscard]] bool isBlock(const std::string& t) const {
        return std::ranges::find(this->blocks, t) != this->blocks.end();
      }
      //! \return true if the type is a preconditioner of Stratimikos
      [[nodiscard]] bool isPreconditioner(const std::string& t) const {
        return std::ranges::find(this->preconditioners, t) !=
               this->preconditioners.end();
      }
      //! \return true if the type is a solver of Stratimikos
      [[nodiscard]] bool isSolver(const std::string& t) const {
        return std::ranges::find(this->solvers, t) != this->solvers.end();
      }
      //! \return true if the type is a solver or a preconditioner of
      //! Stratimikos
      [[nodiscard]] bool isStratimikos(const std::string& t) const {
        return this->isPreconditioner(t) || this->isSolver(t);
      }
      //! \return the names of the solvers and the preconditioners of
      //! Stratimikos
      [[nodiscard]] std::vector<std::string> getStratimikosNames() const {
        auto names = this->preconditioners;
        names.insert(names.end(), this->solvers.begin(), this->solvers.end());
        return names;
      }
      //! \brief names of the block preconditioners of Teko
      std::vector<std::string> blocks;
      //! \brief names of the preconditioners of Stratimikos
      std::vector<std::string> preconditioners;
      //! \brief names of the solvers of Stratimikos
      std::vector<std::string> solvers;
    };  // end of struct InverseTypes

    /*!
     * \return the names of the inverses of the library, separated by
     * commas: the definitions, then the solvers and the preconditioners of
     * Stratimikos which they do not redefine
     * \param[in] l: definitions of the inverses
     * \param[in] types: types known to Teko
     */
    std::string getInverseNames(const Teuchos::ParameterList& l,
                                const InverseTypes& types) {
      auto names = getNames({&l});
      for (const auto& n : types.getStratimikosNames()) {
        if (std::ranges::find(names, n) == names.end()) {
          names.push_back(n);
        }
      }
      return join(names);
    }  // end of getInverseNames

    /*!
     * \return the type of an inverse of the library: the type of its
     * definition, or its name for a solver or a preconditioner of
     * Stratimikos, or an empty string for an unknown inverse
     * \param[in] l: definitions of the inverses
     * \param[in] n: name of the inverse
     * \param[in] types: types known to Teko
     */
    std::string getInverseType(const Teuchos::ParameterList& l,
                               const std::string& n,
                               const InverseTypes& types) {
      if (l.isSublist(n)) {
        return l.sublist(n).get<std::string>("Type");
      }
      return types.isStratimikos(n) ? n : std::string{};
    }  // end of getInverseType

    /*!
     * \return the name of the inverse used as preconditioner: the name
     * given by the user, or the default inverse with the default
     * definitions, or the only block preconditioner of the definitions
     * \param[in] l: definitions of the inverses
     * \param[in] inverse: name given by the user, if any
     * \param[in] defaults: if true, the definitions are the default ones
     * \param[in] types: types known to Teko
     */
    std::string getInverseName(const Teuchos::ParameterList& l,
                               const std::string& inverse,
                               const bool defaults,
                               const InverseTypes& types) {
      if (!inverse.empty()) {
        return inverse;
      }
      if (defaults) {
        return default_inverse;
      }
      auto blocks = std::vector<std::string>{};
      for (auto i = l.begin(); i != l.end(); ++i) {
        const auto& n = l.name(i);
        if (types.isBlock(l.sublist(n).get<std::string>("Type"))) {
          blocks.push_back(n);
        }
      }
      if (blocks.empty()) {
        throw std::runtime_error(
            "the parameter 'Parameters' defines no block preconditioner of "
            "Teko");
      }
      if (blocks.size() != 1) {
        throw std::runtime_error(
            "the parameter 'Inverse' is required when 'Parameters' defines "
            "several block preconditioners of Teko (" +
            join(blocks) + ")");
      }
      return blocks.front();
    }  // end of getInverseName

    /*!
     * \return if the name of a parameter of a block preconditioner of Teko
     * is the given keyword, possibly followed by the index of a block
     * (from 1), as the parameters `Inverse Type` and `Inverse Type 2`
     * \param[in] p: name of the parameter
     * \param[in] k: keyword
     */
    bool isIndexedParameter(const std::string& p, const std::string& k) {
      if (p == k) {
        return true;
      }
      if ((!p.starts_with(k + " ")) || (p.size() == k.size() + 1)) {
        return false;
      }
      const auto index = p.substr(k.size() + 1);
      return std::ranges::all_of(
                 index,
                 [](const unsigned char c) { return std::isdigit(c) != 0; }) &&
             (index.front() != '0');
    }  // end of isIndexedParameter

    /*!
     * \brief check the definitions of the inverses, which Teko 17.2.1
     * reports by exceptions without explanation, or by crashes: each
     * definition is a sublist whose type is a block preconditioner of Teko
     * or a solver or a preconditioner of Stratimikos (a definition named
     * after a solver or a preconditioner of Stratimikos shall have this
     * type, Teko using the one of Stratimikos otherwise); the inverses
     * referenced by the parameters `Inverse Type`, `Inverse Type <i>`,
     * `Preconditioner Type` and `Preconditioner Type <i>` of the block
     * preconditioners, and the preconditioner given by the parameter `Use
     * Preconditioner` of the solvers of Stratimikos, are defined
     * \param[in] l: definitions of the inverses
     * \param[in] types: types known to Teko
     */
    void checkInverses(const Teuchos::ParameterList& l,
                       const InverseTypes& types) {
      const auto unknown = [&l, &types](const std::string& r,
                                        const std::string& p,
                                        const std::string& n) {
        return std::runtime_error(
            "unknown inverse '" + r + "' given by the parameter '" + p +
            "' of the inverse '" + n +
            "' (inverses: " + getInverseNames(l, types) + ")");
      };
      for (auto i = l.begin(); i != l.end(); ++i) {
        const auto& n = l.name(i);
        if (!l.isSublist(n)) {
          throw std::runtime_error("the entry '" + n +
                                   "' of the parameters is not the "
                                   "definition of an inverse (sublist)");
        }
        const auto& d = l.sublist(n);
        if (!d.isType<std::string>("Type")) {
          throw std::runtime_error("the inverse '" + n +
                                   "' has no type (string 'Type')");
        }
        const auto& t = d.get<std::string>("Type");
        if ((!types.isBlock(t)) && (!types.isStratimikos(t))) {
          throw std::runtime_error("unknown type '" + t + "' of the inverse '" +
                                   n + "' (types: " + join(types.blocks) +
                                   ", " + join(types.getStratimikosNames()) +
                                   ")");
        }
        if (types.isStratimikos(n) && (t != n)) {
          throw std::runtime_error(
              "the inverse '" + n +
              "' has the name of a solver or a preconditioner of Stratimikos "
              "but the type '" +
              t + "' (another name is required)");
        }
        if (types.isSolver(t) && d.isType<std::string>("Use Preconditioner")) {
          // the preconditioner of a solver of Stratimikos (Teko 17.2.1
          // crashes when it is not a preconditioner)
          const auto& r = d.get<std::string>("Use Preconditioner");
          if (r != "None") {
            const auto rt = getInverseType(l, r, types);
            if (rt.empty()) {
              throw unknown(r, "Use Preconditioner", n);
            }
            if ((!types.isBlock(rt)) && (!types.isPreconditioner(rt))) {
              throw std::runtime_error(
                  "the inverse '" + r +
                  "' given by the parameter 'Use Preconditioner' of the "
                  "inverse '" +
                  n + "' is not a preconditioner");
            }
          }
        }
        if (!types.isBlock(t)) {
          continue;
        }
        for (auto j = d.begin(); j != d.end(); ++j) {
          const auto& p = d.name(j);
          const auto preconditioner = p.starts_with("Preconditioner Type");
          if ((!p.starts_with("Inverse Type")) && (!preconditioner)) {
            continue;
          }
          // Teko 17.2.1 crashes on an invalid index of a block
          if ((!isIndexedParameter(p, "Inverse Type")) &&
              (!isIndexedParameter(p, "Preconditioner Type"))) {
            throw std::runtime_error(
                "invalid parameter '" + p + "' of the inverse '" + n +
                "' (expected 'Inverse Type' or 'Preconditioner Type', "
                "possibly followed by the index of a block, from 1)");
          }
          if (!d.isType<std::string>(p)) {
            throw std::runtime_error("the parameter '" + p +
                                     "' of the inverse '" + n +
                                     "' is not a string");
          }
          const auto& r = d.get<std::string>(p);
          if (preconditioner && (r == "None")) {
            continue;
          }
          if (getInverseType(l, r, types).empty()) {
            throw unknown(r, p, n);
          }
        }
      }
    }  // end of checkInverses

    /*!
     * \brief set the default parameters of mfem-mgis in the definitions of
     * the inverses of type `MueLu`, unless given, as for the preconditioner
     * `MueLu`: silent setup (verbosity `high` by default in MueLu),
     * reproducible aggregates and, with several equations per node, the
     * factories which are not written with Kokkos (`use kokkos refactor`
     * set to `false`, see the preconditioner `MueLu`). The preconditioner
     * `MueLu` of Stratimikos is defined, if it is not, so that it receives
     * these parameters too.
     * \param[in, out] l: definitions of the inverses
     * \param[in] types: types known to Teko
     */
    void setMueLuDefaultParameters(Teuchos::ParameterList& l,
                                   const InverseTypes& types) {
      if (types.isPreconditioner("MueLu") && (!l.isParameter("MueLu"))) {
        l.sublist("MueLu").set("Type", "MueLu");
      }
      for (auto i = l.begin(); i != l.end(); ++i) {
        auto& d = l.sublist(l.name(i));
        if (d.get<std::string>("Type") != "MueLu") {
          continue;
        }
        if (!d.isParameter("verbosity")) {
          d.set("verbosity", "none");
        }
        if (!d.isParameter("aggregation: deterministic")) {
          d.set("aggregation: deterministic", true);
        }
        const auto* const equations = "number of equations";
        if ((d.isType<int>(equations)) && (d.get<int>(equations) > 1) &&
            (!d.isParameter("use kokkos refactor"))) {
          d.set("use kokkos refactor", false);
        }
      }
    }  // end of setMueLuDefaultParameters

    //! \brief block preconditioner of Teko
    struct TekoPreconditioner final : public Backend::TpetraOperatorBase {
      /*!
       * \param[in] o: definitions of the inverses
       * \param[in] i: name of the inverse used as preconditioner
       * \param[in] b: number of components of each block
       */
      TekoPreconditioner(const Options& o,
                         const std::string& i,
                         std::vector<int> b)
          : inverses(o.empty() ? makeDefaultInverses() : makeParameterList(o)),
            blocks(std::move(b)) {
        // the library is built here, so that the invalid definitions and
        // names are refused at the creation of the preconditioner (the
        // preconditioners and the solvers of Stratimikos are added to the
        // library with their default parameters)
        const auto builder = makeStratimikosBuilder();
        const auto types = InverseTypes{*builder};
        checkInverses(*(this->inverses), types);
        this->inverse = getInverseName(*(this->inverses), i, o.empty(), types);
        setMueLuDefaultParameters(*(this->inverses), types);
        this->library = Teko::InverseLibrary::buildFromParameterList(
            *(this->inverses), builder);
        const auto d = this->library->getParameterList(this->inverse);
        if (d.is_null()) {
          throw std::runtime_error(
              "unknown inverse '" + this->inverse + "' (inverses: " +
              getInverseNames(*(this->inverses), types) + ")");
        }
        // the operator given to the inverse is the blocked matrix, which
        // the solvers and the preconditioners of Stratimikos do not accept
        const auto t = getInverseType(*(this->inverses), this->inverse, types);
        if (!types.isBlock(t)) {
          throw std::runtime_error(
              "the inverse '" + this->inverse +
              "' is not a block preconditioner of Teko (type '" + t +
              "'): the solvers and the preconditioners of Stratimikos are "
              "only used for the blocks (the entries Ifpack2, MueLu and "
              "Amesos2Preconditioner apply them to the whole matrix)");
        }
        this->factory = this->library->getInverseFactory(this->inverse);
      }
      void setMatrix(const DistributedMatrixView& m) override {
        // Teko 17.2.1 fails on a process owning no unknown (the strided
        // maps of the blocks are built from its smallest global index,
        // undefined) while the other ones wait in collective calls: checked
        // by all the processes, which fail together
        auto empty = static_cast<int>(m.local_size == 0);
        MPI_Allreduce(MPI_IN_PLACE, &empty, 1, MPI_INT, MPI_LOR,
                      m.communicator);
        if (empty != 0) {
          throw std::runtime_error(
              "each process must own unknowns (the strided blocks of Teko "
              "do not support empty processes)");
        }
        const auto vdim = m.number_of_components;
        auto b = this->blocks;
        if (b.empty()) {
          b.assign(static_cast<std::size_t>(vdim), 1);
        }
        if (std::reduce(b.begin(), b.end()) != vdim) {
          throw std::runtime_error(
              "the number of components of the blocks does not match the "
              "number of components of the unknowns");
        }
        this->A.set(m);
        // the unknowns of a node are consecutive (see `TpetraMatrix`): the
        // blocks are strided, Teko extracting the blocks of the matrix
        const auto a =
            Teuchos::rcp(new Teko::TpetraHelpers::StridedTpetraOperator(
                b, Teuchos::RCP<const OP>(this->A.matrix)));
        auto p = Teuchos::rcp(
            new Teko::TpetraHelpers::InverseFactoryOperator(this->factory));
        p->buildInverseOperator(Teuchos::RCP<const OP>(a));
        this->preconditioner = p;
      }
      [[nodiscard]] Teuchos::RCP<OP> getTpetraOperator() const override {
        return this->preconditioner;
      }

     private:
      //! \brief definitions of the inverses
      const Teuchos::RCP<Teuchos::ParameterList> inverses;
      //! \brief number of components of each block (empty: one per component)
      const std::vector<int> blocks;
      //! \brief name of the inverse used as preconditioner
      std::string inverse;
      /*!
       * \brief library of inverses, kept alive: the block preconditioners
       * built by its factories only hold a non owning pointer to it
       */
      Teuchos::RCP<Teko::InverseLibrary> library;
      //! \brief factory of the inverse used as preconditioner
      Teuchos::RCP<Teko::InverseFactory> factory;
      //! \brief preconditioner
      Teuchos::RCP<OP> preconditioner;
    };  // end of struct TekoPreconditioner

  }  // end of namespace

  std::unique_ptr<Operator> makeTekoPreconditioner(
      const Options& options,
      const std::string& inverse,
      const std::vector<int>& blocks,
      const bool) {
    return std::make_unique<TekoPreconditioner>(options, inverse, blocks);
  }  // end of makeTekoPreconditioner

}  // end of namespace mfem_mgis::trilinos
