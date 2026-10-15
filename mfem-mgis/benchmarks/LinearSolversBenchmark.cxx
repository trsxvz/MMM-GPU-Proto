/*!
 * \file   benchmarks/LinearSolversBenchmark.cxx
 * \brief  Benchmark of the linear solvers of the factory of linear solvers
 * on assembled systems, in parallel: the Laplace equation on a cube, or the
 * linear elasticity of a cantilever beam (length 4, square section of side
 * 1), possibly heterogeneous (cubic inclusions whose coefficients are
 * multiplied by a contrast). The linear solver and its parameters are given
 * on the command line, the parameters being written in a small text
 * syntax:
 *
 * \code
 * {Solver=CG,RelativeTolerance=1e-8,
 *  Preconditioner={Name=MueLu,Options={Parameters={"smoother:
 * type"=CHEBYSHEV}}}}
 * \endcode
 *
 * (a list of parameters between braces, a list of values between brackets,
 * strings quoted when they contain other characters than letters, digits,
 * `_`, `.`, `+` and `-`; `true` and `false` are booleans, and the numbers
 * are integers or reals). The program prints one line with the time of the
 * setup (`SetOperator`), the mean time of the resolutions (`Mult`, from a
 * null initial guess), the time of the setup for a second operator of the
 * same sparsity (the first one times two), the number of iterations, the
 * relative residual and the convergence reported. The time of the first
 * resolution after each setup is printed apart (`first`, `refirst`), as some
 * preconditioners are built at their first application (the ones of
 * hypre), and the relative residual of the resolution with the second
 * operator (`residual2`).
 * \date   15/10/2026
 */

#include <cmath>
#include <cctype>
#include <chrono>
#include <memory>
#include <string>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <string_view>
#include "mfem.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/LinearSolverFactory.hxx"
#include "MFEMMGIS/Utilities/SolverUtilities.hxx"

namespace {

  //! \brief parser of parameters written in the text syntax of this program
  struct ParametersParser {
    /*!
     * \return the parameters described by a text
     * \param[in] t: text
     */
    static mfem_mgis::Parameters parse(const std::string& t) {
      auto p = ParametersParser{t};
      p.skip();
      auto r = p.readParameters();
      p.skip();
      if (p.position != t.size()) {
        p.raise("unexpected characters after the parameters");
      }
      return r;
    }

   private:
    explicit ParametersParser(const std::string& t) : text(t) {}
    //! \brief report an error at the current position
    [[noreturn]] void raise(const std::string& msg) const {
      throw std::runtime_error("invalid parameters (" + msg + ", at position " +
                               std::to_string(this->position) + " of '" +
                               this->text + "')");
    }
    //! \brief skip the blanks
    void skip() {
      while ((this->position != this->text.size()) &&
             (std::isspace(static_cast<unsigned char>(
                  this->text[this->position])) != 0)) {
        ++(this->position);
      }
    }
    //! \return the current character, a null character at the end
    [[nodiscard]] char current() const {
      return this->position == this->text.size() ? '\0'
                                                 : this->text[this->position];
    }
    //! \brief read the given character
    void expect(const char c) {
      this->skip();
      if (this->current() != c) {
        this->raise(std::string{"'"} + c + "' expected");
      }
      ++(this->position);
    }
    //! \return a word or a quoted string
    std::string readWord() {
      this->skip();
      auto w = std::string{};
      if (this->current() == '"') {
        ++(this->position);
        while ((this->current() != '"') && (this->current() != '\0')) {
          w += this->current();
          ++(this->position);
        }
        this->expect('"');
        return w;
      }
      while ((std::isalnum(static_cast<unsigned char>(this->current())) != 0) ||
             ((this->current() != '\0') &&
              (std::string_view{"_.+-"}.find(this->current()) !=
               std::string_view::npos))) {
        w += this->current();
        ++(this->position);
      }
      if (w.empty()) {
        this->raise("a name or a value expected");
      }
      return w;
    }
    //! \return a value
    mfem_mgis::Parameter readValue() {
      this->skip();
      if (this->current() == '{') {
        return this->readParameters();
      }
      if (this->current() == '[') {
        ++(this->position);
        auto values = std::vector<mfem_mgis::Parameter>{};
        this->skip();
        if (this->current() != ']') {
          values.push_back(this->readValue());
          this->skip();
          while (this->current() == ',') {
            ++(this->position);
            values.push_back(this->readValue());
            this->skip();
          }
        }
        this->expect(']');
        return values;
      }
      const auto quoted = this->current() == '"';
      const auto w = this->readWord();
      if (quoted) {
        return w;
      }
      if (w == "true") {
        return true;
      }
      if (w == "false") {
        return false;
      }
      // integer, real, or string
      auto end = std::size_t{};
      try {
        const auto i = std::stoi(w, &end);
        if (end == w.size()) {
          return i;
        }
      } catch (std::exception&) {
      }
      try {
        const auto d = std::stod(w, &end);
        if (end == w.size()) {
          return d;
        }
      } catch (std::exception&) {
      }
      return w;
    }
    //! \return a list of parameters
    mfem_mgis::Parameters readParameters() {
      auto p = mfem_mgis::Parameters{};
      this->expect('{');
      this->skip();
      if (this->current() == '}') {
        ++(this->position);
        return p;
      }
      while (true) {
        const auto key = this->readWord();
        this->expect('=');
        p.insert(mfem_mgis::may_throw, key, this->readValue());
        this->skip();
        if (this->current() == ',') {
          ++(this->position);
          continue;
        }
        this->expect('}');
        return p;
      }
    }
    //! \brief text
    const std::string& text;
    //! \brief current position
    std::size_t position = 0;
  };  // end of struct ParametersParser

  //! \brief options of the benchmark
  struct Options {
    //! \brief problem: `laplace` or `elasticity`
    const char* problem = "elasticity";
    //! \brief number of elements across the section (cube or beam)
    int size = 8;
    //! \brief polynomial degree
    int order = 1;
    //! \brief if true, the unknowns of the elasticity are ordered by nodes
    bool by_nodes = true;
    //! \brief ratio of the coefficients of the inclusions and of the matrix
    double contrast = 1;
    //! \brief number of inclusions across the section
    int inclusions = 4;
    //! \brief name of the linear solver
    const char* solver = "CGSolver";
    //! \brief parameters of the linear solver
    const char* parameters = "{}";
    //! \brief label printed with the results
    const char* label = "";
    //! \brief number of resolutions
    int repeat = 3;
  };

  /*!
   * \brief system assembled by MFEM: the Laplace equation on the unit cube
   * (Dirichlet conditions on the whole boundary, unit source), or the
   * linear elasticity of a cantilever beam clamped on the face `x = 0`
   * (Lamé coefficients 1.5 and 1, vertical unit volume force)
   */
  struct System {
    explicit System(const Options& o) {
      const auto elasticity = std::string{o.problem} == "elasticity";
      if ((!elasticity) && (std::string{o.problem} != "laplace")) {
        throw std::runtime_error("invalid problem '" + std::string{o.problem} +
                                 "'");
      }
      const auto n = o.size;
      auto smesh =
          elasticity
              ? mfem::Mesh::MakeCartesian3D(
                    4 * n, n, n, mfem::Element::HEXAHEDRON, 4., 1., 1.)
              : mfem::Mesh::MakeCartesian3D(n, n, n, mfem::Element::HEXAHEDRON);
      // inclusions: attribute 2 for the elements in the cells of a
      // checkerboard of `inclusions` cells across the section
      for (int e = 0; e != smesh.GetNE(); ++e) {
        auto c = mfem::Vector(3);
        smesh.GetElementCenter(e, c);
        const auto k = static_cast<double>(o.inclusions);
        const auto i = static_cast<int>(std::floor(c(0) * k)) +
                       static_cast<int>(std::floor(c(1) * k)) +
                       static_cast<int>(std::floor(c(2) * k));
        smesh.SetAttribute(e, (i % 2 == 0) ? 1 : 2);
      }
      smesh.SetAttributes();
      this->mesh = std::make_unique<mfem::ParMesh>(MPI_COMM_WORLD, smesh);
      this->fec = std::make_unique<mfem::H1_FECollection>(o.order, 3);
      const auto vdim = elasticity ? 3 : 1;
      this->fespace = std::make_unique<mfem::ParFiniteElementSpace>(
          this->mesh.get(), this->fec.get(), vdim,
          o.by_nodes ? mfem::Ordering::byNODES : mfem::Ordering::byVDIM);
      auto scaling = mfem::Vector(2);
      scaling(0) = 1;
      scaling(1) = o.contrast;
      auto lambda_values = scaling;
      lambda_values *= 1.5;
      this->lambda.UpdateConstants(lambda_values);
      this->mu.UpdateConstants(scaling);
      auto& lambda = this->lambda;
      auto& mu = this->mu;
      // the matrix refers to the one of the bilinear form, which is kept
      this->form = std::make_unique<mfem::ParBilinearForm>(this->fespace.get());
      auto& a = *(this->form);
      auto b = mfem::ParLinearForm(this->fespace.get());
      auto boundaries = mfem::Array<int>(this->mesh->bdr_attributes.Max());
      boundaries = 0;
      auto f = mfem::Vector(3);
      f = 0.;
      f(2) = -1;
      auto fc = mfem::VectorConstantCoefficient(f);
      auto one = mfem::ConstantCoefficient(1.);
      if (elasticity) {
        a.AddDomainIntegrator(new mfem::ElasticityIntegrator(lambda, mu));
        b.AddDomainIntegrator(new mfem::VectorDomainLFIntegrator(fc));
        boundaries[4] = 1;
      } else {
        a.AddDomainIntegrator(new mfem::DiffusionIntegrator(mu));
        b.AddDomainIntegrator(new mfem::DomainLFIntegrator(one));
        boundaries = 1;
      }
      a.Assemble();
      b.Assemble();
      auto x = mfem::ParGridFunction(this->fespace.get());
      x = 0.;
      auto ess_tdofs = mfem::Array<int>{};
      this->fespace->GetEssentialTrueDofs(boundaries, ess_tdofs);
      auto X = mfem::Vector{};
      a.FormLinearSystem(ess_tdofs, x, b, this->A, X, this->B);
    }
    /*!
     * \return the norm of the residual relative to the right hand side
     * \param[in] m: matrix
     * \param[in] x: solution
     */
    double getRelativeResidual(const mfem::HypreParMatrix& m,
                               const mfem::Vector& x) const {
      auto r = mfem::Vector(this->B.Size());
      r.UseDevice(true);
      m.Mult(x, r);
      r -= this->B;
      const auto comm = this->fespace->GetComm();
      return std::sqrt(mfem::InnerProduct(comm, r, r) /
                       mfem::InnerProduct(comm, this->B, this->B));
    }
    //! \brief mesh
    std::unique_ptr<mfem::ParMesh> mesh;
    //! \brief finite element collection
    std::unique_ptr<mfem::FiniteElementCollection> fec;
    //! \brief finite element space
    std::unique_ptr<mfem::ParFiniteElementSpace> fespace;
    //! \brief first Lamé coefficient (elasticity)
    mfem::PWConstCoefficient lambda;
    //! \brief second Lamé coefficient, or conductivity
    mfem::PWConstCoefficient mu;
    //! \brief bilinear form
    std::unique_ptr<mfem::ParBilinearForm> form;
    //! \brief matrix
    mfem::HypreParMatrix A;
    //! \brief right hand side
    mfem::Vector B;
  };  // end of struct System

  //! \brief wait for the processes and for the device
  void synchronize() {
#if defined(MFEM_USE_CUDA) || defined(MFEM_USE_HIP)
    if (mfem::Device::Allows(mfem::Backend::CUDA_MASK |
                             mfem::Backend::HIP_MASK)) {
      MFEM_DEVICE_SYNC;
    }
#endif /* defined(MFEM_USE_CUDA) || defined(MFEM_USE_HIP) */
    MPI_Barrier(MPI_COMM_WORLD);
  }  // end of synchronize

  /*!
   * \return the duration of a function, in seconds, the processes and the
   * device being synchronized before and after
   * \param[in] f: function
   */
  template <typename Function>
  double measure(const Function& f) {
    synchronize();
    const auto start = std::chrono::steady_clock::now();
    f();
    synchronize();
    const auto end = std::chrono::steady_clock::now();
    return std::chrono::duration<double>(end - start).count();
  }  // end of measure

}  // end of namespace

int main(int argc, char** argv) {
  using namespace mfem_mgis;
  initialize(argc, argv);
  auto o = Options{};
  auto args = mfem::OptionsParser(argc, argv);
  declareDefaultOptions(args);
  args.AddOption(&o.problem, "-p", "--problem",
                 "Problem: laplace or elasticity.");
  args.AddOption(&o.size, "-n", "--size",
                 "Number of elements across the section.");
  args.AddOption(&o.order, "-o", "--order", "Polynomial degree.");
  args.AddOption(&o.by_nodes, "-bn", "--by-nodes", "-bv", "--by-vdim",
                 "Ordering of the unknowns of the elasticity.");
  args.AddOption(&o.contrast, "-c", "--contrast",
                 "Ratio of the coefficients of the inclusions.");
  args.AddOption(&o.inclusions, "-i", "--inclusions",
                 "Number of inclusions across the section.");
  args.AddOption(&o.solver, "-s", "--solver", "Name of the linear solver.");
  args.AddOption(&o.parameters, "-sp", "--solver-parameters",
                 "Parameters of the linear solver.");
  args.AddOption(&o.label, "-lb", "--label", "Label printed with the results.");
  args.AddOption(&o.repeat, "-rp", "--repeat", "Number of resolutions.");
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(std::cerr);
    return EXIT_FAILURE;
  }
  auto rank = 0;
  auto nprocs = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  try {
    const auto parameters = ParametersParser::parse(o.parameters);
    auto s = System(o);
    auto ctx = Context{};
    auto ls = LinearSolverFactory<true>::getFactory().generate(
        ctx, o.solver, *(s.fespace), parameters);
    if (isInvalid(ls)) {
      throw std::runtime_error(ctx.getErrorMessage());
    }
    auto& solver = *(ls.linear_solver);
    auto x = mfem::Vector(s.B.Size());
    x.UseDevice(true);
    const auto resolve = [&](const mfem::HypreParMatrix& m) {
      x = 0.;
      return measure([&] { solver.Mult(s.B, x); });
    };
    // the first resolution may include a part of the setup (the
    // preconditioners of hypre are built at their first application)
    const auto setup = measure([&] { solver.SetOperator(s.A); });
    const auto first = resolve(s.A);
    auto solve = 0.;
    for (int i = 1; i < o.repeat; ++i) {
      solve += resolve(s.A);
    }
    solve = o.repeat > 1 ? solve / (o.repeat - 1) : first;
    const auto converged = hasConverged(solver);
    const auto iterations =
        getNumberOfIterationsAtConvergence(solver).value_or(-1);
    const auto residual = s.getRelativeResidual(s.A, x);
    // second operator of the same sparsity
    auto a2 = mfem::HypreParMatrix(s.A);
    a2 *= 2.;
    const auto resetup = measure([&] { solver.SetOperator(a2); });
    const auto refirst = resolve(a2);
    const auto residual2 = s.getRelativeResidual(a2, x);
    if (rank == 0) {
      std::cout << "label=" << o.label << " problem=" << o.problem
                << " unknowns=" << s.fespace->GlobalTrueVSize()
                << " processes=" << nprocs << " solver=" << o.solver
                << " setup=" << setup << " first=" << first
                << " solve=" << solve << " resetup=" << resetup
                << " refirst=" << refirst << " iterations=" << iterations
                << " residual=" << residual << " residual2=" << residual2
                << " converged=" << (converged ? 1 : 0) << '\n';
    }
  } catch (std::exception& e) {
    std::cerr << "error: " << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}  // end of main
