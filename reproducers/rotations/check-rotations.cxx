/*!
 * \file   check-rotations.cxx
 * \brief  Independent check of the rotation of the tangent operator blocks
 * of an orthotropic behaviour of the `generic` interface: for any tangent
 * operator K in the material frame and any gradients x in the global frame,
 *
 *   rotateTangentOperatorBlocks(K) . x
 *     == rotateThermodynamicForces(K . rotateGradients(x)),
 *
 * since the gradients are rotated to the material frame and the
 * thermodynamic forces and the tangent operator to the global frame. The
 * rotations of the gradients and of the thermodynamic forces only use the
 * rotations of vectors, symmetric tensors and tensors, independently of the
 * code rotating the tangent operator blocks. The in place rotation
 * (destination equal to the source) is also checked.
 *
 * usage: check-rotations library behaviour hypothesis...
 */

#include <dlfcn.h>
#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <vector>
#include <algorithm>

using RotationFunction = void (*)(double*, const double*, const double*);

static void* library = nullptr;

static const void* symbol(const std::string& n) {
  return dlsym(library, n.c_str());
}

//! \return the number of components of a variable (see compare.cxx)
static std::size_t size(const int type, const std::string& h) {
  const auto d = (((type >> 3) & 3) != 0) ? ((type >> 3) & 3)
                 : (h == "Tridimensional")
                     ? 3
                     : (h.starts_with("AxisymmetricalGeneralised") ? 1 : 2);
  switch (type & 7) {
    case 0:
      return 1;
    case 1:
      return d == 1 ? 3 : (d == 2 ? 4 : 6);
    case 2:
      return static_cast<std::size_t>(d);
    case 3:
      return d == 1 ? 3 : (d == 2 ? 5 : 9);
  }
  std::fprintf(stderr, "unsupported type\n");
  std::exit(EXIT_FAILURE);
}

struct Variable {
  std::string name;
  std::size_t size, offset;
};

static std::vector<Variable> variables(const std::string& b,
                                       const std::string& g,
                                       const std::string& h) {
  const auto n = *static_cast<const unsigned short*>(symbol(b + "_n" + g));
  const auto* const names =
      static_cast<const char* const*>(symbol(b + "_" + g));
  const auto* const types = static_cast<const int*>(symbol(b + "_" + g + "Types"));
  auto r = std::vector<Variable>{};
  auto offset = std::size_t{};
  for (unsigned short i = 0; i != n; ++i) {
    r.push_back({names[i], size(types[i], h), offset});
    offset += r.back().size;
  }
  return r;
}

static const Variable& find(const std::vector<Variable>& vs,
                            const std::string& n) {
  const auto p = std::find_if(vs.begin(), vs.end(),
                              [&n](const Variable& v) { return v.name == n; });
  if (p == vs.end()) {
    std::fprintf(stderr, "unknown variable '%s'\n", n.c_str());
    std::exit(EXIT_FAILURE);
  }
  return *p;
}

//! \return a rotation matrix stored in column-major format
static std::vector<double> rotation(std::mt19937& g, const bool in_plane) {
  auto angle = std::uniform_real_distribution<double>(-M_PI, M_PI);
  const auto a = angle(g), b = in_plane ? 0. : angle(g),
             c = in_plane ? 0. : angle(g);
  // rotations about z, then x, then z
  const auto rz = [](const double t) {
    return std::vector<double>{std::cos(t), -std::sin(t), 0,
                               std::sin(t), std::cos(t),  0,
                               0,           0,            1};
  };
  const auto rx = [](const double t) {
    return std::vector<double>{1, 0,           0,
                               0, std::cos(t), -std::sin(t),
                               0, std::sin(t), std::cos(t)};
  };
  const auto product = [](const std::vector<double>& m1,
                          const std::vector<double>& m2) {
    auto m = std::vector<double>(9, 0.);
    for (int i = 0; i != 3; ++i) {
      for (int j = 0; j != 3; ++j) {
        for (int k = 0; k != 3; ++k) {
          m[3 * i + j] += m1[3 * i + k] * m2[3 * k + j];
        }
      }
    }
    return m;
  };
  const auto m = product(product(rz(a), rx(b)), rz(c));
  auto r = std::vector<double>(9);
  for (int i = 0; i != 3; ++i) {
    for (int j = 0; j != 3; ++j) {
      r[3 * j + i] = m[3 * i + j];  // column-major
    }
  }
  return r;
}

static bool check(const std::string& b, const std::string& h) {
  const auto f = b + "_" + h;
  const auto rg = reinterpret_cast<RotationFunction>(
      const_cast<void*>(symbol(f + "_rotateGradients")));
  const auto rf = reinterpret_cast<RotationFunction>(
      const_cast<void*>(symbol(f + "_rotateThermodynamicForces")));
  const auto rk = reinterpret_cast<RotationFunction>(
      const_cast<void*>(symbol(f + "_rotateTangentOperatorBlocks")));
  if ((rg == nullptr) || (rf == nullptr) || (rk == nullptr)) {
    std::printf("%s: missing rotation functions\n", f.c_str());
    return false;
  }
  const auto g = variables(b, "Gradients", h);
  const auto th = variables(b, "ThermodynamicForces", h);
  const auto ng = g.back().offset + g.back().size;
  const auto nf = th.back().offset + th.back().size;
  // tangent operator blocks, given as pairs of variables
  const auto nb =
      *static_cast<const unsigned short*>(symbol(b + "_nTangentOperatorBlocks"));
  const auto* const blocks =
      static_cast<const char* const*>(symbol(b + "_TangentOperatorBlocks"));
  auto layout = std::vector<std::pair<const Variable*, const Variable*>>{};
  auto nk = std::size_t{};
  for (unsigned short i = 0; i != nb / 2; ++i) {
    layout.emplace_back(&find(th, blocks[2 * i]), &find(g, blocks[2 * i + 1]));
    nk += layout.back().first->size * layout.back().second->size;
  }
  // full matrix of the tangent operator
  const auto assemble = [&](const std::vector<double>& k) {
    auto m = std::vector<double>(nf * ng, 0.);
    auto o = std::size_t{};
    for (const auto& [v1, v2] : layout) {
      for (std::size_t i = 0; i != v1->size; ++i) {
        for (std::size_t j = 0; j != v2->size; ++j) {
          m[(v1->offset + i) * ng + v2->offset + j] = k[o + i * v2->size + j];
        }
      }
      o += v1->size * v2->size;
    }
    return m;
  };
  const auto product = [&](const std::vector<double>& m,
                           const std::vector<double>& x) {
    auto y = std::vector<double>(nf, 0.);
    for (std::size_t i = 0; i != nf; ++i) {
      for (std::size_t j = 0; j != ng; ++j) {
        y[i] += m[i * ng + j] * x[j];
      }
    }
    return y;
  };
  auto generator = std::mt19937{42};
  auto value = std::uniform_real_distribution<double>(-1, 1);
  auto error = 0., in_place_error = 0.;
  for (int test = 0; test != 20; ++test) {
    const auto r = rotation(generator, h != "Tridimensional");
    auto km = std::vector<double>(nk);
    auto x = std::vector<double>(ng);
    std::generate(km.begin(), km.end(), [&] { return value(generator); });
    std::generate(x.begin(), x.end(), [&] { return value(generator); });
    // rotated tangent operator, out of place and in place
    auto kg = std::vector<double>(nk);
    rk(kg.data(), km.data(), r.data());
    auto kg2 = km;
    rk(kg2.data(), kg2.data(), r.data());
    // lhs = kg . x ; rhs = rf(km . rg(x))
    const auto lhs = product(assemble(kg), x);
    auto xm = std::vector<double>(ng);
    rg(xm.data(), x.data(), r.data());
    const auto fm = product(assemble(km), xm);
    auto rhs = std::vector<double>(nf);
    rf(rhs.data(), fm.data(), r.data());
    for (std::size_t i = 0; i != nf; ++i) {
      error = std::max(error, std::abs(lhs[i] - rhs[i]));
    }
    for (std::size_t i = 0; i != nk; ++i) {
      in_place_error = std::max(in_place_error, std::abs(kg[i] - kg2[i]));
    }
  }
  const auto ok = (error < 1e-12) && (in_place_error == 0);
  std::printf("%-60s %zu blocks | error %.1e, in place %.1e | %s\n",
              f.c_str(), layout.size(), error, in_place_error,
              ok ? "OK" : "FAILED");
  return ok;
}

int main(const int argc, const char* const* const argv) {
  if (argc < 4) {
    std::fprintf(stderr, "usage: %s library behaviour hypothesis...\n",
                 argv[0]);
    return EXIT_FAILURE;
  }
  library = dlopen(argv[1], RTLD_NOW);
  if (library == nullptr) {
    std::fprintf(stderr, "%s\n", dlerror());
    return EXIT_FAILURE;
  }
  auto ok = true;
  for (int i = 3; i != argc; ++i) {
    ok = check(argv[2], argv[i]) && ok;
  }
  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
