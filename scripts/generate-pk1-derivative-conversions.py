"""Code of the conversions of the derivative of the Kirchhoff stress with
respect to the deformation gradient (DTAU_DF) or to its increment
(DTAU_DDF) to the derivative of the first Piola-Kirchhoff stress with
respect to the deformation gradient (DPK1_DF), in the series of TFEL:

- DTAU_DF: structure DTAU_DF_to_DPK1_DF_Converter of
  include/TFEL/Math/T2toT2/ConvertToPK1Derivative.ixx;
- DTAU_DDF: structure DTAU_DDF_to_DPK1_DF_Converter of
  include/TFEL/Material/FiniteStrainBehaviourTangentOperator.ixx.

The code is explicit (no loop, constant indices) and computed column by
column, one column per component of the deformation gradient, so that each
column only needs a few intermediate values (6 in 3D) and that no
intermediate tangent operator is stored:

  P = tau . F^-T,   DF = F . F0^-1
  dtau/dF[ik,pq] = sum_m DTAU_DDF[ik,pm] F0^-1[q,m]
  DPK1_DF[iJ,pq] = sum_k dtau/dF[ik,pq] F^-1[J,k] - F^-1[J,p] P[i,q]

The bodies are formatted by clang-format with the style of TFEL.

usage: python3 scripts/generate-pk1-derivative-conversions.py DTAU_DF|DTAU_DDF
"""
import os
import subprocess
import sys

STYLE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tfel",
                     ".clang-format")
NAMES = {1: "1D", 2: "2D", 3: "3D"}


def ti(N, i, j):
    """index of the component (i, j) of a tensor, None if it does not exist
    in this space dimension"""
    t = [[0, 3, 5], [4, 1, 7], [6, 8, 2]][i][j]
    return t if t < {1: 3, 2: 5, 3: 9}[N] else None


def si(N, i, j):
    """index of the component (i, j) of a symmetric tensor"""
    t = [[0, 3, 4], [3, 1, 5], [4, 5, 2]][i][j]
    return t if t < {1: 3, 2: 4, 3: 6}[N] else None


def factor(i, k):
    """factor of the component (i, k) of a symmetric tensor"""
    return "" if i == k else "cste * "


def body(source, N):
    """explicit code of the conversion in the given space dimension"""
    columns = sorted(((p, q) for p in range(3) for q in range(3)
                      if ti(N, p, q) is not None), key=lambda pq: ti(N, *pq))
    srows = sorted({si(N, i, k) for i in range(3) for k in range(3)
                    if si(N, i, k) is not None})
    lines = []
    for p, q in columns:
        c = ti(N, p, q)
        if source == "DTAU_DDF":
            lines.append("{")
            lines.append("// derivative of the Kirchhoff stress with respect "
                         f"to F({p}, {q})")
            for r in srows:
                terms = [f"dtau({r}, {ti(N, p, m)}) * iF0[{ti(N, q, m)}]"
                         for m in range(3) if ti(N, p, m) is not None
                         and ti(N, q, m) is not None]
                lines.append(f"const auto a{r} = {' + '.join(terms)};")
            column = lambda r: f"a{r}"
        else:
            column = lambda r, c=c: f"dtau({r}, {c})"
        for i, J in columns:
            terms = [f"{factor(i, k)}{column(si(N, i, k))} * iF[{ti(N, J, k)}]"
                     for k in range(3) if si(N, i, k) is not None
                     and ti(N, J, k) is not None]
            if ti(N, J, p) is not None and ti(N, i, q) is not None:
                terms.append(f"- iF[{ti(N, J, p)}] * P[{ti(N, i, q)}]")
            e = " + ".join(terms).replace("+ -", "-")
            lines.append(f"dP({ti(N, i, J)}, {c}) = {e};")
        if source == "DTAU_DDF":
            lines.append("}")
    # formatted at the depth of the body of the exe method
    wrapped = ("namespace a {\nstruct b {\nvoid c() {\n// BEGIN\n" +
               "\n".join(lines) + "\n// END\n}\n};\n}\n")
    f = subprocess.run(["clang-format", f"--style=file:{STYLE}"],
                       input=wrapped, capture_output=True, text=True,
                       check=True).stdout.splitlines()
    i0 = next(i for i, l in enumerate(f) if "// BEGIN" in l) + 1
    i1 = next(i for i, l in enumerate(f) if "// END" in l)
    return "".join(l + "\n" for l in f[i0:i1])


def dtau_df():
    struct = "DTAU_DF_to_DPK1_DF_Converter"
    out = [f"""  /*!
   * \\brief an helper structure meant to be specialized in each spatial
   * dimension. The derivative of the first Piola-Kirchhoff stress
   * \\f$P=\\tau\\,F^{{-T}}\\f$ is computed column by column, i.e. component by
   * component of the deformation gradient:
   * \\f[
   * \\frac{{\\partial P_{{iJ}}}}{{\\partial F_{{pq}}}}=
   * \\frac{{\\partial \\tau_{{ik}}}}{{\\partial F_{{pq}}}}\\,F^{{-1}}_{{Jk}}-
   * F^{{-1}}_{{Jp}}\\,P_{{iq}}
   * \\f]
   * \\tparam N: spatial dimension
   */
  template <unsigned short N>
  struct {struct};
"""]
    for n in (1, 2, 3):
        cste = "" if n == 1 else (
            "      using real = numeric_type<DeformationGradientType>;\n"
            "      constexpr auto cste = Cste<real>::isqrt2;\n")
        out.append(f"""  /*!
   * \\brief partial specialisation in {NAMES[n]}
   */
  template <>
  struct {struct}<{n}u> {{
    /*!
     * \\param[out] dP: derivative of the first Piola-Kirchhoff stress with
     *                 respect to the deformation gradient
     * \\param[in] dtau: derivative of the Kirchhoff stress with respect to the
     *                  deformation gradient
     * \\param[in] F: deformation gradient
     * \\param[in] s: Cauchy stress
     */
    template <StressT2toT2Concept PK1StressDerivativeType,
              StressT2toST2Concept KirchhoffStressDerivativeType,
              NoUnitTensorConcept DeformationGradientType,
              StressStensorConcept CauchyStressType>
    requires(checkThatAllSpaceDimensionsAreEqual<PK1StressDerivativeType,
                                                 KirchhoffStressDerivativeType,
                                                 DeformationGradientType,
                                                 CauchyStressType>())  //
        TFEL_HOST_DEVICE static constexpr void                         //
        exe(PK1StressDerivativeType& dP,
            const KirchhoffStressDerivativeType& dtau,
            const DeformationGradientType& F,
            const CauchyStressType& s) noexcept {{
{cste}      const auto iF = invert(F);
      const auto P = convertCauchyStressToFirstPiolaKirchhoffStress(s, F);
{body("DTAU_DF", n)}    }}
  }};  // end of struct {struct}<{n}u>
""")
    return "\n".join(out)


def dtau_ddf():
    struct = "DTAU_DDF_to_DPK1_DF_Converter"
    out = [f"""  /*!
   * \\brief an helper structure meant to be specialized in each spatial
   * dimension, converting the derivative of the Kirchhoff stress with
   * respect to the deformation gradient increment
   * \\f$\\Delta\\,F=F\\,F_{{0}}^{{-1}}\\f$ to the derivative of the first
   * Piola-Kirchhoff stress \\f$P=\\tau\\,F^{{-T}}\\f$ with respect to the
   * deformation gradient. The latter is computed column by column, i.e.
   * component by component of the deformation gradient:
   * \\f[
   * \\frac{{\\partial \\tau_{{ik}}}}{{\\partial F_{{pq}}}}=
   * \\frac{{\\partial \\tau_{{ik}}}}{{\\partial \\Delta\\,F_{{pm}}}}\\,
   * F^{{-1}}_{{0\\,qm}}
   * \\quad\\mathrm{{and}}\\quad
   * \\frac{{\\partial P_{{iJ}}}}{{\\partial F_{{pq}}}}=
   * \\frac{{\\partial \\tau_{{ik}}}}{{\\partial F_{{pq}}}}\\,F^{{-1}}_{{Jk}}-
   * F^{{-1}}_{{Jp}}\\,P_{{iq}}
   * \\f]
   * so that the derivative of the Kirchhoff stress with respect to the
   * deformation gradient is never stored.
   * \\tparam N: spatial dimension
   */
  template <unsigned short N>
  struct {struct};
"""]
    for n in (1, 2, 3):
        cste = "" if n == 1 else (
            "      constexpr auto cste = tfel::math::Cste<real>::isqrt2;\n")
        out.append(f"""  /*!
   * \\brief partial specialisation in {NAMES[n]}
   */
  template <>
  struct {struct}<{n}u> {{
    /*!
     * \\tparam stress: stress type
     * \\tparam real: numeric type
     * \\param[out] dP: derivative of the first Piola-Kirchhoff stress with
     *                 respect to the deformation gradient
     * \\param[in] dtau: derivative of the Kirchhoff stress with respect to the
     *                  deformation gradient increment
     * \\param[in] F0: deformation gradient at the beginning of the time step
     * \\param[in] F: deformation gradient at the end of the time step
     * \\param[in] s: Cauchy stress
     */
    template <typename stress, typename real>
    TFEL_HOST_DEVICE static constexpr void exe(
        tfel::math::t2tot2<{n}u, stress>& dP,
        const tfel::math::t2tost2<{n}u, stress>& dtau,
        const tfel::math::tensor<{n}u, real>& F0,
        const tfel::math::tensor<{n}u, real>& F,
        const tfel::math::stensor<{n}u, stress>& s) noexcept {{
{cste}      const auto iF0 = tfel::math::invert(F0);
      const auto iF = tfel::math::invert(F);
      const auto P =
          tfel::math::convertCauchyStressToFirstPiolaKirchhoffStress(s, F);
{body("DTAU_DDF", n)}    }}
  }};  // end of struct {struct}<{n}u>
""")
    return "\n".join(out)


if __name__ == "__main__":
    if len(sys.argv) != 2 or sys.argv[1] not in ("DTAU_DF", "DTAU_DDF"):
        sys.exit(__doc__)
    print(dtau_df() if sys.argv[1] == "DTAU_DF" else dtau_ddf(), end="")
