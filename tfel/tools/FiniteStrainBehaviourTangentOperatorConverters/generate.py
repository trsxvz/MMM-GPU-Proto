"""Generates the explicit code of the conversions of the finite strain
tangent operators (structures `<source>_to_<target>_Converter` of
include/TFEL/Material/FiniteStrainBehaviourTangentOperatorConverters.ixx).

The conversions are described in conversions.py as the composition of
elementary relations between the strain measures and between the stress
measures. The code generated is explicit (no loop, constant indices) and
computes the result column by column.

usage:
  generate.py                    writes the header on the standard output
  generate.py --statistics       number of multiplications and of local
                                 variables of each conversion
  generate.py <source> <target>  code of one conversion
"""
import os
import subprocess
import sys

from conversions import OPERATORS, Conversion

# conversions generated (source, target)
CONVERSIONS = [
    ("DSIG_DF", "DSIG_DDF"),
    ("DSIG_DF", "DTAU_DDF"),
    ("DSIG_DF", "DPK1_DF"),
    ("DSIG_DF", "DS_DEGL"),
    ("DSIG_DF", "ABAQUS"),
    ("DSIG_DF", "C_TRUESDELL"),
    ("DSIG_DDF", "DSIG_DF"),
    ("DTAU_DF", "DTAU_DDF"),
    ("DTAU_DF", "DS_DEGL"),
    ("DTAU_DF", "SPATIAL_MODULI"),
    ("DTAU_DF", "C_TRUESDELL"),
    ("DTAU_DF", "ABAQUS"),
    ("DTAU_DDF", "DTAU_DF"),
    ("DTAU_DDF", "DSIG_DF"),
    ("DTAU_DDF", "DS_DEGL"),
    ("DTAU_DDF", "ABAQUS"),
    ("DTAU_DDF", "C_TRUESDELL"),
    ("DS_DC", "DS_DF"),
    ("DS_DC", "DSIG_DF"),
    ("DS_DC", "DTAU_DDF"),
    ("DS_DC", "DPK1_DF"),
    ("DS_DC", "ABAQUS"),
    ("DS_DC", "C_TRUESDELL"),
    ("DS_DEGL", "DS_DF"),
    ("DS_DEGL", "DSIG_DF"),
    ("DS_DEGL", "DTAU_DDF"),
    ("DS_DEGL", "DPK1_DF"),
    ("DS_DEGL", "SPATIAL_MODULI"),
    ("DS_DEGL", "ABAQUS"),
    ("DS_DEGL", "C_TRUESDELL"),
    ("SPATIAL_MODULI", "DTAU_DF"),
    ("SPATIAL_MODULI", "DTAU_DDF"),
    ("SPATIAL_MODULI", "DSIG_DF"),
    ("SPATIAL_MODULI", "DPK1_DF"),
    ("SPATIAL_MODULI", "DS_DEGL"),
    ("C_TRUESDELL", "DSIG_DF"),
    ("C_TAU_JAUMANN", "DTAU_DF"),
    ("ABAQUS", "DTAU_DF"),
    ("ABAQUS", "DSIG_DF"),
]

# clang-format style of TFEL
STYLE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..",
                     ".clang-format")

STRAIN_NAMES = {
    "F": "deformation gradient",
    "dF": "deformation gradient increment "
          "\\f$\\Delta\\,F=F\\,F_{0}^{-1}\\f$",
    "E": "Green-Lagrange strain",
    "C": "right Cauchy-Green tensor",
    "d": "rate of deformation",
}

# variation of the strain measures, as a function of the variation of the
# deformation gradient
STRAIN_VARIATIONS = {
    "F": "\\delta\\,F",
    "dF": "\\delta\\,\\Delta\\,F=\\delta\\,F\\,F_{0}^{-1}",
    "E": "\\delta\\,E=\\mathrm{sym}\\left(F^{T}\\,\\delta\\,F\\right)",
    "C": "\\delta\\,C=2\\,\\mathrm{sym}\\left(F^{T}\\,\\delta\\,F\\right)",
    "d": "d=\\mathrm{sym}\\left(l\\right)",
}

# variation of the stress measures, as a function of the variation of the
# Kirchhoff stress
STRESS_VARIATIONS = {
    "tau": "\\delta\\,\\tau",
    "sig": "\\delta\\,\\sigma=\\frac{1}{J}\\,\\delta\\,\\tau-"
           "\\mathrm{tr}\\left(l\\right)\\,\\sigma",
    "S": "\\delta\\,S=F^{-1}\\,\\left(\\delta\\,\\tau-l\\,\\tau-"
         "\\tau\\,l^{T}\\right)\\,F^{-T}",
    "P": "\\delta\\,P=\\left(\\delta\\,\\tau-\\tau\\,l^{T}\\right)\\,F^{-T}",
    "lie": "\\delta\\,\\tau-l\\,\\tau-\\tau\\,l^{T}",
    "lieJ": "\\frac{1}{J}\\,\\left(\\delta\\,\\tau-l\\,\\tau-"
            "\\tau\\,l^{T}\\right)",
    "jau": "\\delta\\,\\tau-w\\,\\tau+\\tau\\,w",
    "jauJ": "\\frac{1}{J}\\,\\left(\\delta\\,\\tau-w\\,\\tau+"
            "\\tau\\,w\\right)",
}

DIMENSIONS = {1: "1D", 2: "2D", 3: "3D"}


def wrap(text, prefix, width=78):
    """wraps a text in a comment"""
    lines = []
    line = prefix
    for word in text.split():
        if len(line) + 1 + len(word) > width and line != prefix:
            lines.append(line)
            line = prefix
        line += (" " if line != prefix or prefix.endswith("*") else "") + word
    lines.append(line)
    return "\n".join(l.rstrip() for l in lines)


def structure_name(source, target):
    return f"{source}_to_{target}_Converter"


def specialisation(source, target, N):
    """partial specialisation of the conversion structure"""
    c = Conversion(source, target, N)
    A, B = OPERATORS[source], OPERATORS[target]
    args = c.arguments()
    docs = {
        "F0": "\\param[in] F0: deformation gradient at the beginning of the "
              "time step",
        "F": "\\param[in] F: deformation gradient at the end of the time "
             "step",
        "s": "\\param[in] s: Cauchy stress",
    }
    types = {"F0": f"const tfel::math::tensor<{N}u, real>& F0",
             "F": f"const tfel::math::tensor<{N}u, real>& F",
             "s": f"const tfel::math::stensor<{N}u, stress>& s"}
    parameters = [f"tfel::math::{B.kind}<{N}u, stress>& Kr",
                  f"const tfel::math::{A.kind}<{N}u, stress>& Ks"] + \
        [types[a] for a in args]
    doc = ["\\tparam stress: stress type", "\\tparam real: numeric type",
           f"\\param[out] Kr: {B.description}",
           f"\\param[in] Ks: {A.description}"] + [docs[a] for a in args]
    out = [f"/*!\n * \\brief partial specialisation in {DIMENSIONS[N]}\n */",
           "template <>",
           f"struct {structure_name(source, target)}<{N}u> {{",
           "/*!"] + [wrap(d, " *") for d in doc] + [" */"]
    out += ["template <typename stress, typename real>",
            "TFEL_HOST_DEVICE static constexpr void exe(" +
            ", ".join(parameters) + ") noexcept {"]
    out += c.body()
    out += ["}", f"}};  // end of struct {structure_name(source, target)}"
            f"<{N}u>", ""]
    return "\n".join(out)


def structure(source, target):
    """declaration and partial specialisations of the conversion structure"""
    A, B = OPERATORS[source], OPERATORS[target]
    brief = (f"\\brief an helper structure meant to be specialized in each "
             f"spatial dimension, converting the {A.description} "
             f"({source}) to the {B.description} ({target}).")
    relations = [STRAIN_VARIATIONS[A.strain], STRESS_VARIATIONS[A.stress]]
    target_relations = [STRAIN_VARIATIONS[B.strain],
                        STRESS_VARIATIONS[B.stress]]
    text = (f"The result is computed column by column, i.e. component by "
            f"component of the {STRAIN_NAMES[B.strain]}. Let "
            f"\\f$\\delta\\,F\\f$ be the variation of the deformation "
            f"gradient associated with this component, "
            f"\\f$l=\\delta\\,F\\,F^{{-1}}\\f$, \\f$w\\f$ the skew part of "
            f"\\f$l\\f$ and \\f$\\delta\\,\\tau\\f$ the variation of the "
            f"Kirchhoff stress. {source} relates")
    out = ["/*!", wrap(brief, " *"), wrap(text, " *"),
           " * \\f[", wrap(relations[1], " *"), " * \\quad\\mathrm{to}\\quad",
           wrap(relations[0], " *"), " * \\f]",
           wrap(f"and {target} relates", " *"),
           " * \\f[", wrap(target_relations[1], " *"),
           " * \\quad\\mathrm{to}\\quad",
           wrap(target_relations[0], " *"), " * \\f]",
           wrap("so that no intermediate tangent operator is stored.", " *"),
           " * \\tparam N: spatial dimension", " */",
           "template <unsigned short N>",
           f"struct {structure_name(source, target)};", ""]
    return "\n".join(out) + "\n" + "\n".join(
        specialisation(source, target, N) for N in (1, 2, 3))


def clang_format(text):
    return subprocess.run(["clang-format", f"--style=file:{STYLE}"],
                          input=text, capture_output=True, text=True,
                          check=True).stdout


def header():
    """content of the header"""
    guard = ("LIB_TFEL_MATERIAL_FINITESTRAINBEHAVIOURTANGENTOPERATOR"
             "CONVERTERS_IXX")
    namespace = "tfel::material::internals"
    out = ["""/*!
 * \\file   FiniteStrainBehaviourTangentOperatorConverters.ixx
 * \\brief  conversions of the finite strain tangent operators computed
 * column by column. This file is generated by the script
 * tools/FiniteStrainBehaviourTangentOperatorConverters/generate.py.
 * \\date   22/11/2026
 * \\copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */
""", f"#ifndef {guard}\n#define {guard}\n",
           """#include "TFEL/Math/General/MathConstants.hxx"
#include "TFEL/Math/tensor.hxx"
#include "TFEL/Math/stensor.hxx"
#include "TFEL/Math/t2tot2.hxx"
#include "TFEL/Math/t2tost2.hxx"
#include "TFEL/Math/st2tost2.hxx"
""", f"namespace {namespace} {{\n"]
    for source, target in CONVERSIONS:
        out.append(structure(source, target))
    out.append(f"}}  // end of namespace {namespace}\n")
    out.append(f"#endif /* {guard} */")
    return clang_format("\n".join(out) + "\n")


def statistics():
    print(f"{'conversion':34s} {'1D':>10s} {'2D':>10s} {'3D':>10s}"
          "   (multiplications / local variables per column)")
    for source, target in CONVERSIONS:
        r = [Conversion(source, target, N).operations() for N in (1, 2, 3)]
        print(f"{source + ' -> ' + target:34s} " +
              " ".join(f"{m:>6d}/{v:<3d}" for m, v in r))


if __name__ == "__main__":
    if len(sys.argv) == 1:
        print(header(), end="")
    elif sys.argv[1:] == ["--statistics"]:
        statistics()
    elif len(sys.argv) == 3:
        print(clang_format(structure(sys.argv[1], sys.argv[2])), end="")
    else:
        sys.exit(__doc__)
