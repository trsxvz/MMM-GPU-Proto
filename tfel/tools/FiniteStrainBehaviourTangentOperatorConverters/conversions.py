"""Conversions of the finite strain tangent operators, expressed as the
composition of elementary relations, and generation of their explicit code.

A tangent operator relates the variation of a stress measure to the
variation of a strain measure. Let dF be a variation of the deformation
gradient F, l = dF.F^-1, d and w the symmetric and skew parts of l,
tau = J sig the Kirchhoff stress. The elementary relations are:

- strain measures (columns of the operators):
    F     dF
    dF    variation of the increment DF = F.F0^-1: dF.F0^-1
    E     Green-Lagrange strain: sym(F^T.dF) (C = 2 E + I)
    d     rate of deformation: sym(dF.F^-1)
- stress measures (rows of the operators), as functions of the variation
  dtau of the Kirchhoff stress:
    tau   dtau
    sig   dtau / J - tr(l) sig
    S     F^-1.(dtau - l.tau - tau.l^T).F^-T
    P     (dtau - tau.l^T).F^-T
    lie   dtau - l.tau - tau.l^T      (spatial moduli; / J: Truesdell rate
                                       of the Cauchy stress)
    jau   dtau - w.tau + tau.w        (Jaumann rate of the Kirchhoff
                                       stress; / J: Abaqus)

The conversion of the operator Ks to the operator Kr is computed column by
column: for each component of the strain measure of Kr,
1. the variation of the strain measure of Ks is computed (`strain_variation`);
2. Ks is applied to it (`apply`), which gives the variation of the stress
   measure of Ks;
3. the latter is converted to the variation of the stress measure of Kr
   through the variation of the Kirchhoff stress (`convert_stress`), and
   stored in the column of Kr.
Each column only needs a few local variables and no intermediate tangent
operator is stored.

Three choices reduce the number of operations:
- the components of the variation of the strain measure of Ks which are
  proportional to a component of a tensor are stored in local variables
  when this saves multiplications (see `apply`);
- tr(l) is computed from the derivative of the determinant of F when the
  inverse of F is not required by other terms (see
  `Direction.trace_stress`);
- in 3D, if the strain measure of Ks is a function of F^T.dF (E, C) and
  the one of Kr is DF, the variations of the stress measure of Ks
  associated with the components of a row of F are shared by the columns
  associated with the same row of DF (see `Conversion._rows`).

If Ks and Kr have the same stress measure, a row of Kr only depends on the
same row of Ks: the statements are ordered row by row, so that the
components of a row of Ks are only used by consecutive statements, which
requires fewer registers (see `Conversion._ordered_by_rows`).
"""
import re
from fractions import Fraction

from symbolic import (Block, Poly, T2, ISQRT2, SQRT2, code, stensor,
                      stensor_components, stensor_index, tensor,
                      tensor_components, tensor_index, unit_stensor,
                      unit_tensor)


class Operator:
    def __init__(self, stress, strain, kind, description):
        self.stress = stress
        self.strain = strain
        # TFEL type: t2tost2, st2tost2, t2tot2
        self.kind = kind
        self.description = description
        self.rows = "full" if kind == "t2tot2" else "sym"
        self.columns = "sym" if kind == "st2tost2" else "full"


OPERATORS = {
    "DSIG_DF": Operator(
        "sig", "F", "t2tost2",
        "derivative of the Cauchy stress with respect to the deformation "
        "gradient"),
    "DSIG_DDF": Operator(
        "sig", "dF", "t2tost2",
        "derivative of the Cauchy stress with respect to the deformation "
        "gradient increment"),
    "DTAU_DF": Operator(
        "tau", "F", "t2tost2",
        "derivative of the Kirchhoff stress with respect to the deformation "
        "gradient"),
    "DTAU_DDF": Operator(
        "tau", "dF", "t2tost2",
        "derivative of the Kirchhoff stress with respect to the deformation "
        "gradient increment"),
    "DS_DF": Operator(
        "S", "F", "t2tost2",
        "derivative of the second Piola-Kirchhoff stress with respect to the "
        "deformation gradient"),
    "DS_DC": Operator(
        "S", "C", "st2tost2",
        "derivative of the second Piola-Kirchhoff stress with respect to the "
        "right Cauchy-Green tensor"),
    "DS_DEGL": Operator(
        "S", "E", "st2tost2",
        "derivative of the second Piola-Kirchhoff stress with respect to the "
        "Green-Lagrange strain"),
    "SPATIAL_MODULI": Operator("lie", "d", "st2tost2", "spatial moduli"),
    "C_TRUESDELL": Operator(
        "lieJ", "d", "st2tost2",
        "moduli associated with the Truesdell rate of the Cauchy stress"),
    "C_TAU_JAUMANN": Operator(
        "jau", "d", "st2tost2",
        "moduli associated with the Jaumann rate of the Kirchhoff stress"),
    "ABAQUS": Operator(
        "jauJ", "d", "st2tost2",
        "moduli associated with the Jaumann rate of the Kirchhoff stress "
        "divided by the change of volume"),
    "DPK1_DF": Operator(
        "P", "F", "t2tot2",
        "derivative of the first Piola-Kirchhoff stress with respect to the "
        "deformation gradient"),
}

# stress measures whose variation is homogeneous to the one of the Kirchhoff
# stress (tau) or of the Cauchy stress (sig)
LEVEL = {"tau": "tau", "lie": "tau", "jau": "tau", "S": "tau", "P": "tau",
         "sig": "sig", "lieJ": "sig", "jauJ": "sig"}

class Environment:
    """quantities shared by the columns, defined before them"""

    def __init__(self, N, inverse):
        self.N = N
        # true if the inverse of the deformation gradient may be used to
        # compute the trace of l
        self.inverse = inverse
        self.lines = []
        self.arguments = set()
        self.defined = set()

    def _define(self, name, expression, type_=None):
        if name not in self.defined:
            self.defined.add(name)
            t = f"const tfel::math::{type_}" if type_ else "const auto"
            self.lines.append(f"{t} {name} = {expression};")

    def _tensor(self, real="real"):
        return f"tensor<{self.N}u, {real}>"

    def F(self):
        self.arguments.add("F")
        return tensor("F", self.N)

    def F0(self):
        self.arguments.add("F0")
        return tensor("F0", self.N)

    def sig(self):
        self.arguments.add("s")
        return stensor("s", self.N)

    def J(self):
        self.F()
        self._define("J", "tfel::math::det(F)")
        return Poly.atom("J")

    def iJ(self):
        self.J()
        self._define("iJ", "1 / J")
        return Poly.atom("iJ")

    def iF(self):
        self.F()
        self._define("iF", "tfel::math::invert(F)")
        return tensor("iF", self.N)

    def iF0(self):
        self.F0()
        self._define("iF0", "tfel::math::invert(F0)")
        return tensor("iF0", self.N)

    def dF(self):
        """deformation gradient increment"""
        self.F()
        self.iF0()
        self._define("dF", "F * iF0", self._tensor())
        return tensor("dF", self.N)

    def idF(self):
        """inverse of the deformation gradient increment"""
        self.F0()
        self.iF()
        self._define("idF", "F0 * iF", self._tensor())
        return tensor("idF", self.N)

    def iC(self):
        """inverse of the right Cauchy-Green tensor"""
        self.iF()
        self._define("iC", "tfel::math::computeLeftCauchyGreenTensor(iF)")
        return stensor("iC", self.N)

    def dJ(self):
        """derivative of the determinant: J F^-T"""
        self.F()
        self._define("dJ", "tfel::math::computeDeterminantDerivative(F)")
        return tensor("dJ", self.N)

    def dJF0(self):
        """product dJ.F0^T of the derivative of the determinant and of the
        transpose of the deformation gradient at the beginning of the time
        step"""
        self.dJ()
        self.F0()
        self._define("dJF0", "dJ * tfel::math::transpose(F0)", self._tensor())
        return tensor("dJF0", self.N)

    def tau(self):
        self.J()
        self.sig()
        self._define("tau", "J * s", f"stensor<{self.N}u, stress>")
        return stensor("tau", self.N)

    def P(self):
        self.F()
        self.sig()
        self._define(
            "P",
            "tfel::math::convertCauchyStressToFirstPiolaKirchhoffStress(s, F)")
        return tensor("P", self.N)

    def S(self):
        self.F()
        self.sig()
        self._define(
            "S",
            "tfel::math::convertCauchyStressToSecondPiolaKirchhoffStress(s, F)")
        return stensor("S", self.N)

    def stress(self, level):
        return self.tau() if level == "tau" else self.sig()

    def stress_G(self, level, G):
        """product X.G^T of the stress X (Kirchhoff or Cauchy) and of
        G = F^-1 or DF^-1"""
        self.P()
        if (level, G) == ("tau", "iF"):
            return tensor("P", self.N)
        if (level, G) == ("sig", "iF"):
            self.iJ()
            self._define("siF", "iJ * P", self._tensor("stress"))
            return tensor("siF", self.N)
        self.F0()
        if level == "tau":
            self._define("tidF", "P * tfel::math::transpose(F0)",
                         self._tensor("stress"))
            return tensor("tidF", self.N)
        self.iJ()
        self._define("sidF", "iJ * (P * tfel::math::transpose(F0))",
                     self._tensor("stress"))
        return tensor("sidF", self.N)


class Direction:
    """component of the strain measure of the result"""

    def __init__(self, strain, i, j, env):
        N = env.N
        self.strain = strain
        self.env = env
        self.rank_one = strain in ("F", "dF")
        if self.rank_one:
            self.column = tensor_index(N, i, j)
            self.E = unit_tensor(i, j)
            self.G = "iF" if strain == "F" else "idF"
            what = ("the deformation gradient" if strain == "F"
                    else "the deformation gradient increment")
            self.comment = (f"derivative with respect to the component "
                            f"({i}, {j}) of {what}")
        else:
            self.column = stensor_index(N, i, j)
            self.D = unit_stensor(i, j)
            if strain == "C":
                # the Green-Lagrange strain is half the right Cauchy-Green
                # tensor
                self.D = self.D * Fraction(1, 2)
            what = {"d": "the rate of deformation",
                    "E": "the Green-Lagrange strain",
                    "C": "the right Cauchy-Green tensor"}[strain]
            self.comment = (f"column associated with the component "
                            f"({i}, {j}) of {what}")

    def l(self):
        """l = dF.F^-1"""
        env = self.env
        if self.strain == "F":
            return self.E @ env.iF()
        if self.strain == "dF":
            return self.E @ env.idF()
        if self.strain == "d":
            return self.D
        iF = env.iF()
        return iF.T @ self.D @ iF

    def trace_stress(self, level):
        """tr(l) X, X being the Kirchhoff or the Cauchy stress"""
        env = self.env
        if not self.rank_one or env.inverse:
            return self.l().tr() * env.stress(level)
        # tr(l) = dJ : dF / J, where dJ = J F^-T is the derivative of the
        # determinant, which does not require the inverse of F
        dJ = env.dJ() if self.strain == "F" else env.dJF0()
        t = (self.E.T @ dJ).tr() * env.sig()
        return t if level == "tau" else env.iJ() * t

    def l_stress(self, level):
        """l.X, X being the Kirchhoff or the Cauchy stress"""
        if self.rank_one:
            return self.E @ self.env.stress_G(level, self.G).T
        return self.l() @ self.env.stress(level)


def directions(strain, env):
    c = tensor_components(env.N) if strain in ("F", "dF") \
        else stensor_components(env.N)
    return [Direction(strain, i, j, env) for i, j in c]


def sym(t):
    return (t + t.T) * Fraction(1, 2)


def strain_variation(source, d):
    """variation of the strain measure of the source in the direction d"""
    env = d.env
    target = d.strain
    if target in ("E", "C"):
        if source == "F":
            return env.iF().T @ d.D
        if source == "dF":
            return env.iF().T @ d.D @ env.iF0()
        if source in ("E", "C"):
            x = d.D
        else:
            x = d.l()
    elif target == "d":
        if source == "F":
            return d.D @ env.F()
        if source == "dF":
            return d.D @ env.dF()
        if source in ("E", "C"):
            F = env.F()
            x = F.T @ d.D @ F
        else:
            x = d.D
    else:
        if source == target:
            return d.E
        if source == "F":
            # target: dF
            return d.E @ env.F0()
        if source == "dF":
            return d.E @ env.iF0()
        if source in ("E", "C"):
            dF = d.E if target == "F" else d.E @ env.F0()
            x = sym(env.F().T @ dF)
        else:
            x = sym(d.l())
    return 2 * x if source == "C" else x


def apply(op, x, define, N):
    """variation of the stress measure of the operator `op` associated with
    the variation x of its strain measure. The components of x which are
    not reduced to an atom are stored in local variables by
    `define(c, v)`, c being the index of the component and v its value."""
    if op.columns == "full":
        xs = {tensor_index(N, i, j): x[i, j] for i, j in tensor_components(N)}
        existing = tensor_components(N)
    else:
        xs = {stensor_index(N, i, j): x[i, j] * (1 if i == j else SQRT2)
              for i, j in stensor_components(N)}
        existing = [(i, j) for i in range(3) for j in range(3)
                    if stensor_index(N, i, j) is not None]
    for i in range(3):
        for j in range(3):
            if (i, j) not in existing and not x[i, j].is_null():
                raise ValueError("invalid strain variation")
    # the components proportional to an atom are stored if this requires
    # less multiplications than the factorisation of their coefficients in
    # each row
    scalings = [v.scaling() for v in xs.values() if v.scaling() is not None]
    rows = len(tensor_components(N) if op.rows == "full"
               else stensor_components(N))
    scaled = len(scalings) < len(set(scalings)) * rows
    xs = {c: v if v.is_atomic() and not (scaled and v.scaling() is not None)
          else define(c, v) for c, v in xs.items()}
    a = T2.null()

    def row(r):
        v = Poly()
        for c, xc in xs.items():
            if not xc.is_null():
                v = v + Poly.atom(f"Ks({r}, {c})") * xc
        return v

    if op.rows == "full":
        for i, j in tensor_components(N):
            a.m[i][j] = row(tensor_index(N, i, j))
    else:
        for i, j in stensor_components(N):
            a.m[i][j] = row(stensor_index(N, i, j)) * (1 if i == j else ISQRT2)
            a.m[j][i] = a.m[i][j]
    return a


def sandwich(M, a, block, N):
    """M.a.M^T, a being symmetric. The product M.a is stored in local
    variables."""
    return block.tensor(M @ block.stensor(a, "a", N), "u", N) @ M.T


def convert_stress(source, target, a, d, block):
    """variation of the stress measure `target`, the variation `a` of the
    stress measure `source` being known in the direction d"""
    env = d.env
    N = env.N
    if target == "S":
        if source == "S":
            return a
        if source not in ("tau", "sig", "lie", "lieJ"):
            raise NotImplementedError(f"{source} -> S")
        if d.strain not in ("E", "C"):
            raise NotImplementedError(f"{source} -> S, {d.strain}")
        b = sandwich(env.iF(), a, block, N)
        if LEVEL[source] == "sig":
            b = env.J() * b
        if source in ("tau", "sig"):
            # F^-1.l.F = C^-1.dE
            M = env.iC() @ d.D
            S = env.S()
            MS = M @ S
            b = b - MS - MS.T
            if source == "sig":
                b = b + S * M.tr()
        return b
    if target == "P":
        if d.strain != "F":
            raise NotImplementedError(f"P, {d.strain}")
        if source == "tau":
            iF = env.iF()
            return block.stensor(a, "a", N) @ iF.T - env.P() @ d.E.T @ iF.T
        if source == "sig":
            # P = sig.dJ, where dJ = J F^-T is the derivative of the
            # determinant, whose derivative is linear in F
            F = env.F()
            ddJ = T2.null()
            for i in range(3):
                for j in range(3):
                    for p in range(3):
                        for q in range(3):
                            e = levi_civita(i, p, 3 - i - p) * \
                                levi_civita(j, q, 3 - j - q) \
                                if i != p and j != q else 0
                            if e != 0:
                                ddJ.m[i][j] = ddJ.m[i][j] + e * \
                                    F[3 - i - p, 3 - j - q] * d.E[p, q]
            return block.stensor(a, "a", N) @ env.dJ() + env.sig() @ ddJ
        if source == "S":
            return env.F() @ block.stensor(a, "a", N) + d.E @ env.S()
        if source in ("lie", "lieJ"):
            b = block.stensor(a, "a", N) @ env.iF().T
            if source == "lieJ":
                b = env.J() * b
            return b + d.E @ env.S()
        raise NotImplementedError(f"{source} -> P")
    # variation of the Kirchhoff stress (or of the Kirchhoff stress divided
    # by J), then of the target
    level = LEVEL[target]
    if source == "S":
        b = sandwich(env.F(), a, block, N)
    elif source == "P":
        b = a @ env.F().T
    else:
        b = a
    if LEVEL[source] != level:
        b = (env.iJ() if level == "sig" else env.J()) * b
    # terms proportional to the stress X (Kirchhoff or Cauchy): tr(l) X,
    # l.X + X.l^T and w.X - X.w
    trace = (source == "sig") - (target == "sig")
    lie = (source in ("S", "lie", "lieJ")) - (target in ("lie", "lieJ"))
    jaumann = (source in ("jau", "jauJ")) - (target in ("jau", "jauJ"))
    if trace != 0:
        b = b + trace * d.trace_stress(level)
    if lie != 0:
        lX = d.l_stress(level)
        b = b + lie * (lX + lX.T)
    if jaumann != 0 and d.rank_one:
        t = d.l_stress(level) - d.l().T @ env.stress(level)
        b = b + jaumann * (t + t.T) * Fraction(1, 2)
    if source == "P":
        if not d.rank_one:
            raise NotImplementedError("P")
        # tau.l^T
        b = b + env.stress_G(level, d.G) @ d.E.T
    return b


def levi_civita(i, j, k):
    return (i - j) * (j - k) * (k - i) // 2


def rows(op, b, N, lower):
    """values stored in the rows of a column of the result for the
    variation b of the stress measure. If `lower` is true, the lower part
    of b is used."""
    if op.rows == "full":
        r = [b[i, j] for i, j in tensor_components(N)]
        existing = tensor_components(N)
    else:
        r = [(b[j, i] if lower else b[i, j]) * (1 if i == j else SQRT2)
             for i, j in stensor_components(N)]
        existing = [(i, j) for i in range(3) for j in range(3)
                    if stensor_index(N, i, j) is not None]
    for i in range(3):
        for j in range(3):
            if (i, j) not in existing and not b[i, j].is_null():
                raise ValueError("invalid stress variation")
    return r


def store(op, b, column, block, N, lower):
    """stores the variation b of the stress measure in a column of the
    result"""
    for r, v in enumerate(rows(op, b, N, lower)):
        block.assign(f"Kr({r}, {column})", v)


class Conversion:
    """explicit code of the conversion of the operator `source` to the
    operator `target` in the space dimension N"""

    def __init__(self, source, target, N):
        self.source = source
        self.target = target
        self.N = N
        # the trace of l is computed from the derivative of the determinant,
        # unless the inverse of the deformation gradient is required by
        # other terms
        self._build(inverse=False)
        if "iF" in self.env.defined:
            self._build(inverse=True)

    def _build(self, inverse):
        N = self.N
        A, B = OPERATORS[self.source], OPERATORS[self.target]
        self.env = Environment(N, inverse)
        self.blocks = []
        if A.stress == B.stress:
            self._ordered_by_rows(A, B)
            return
        if N == 3 and B.strain == "dF" and A.strain in ("E", "C"):
            self._rows(A, B)
            return
        for d in directions(B.strain, self.env):
            block = Block()
            block.comment(d.comment)
            x = strain_variation(A.strain, d)
            a = apply(A, x, lambda c, v: block.define(f"x{c}", v), N)
            b = convert_stress(A.stress, B.stress, a, d, block)
            store(B, b, d.column, block, N, lower=(A.stress == "P"))
            self.blocks.append(block)

    def _ordered_by_rows(self, A, B):
        """statements ordered row by row, Ks and Kr having the same stress
        measure. The components of the variations of the strain measure of
        Ks which are stored in local variables are defined once, before
        the rows."""
        N = self.N
        block = Block()
        variables = {}

        def define(c, v):
            if code(v) not in variables:
                variables[code(v)] = block.define(f"x{len(variables)}", v)
            return variables[code(v)]

        columns = []
        for d in directions(B.strain, self.env):
            x = strain_variation(A.strain, d)
            a = apply(A, x, define, N)
            columns.append((d.column, rows(B, a, N, lower=False)))
        components = tensor_components(N) if B.rows == "full" \
            else stensor_components(N)
        for r, (i, j) in enumerate(components):
            block.comment(f"derivatives of the component ({i}, {j}) of the "
                          "stress")
            for column, values in columns:
                block.assign(f"Kr({r}, {column})", values[r])
        self.blocks.append(block)

    def _rows(self, A, B):
        """columns associated with the components of the deformation
        gradient increment DF, the strain measure of the source being a
        function of F^T.dF. As dF = dDF.F0, the variation of the stress
        measure of the source associated with a component of DF is a
        combination of the ones associated with the components of the same
        row of F, which are computed once for the three columns."""
        env, N = self.env, self.N
        F0 = env.F0()
        for i in range(3):
            columns = [k for p, k in tensor_components(N) if p == i]
            block = Block()
            block.comment(f"derivatives with respect to the components of "
                          f"the row {i} of the deformation gradient")
            r = {}
            for k in columns:
                x = strain_variation(A.strain, Direction("F", i, k, env))
                a = apply(A, x, lambda c, v: block.define(f"x{k}{c}", v), N)
                r[k] = block.stensor(a, f"r{k}", N)
            for j in columns:
                d = Direction(B.strain, i, j, env)
                column = Block()
                column.comment(d.comment)
                a = T2.null()
                for k in columns:
                    a = a + F0[j, k] * r[k]
                b = convert_stress(A.stress, B.stress, a, d, column)
                store(B, b, d.column, column, N, lower=False)
                block.scope(column)
            self.blocks.append(block)

    def arguments(self):
        """arguments used, among F0, F and s"""
        return [a for a in ("F0", "F", "s") if a in self.env.arguments]

    def uses(self, constant):
        return any(re.search(rf"\b{constant}\b", l)
                   for b in self.blocks for l in b.lines)

    def body(self):
        lines = []
        if self.uses("cste"):
            lines.append("constexpr auto cste = tfel::math::Cste<real>::sqrt2;")
        if self.uses("icste"):
            lines.append(
                "constexpr auto icste = tfel::math::Cste<real>::isqrt2;")
        lines += self.env.lines
        for b in self.blocks:
            # the local variables of a block are scoped, unless it is the
            # only one
            scoped = len(self.blocks) > 1 and any(
                l.startswith("const auto") for l in b.lines)
            lines += (["{"] if scoped else []) + b.lines + \
                (["}"] if scoped else [])
        return lines

    def operations(self):
        """number of multiplications and of local variables of the columns"""
        mults = 0
        variables = 0
        for b in self.blocks:
            # variables defined in the scopes opened
            defined = [0]
            for l in b.lines:
                if l == "{":
                    defined.append(0)
                elif l == "}":
                    defined.pop()
                else:
                    mults += l.count("*")
                    defined[-1] += l.startswith("const auto")
                    variables = max(variables, sum(defined))
        return mults, variables
