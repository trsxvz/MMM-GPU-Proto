"""Small symbolic layer used to generate the explicit code of the
conversions of the finite strain tangent operators.

A second order tensor (`T2`) is a 3x3 array of polynomials (`Poly`) of
atoms, an atom being a piece of C++ code (a component of a tensor of TFEL,
a local variable). The coefficients of the polynomials are numbers of the
form a + b sqrt(2) (`Coef`), which is what the storage of the symmetric
tensors of TFEL requires. Components which do not exist in the space
dimension considered are null polynomials, so that the usual tensor
algebra (products, transposition, symmetric part...) directly gives the
sums restricted to the existing components.
"""
import re
from fractions import Fraction

# atoms which are factorised as scalar factors when writing a polynomial
SCALARS = ("J", "iJ")


def tensor_index(N, i, j):
    """index of the component (i, j) of a tensor, None if it does not exist
    in this space dimension"""
    t = [[0, 3, 5], [4, 1, 7], [6, 8, 2]][i][j]
    return t if t < {1: 3, 2: 5, 3: 9}[N] else None


def stensor_index(N, i, j):
    """index of the component (i, j) of a symmetric tensor"""
    t = [[0, 3, 4], [3, 1, 5], [4, 5, 2]][i][j]
    return t if t < {1: 3, 2: 4, 3: 6}[N] else None


def tensor_components(N):
    """components (i, j) of a tensor, in the order of their storage"""
    c = [(i, j) for i in range(3) for j in range(3)
         if tensor_index(N, i, j) is not None]
    return sorted(c, key=lambda ij: tensor_index(N, *ij))


def stensor_components(N):
    """components (i, j), i <= j, of a symmetric tensor, in the order of
    their storage"""
    c = [(i, j) for i in range(3) for j in range(i, 3)
         if stensor_index(N, i, j) is not None]
    return sorted(c, key=lambda ij: stensor_index(N, *ij))


class Coef:
    """number of the form a + b sqrt(2), a and b being rational"""
    __slots__ = ("a", "b")

    def __init__(self, a=0, b=0):
        self.a = Fraction(a)
        self.b = Fraction(b)

    def __add__(self, o):
        return Coef(self.a + o.a, self.b + o.b)

    def __mul__(self, o):
        return Coef(self.a * o.a + 2 * self.b * o.b,
                    self.a * o.b + self.b * o.a)

    def __neg__(self):
        return Coef(-self.a, -self.b)

    def __eq__(self, o):
        return self.a == o.a and self.b == o.b

    def is_null(self):
        return self.a == 0 and self.b == 0

    def split(self):
        """(sign, q, k) such that the coefficient is sign * q * sqrt(2)^k"""
        if self.a != 0 and self.b != 0:
            raise ValueError("unsupported coefficient")
        q, k = (self.b, 1) if self.a == 0 else (self.a, 0)
        return (-1 if q < 0 else 1), abs(q), k


ONE = Coef(1)
SQRT2 = Coef(0, 1)
ISQRT2 = Coef(0, Fraction(1, 2))


def as_coef(c):
    return c if isinstance(c, Coef) else Coef(c)


class Poly:
    """polynomial of atoms"""
    __slots__ = ("t",)

    def __init__(self, terms=None):
        self.t = {m: c for m, c in (terms or {}).items() if not c.is_null()}

    @staticmethod
    def atom(name):
        return Poly({(name,): ONE})

    def is_null(self):
        return not self.t

    def is_atomic(self):
        """null, or a single monomial made of at most one atom"""
        return len(self.t) <= 1 and all(len(m) <= 1 for m in self.t)

    def scaling(self):
        """(q, k) if the polynomial is an atom multiplied by a coefficient
        q * sqrt(2)^k other than 1 and -1, None otherwise"""
        if len(self.t) != 1:
            return None
        (m, c), = self.t.items()
        _, q, k = c.split()
        return None if len(m) != 1 or (q, k) == (1, 0) else (q, k)

    def __add__(self, o):
        r = dict(self.t)
        for m, c in o.t.items():
            r[m] = r[m] + c if m in r else c
        return Poly(r)

    def __neg__(self):
        return Poly({m: -c for m, c in self.t.items()})

    def __sub__(self, o):
        return self + (-o)

    def __mul__(self, o):
        if isinstance(o, T2):
            return NotImplemented
        if not isinstance(o, Poly):
            o = as_coef(o)
            return Poly({m: c * o for m, c in self.t.items()})
        r = {}
        for m1, c1 in self.t.items():
            for m2, c2 in o.t.items():
                m = tuple(sorted(m1 + m2))
                c = c1 * c2
                r[m] = r[m] + c if m in r else c
        return Poly(r)

    __rmul__ = __mul__



def _factor(q, k):
    """C++ factors (prefix, suffix) of the positive coefficient
    q * sqrt(2)^k"""
    name = ""
    if k == 1:
        # q sqrt(2) = (2 q) / sqrt(2)
        if q.denominator % 2 == 0:
            name, q = "icste", 2 * q
        else:
            name = "cste"
    prefix = [str(q.numerator)] if q.numerator != 1 else []
    if name:
        prefix.append(name)
    suffix = f" / {q.denominator}" if q.denominator != 1 else ""
    return prefix, suffix


def _atom_key(a):
    # the components of the tangent operator converted first, then the
    # local variables of the columns, then the other quantities
    if a.startswith("Ks("):
        return (0, a)
    return (1 if re.fullmatch(r"[a-z]\d+", a) else 2, a)


def code(p):
    """C++ expression of a polynomial. The monomials sharing the same
    scalar factors (J, 1 / J), then the same coefficient, are grouped."""
    if p.is_null():
        raise ValueError("null polynomial")
    monomials = []
    for m, c in p.t.items():
        scalars = tuple(a for a in m if a in SCALARS)
        atoms = sorted((a for a in m if a not in SCALARS), key=_atom_key)
        monomials.append((scalars, atoms, c))
    monomials.sort(key=lambda m: (min([_atom_key(a)[0] for a in m[1]],
                                      default=3),
                                  [_atom_key(a) for a in m[1]]))
    groups = {}
    for scalars, atoms, c in monomials:
        groups.setdefault(scalars, []).append((atoms, c))
    terms = []  # (sign, text)
    for scalars, ms in groups.items():
        inner = _coefficient_groups(ms)
        if not scalars:
            terms += inner
        elif len(inner) == 1:
            sign, text = inner[0]
            terms.append((sign, " * ".join(
                list(scalars) + ([text] if text != "1" else []))))
        else:
            sign = inner[0][0]
            text = _join([(s * sign, t) for s, t in inner])
            terms.append((sign, " * ".join(list(scalars) + [f"({text})"])))
    return _join(terms)


def _coefficient_groups(monomials):
    """terms (sign, text) of a list of monomials (atoms, coefficient), the
    monomials sharing the same coefficient being grouped"""
    groups = {}
    for atoms, c in monomials:
        sign, q, k = c.split()
        groups.setdefault((q, k), []).append((sign, atoms))
    terms = []
    for (q, k), ms in groups.items():
        prefix, suffix = _factor(q, k)
        if not prefix and not suffix:
            for sign, atoms in ms:
                terms.append((sign, " * ".join(atoms) if atoms else "1"))
            continue
        # the sign of the group is the one of its first monomial
        sign = ms[0][0]
        if len(ms) == 1:
            atoms = ms[0][1]
            text = " * ".join(prefix + atoms) if (prefix or atoms) else "1"
        else:
            inner = _join([(s * sign, " * ".join(a) if a else "1")
                           for s, a in ms])
            text = " * ".join(prefix + [f"({inner})"])
        terms.append((sign, text + suffix))
    return terms


def _join(terms):
    r = ""
    for i, (sign, text) in enumerate(terms):
        if i == 0:
            r = text if sign > 0 else f"-{text}"
        else:
            r += f" + {text}" if sign > 0 else f" - {text}"
    return r


class T2:
    """second order tensor: 3x3 array of polynomials"""

    def __init__(self, m):
        self.m = m

    @staticmethod
    def null():
        return T2([[Poly() for _ in range(3)] for _ in range(3)])

    def __getitem__(self, ij):
        return self.m[ij[0]][ij[1]]

    def __matmul__(self, o):
        r = T2.null()
        for i in range(3):
            for j in range(3):
                v = Poly()
                for k in range(3):
                    if not (self.m[i][k].is_null() or o.m[k][j].is_null()):
                        v = v + self.m[i][k] * o.m[k][j]
                r.m[i][j] = v
        return r

    def _map(self, f):
        return T2([[f(i, j) for j in range(3)] for i in range(3)])

    @property
    def T(self):
        return self._map(lambda i, j: self.m[j][i])

    def __add__(self, o):
        return self._map(lambda i, j: self.m[i][j] + o.m[i][j])

    def __sub__(self, o):
        return self._map(lambda i, j: self.m[i][j] - o.m[i][j])

    def __neg__(self):
        return self._map(lambda i, j: -self.m[i][j])

    def __mul__(self, s):
        return self._map(lambda i, j: self.m[i][j] * s)

    __rmul__ = __mul__

    def tr(self):
        return self.m[0][0] + self.m[1][1] + self.m[2][2]



def tensor(name, N):
    """tensor of TFEL"""
    r = T2.null()
    for i, j in tensor_components(N):
        r.m[i][j] = Poly.atom(f"{name}[{tensor_index(N, i, j)}]")
    return r


def stensor(name, N):
    """symmetric tensor of TFEL"""
    r = T2.null()
    for i, j in stensor_components(N):
        v = Poly.atom(f"{name}[{stensor_index(N, i, j)}]")
        r.m[i][j] = v if i == j else v * ISQRT2
        r.m[j][i] = r.m[i][j]
    return r


def unit_tensor(p, q):
    """tensor whose only non zero component is (p, q)"""
    r = T2.null()
    r.m[p][q] = Poly({(): ONE})
    return r


def unit_stensor(k, l):
    """symmetric tensor associated with the component (k, l) of the storage
    of the symmetric tensors of TFEL"""
    r = T2.null()
    if k == l:
        r.m[k][k] = Poly({(): ONE})
    else:
        r.m[k][l] = r.m[l][k] = Poly({(): ISQRT2})
    return r


class Block:
    """sequence of statements"""

    def __init__(self):
        self.lines = []

    def comment(self, text):
        self.lines.append(f"// {text}")

    def define(self, name, p):
        """defines a local variable and returns the associated atom"""
        self.lines.append(f"const auto {name} = {code(p)};")
        return Poly.atom(name)

    def assign(self, lhs, p):
        self.lines.append(f"{lhs} = {code(p)};")

    def scope(self, block):
        """appends the statements of a block in a nested scope"""
        self.lines += ["{"] + block.lines + ["}"]

    def tensor(self, t, prefix, N):
        """stores the components of a tensor in local variables, except the
        ones which are null or reduced to an atom"""
        r = T2.null()
        for i, j in tensor_components(N):
            v = t[i, j]
            r.m[i][j] = v if v.is_atomic() else self.define(
                f"{prefix}{tensor_index(N, i, j)}", v)
        return r

    def stensor(self, t, prefix, N):
        """stores the components of a symmetric tensor in local variables"""
        r = T2.null()
        for i, j in stensor_components(N):
            v = t[i, j]
            r.m[i][j] = v if v.is_atomic() else self.define(
                f"{prefix}{stensor_index(N, i, j)}", v)
            r.m[j][i] = r.m[i][j]
        return r
