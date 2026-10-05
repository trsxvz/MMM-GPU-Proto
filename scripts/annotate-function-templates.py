#!/usr/bin/env python3
"""Add TFEL_HOST_DEVICE to the function templates of the given headers.

A declaration or definition is annotated if it follows a single line
`template <...>` and starts with a return type followed by a function name
and an opening parenthesis. Lines already annotated, classes, aliases and
friend declarations are left untouched.
"""
import re
import sys

decl = re.compile(
    r'^(\s*)((?:\[\[nodiscard\]\]\s*)?(?:inline\s+|constexpr\s+|static\s+)*'
    r'(?:const\s+)?[A-Za-z_][A-Za-z0-9_:<>,\s]*?[&*]?\s+[A-Za-z_][A-Za-z0-9_:]*\s*\()(.*)$')
for p in sys.argv[1:]:
    lines = open(p).read().split('\n')
    n = 0
    for i in range(1, len(lines)):
        l = lines[i]
        prev = lines[i - 1].strip()
        if 'TFEL_HOST_DEVICE' in l or not prev.startswith('template'):
            continue
        if not prev.endswith('>'):
            continue
        s = l.strip()
        if s.startswith(('struct', 'class', 'using', 'typedef', 'return', 'friend', '//', 'requires')):
            continue
        m = decl.match(l)
        if m:
            lines[i] = m.group(1) + 'TFEL_HOST_DEVICE ' + m.group(2) + m.group(3)
            n += 1
    if n:
        open(p, 'w').write('\n'.join(lines))
    print(p, n)
