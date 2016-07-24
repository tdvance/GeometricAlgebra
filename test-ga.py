#!/usr/bin/env python3
from ga import GA
from multivector import MultiVector

ga = GA(3)
assert ga.n == 3

x = 4 + 3*ga[1] + 4*ga[2] + 5*ga[1,2]
y = 3 + 2*ga[1] + 3*ga[2] + 4*ga[1,2]
z = 10 + 16*ga[1] + 26*ga[2] + 32*ga[1,2]

assert x*y == z

ga2 = GA(2)
xx =ga2.coerce_multivector(x)

assert x != xx
assert xx == ga2.scalar(4) + ga2.blade(3, 1) + ga2.blade(4, 2) + ga2.blade(5, 1, 2)

