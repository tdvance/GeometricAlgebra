#!/usr/bin/env python3
from ga import GA
from multivector import MultiVector
import random

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

def random_scalar():
    if random.randint(0,100) < 33:
        return ga.scalar(random.randint(-5, 5))
    return ga.scalar(random.gauss(0, 20))

def random_rank():
    return random.randint(0, 3)

def random_blade():
    x = ga.scalar(0)
    i = random.randint(0, 2**3 - 1)
    x._data[i] = 1.0
    return x

def random_len():
    if random.randint(0,100) < 2:
        return 0
    return random.randint(1, 2**3 - 1)    

def random_multivector():
    x = ga.scalar(0)
    l = random_len()
    for i  in range(l):
        x += random_scalar()*random_blade()
    return x

def asserteq(a, b):
    if a!=b:
        if abs(a - b) > 1E-10:
            print('a=', a)
            print('b=', b)
            print('diff=', abs(a-b))
            assert(a==b)

i = ga.blade(1,2,3)
j = ga.blade(1,1,3)
k = ga.blade(1,1,2)

asserteq(i*i,-1)
asserteq(j*j,-1)
asserteq(k*k,-1)
asserteq(i*j*k,-1)

one = ga.scalar(1)

asserteq(ga.I.left_inv(), ga.I*ga.I*ga.I)
asserteq(ga.I.right_inv(), ga.I*ga.I*ga.I)
asserteq(i.left_inv(), -i)
asserteq(j.left_inv(), -j)
asserteq(k.left_inv(), -k)
asserteq(i.right_inv(), -i)
asserteq(j.right_inv(), -j)
asserteq(k.right_inv(), -k)        

asserteq(i*j, k)
asserteq(j*k, i)
asserteq(k*i, j)

for iter in range(100):
    x = random_multivector()
    y = random_multivector()
    z = random_multivector()
    asserteq((x*y)*z, x*(y*z))
    asserteq(x*y, x @ y + (x & y))
    asserteq(2*x, x+x)

    a = random_scalar()
    b = random_scalar()
    c = random_scalar()
    d = random_scalar()
    q = a+b*i+c*j+d*k
    r =  a-b*i-c*j-d*k
    asserteq(q*r, abs(q)**2)

    try:
        xx = x.left_inv()
        if xx is not NotImplemented:
            asserteq(xx*x, one)
    except(TypeError, ZeroDivisionError):
        pass

    try:
        xx = x.right_inv()
        if xx is not NotImplemented:
            asserteq(x*xx, one)
    except(TypeError, ZeroDivisionError):
        pass    
    
   
