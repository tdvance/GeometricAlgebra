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

e1 = ga[1]
e2 = ga[2]
e3 = ga[3]
assert e1.cross_product(e2) == e3
assert e2.cross_product(e3) == e1
assert e3.cross_product(e1) == e2

def random_scalar(n=3):
    if random.randint(0,100) < 33:
        return GA(n).scalar(random.randint(-5, 5))
    return GA(n).scalar(random.gauss(0, 20))

def random_vector(n=3):
    bit = 1
    x = GA(n).scalar(0)
    while bit < len(x):
        x[bit] = random_scalar(n)
        bit *= 2
    return x

def random_rank():
    return random.randint(0, 3)

def random_blade(n=3):
    x = GA(n).scalar(0)
    i = random.randint(0, 2**n - 1)
    x._data[i] = 1.0
    return x

def random_len():
    if random.randint(0,100) < 2:
        return 0
    return random.randint(1, 2**3 - 1)    

def random_multivector(n=3):
    x = GA(n).scalar(0)
    l = random_len()
    for i  in range(l):
        x += random_scalar(n)*random_blade(n)
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
    asserteq(x+y, y+x)
    asserteq(x+(y+z), (x+y)+z)
    asserteq(x+y - y, x)
    asserteq(2*x, x+x)
    asserteq(0*x, ga.scalar(0))
    asserteq(1*x, x)
    asserteq(-1*x, -x)
    asserteq(x-y, -(y-x))
    asserteq(x+0, x)
    asserteq(0+x, x)    
    a = random_scalar()
    asserteq(a*x, x*a)
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
    

def signum(l):
    s = 1
    for i in range(len(l)):
        for j in range(len(l)-1, i, -1):
            if l[j] < l[j-1]:
                x = l[j-1]
                l[j-1] = l[j]
                l[j] = x
                s = -s
    return s
                        

assert signum([1,2,3,6,12]) == 1
assert signum([1,2,6,3,12]) == -1

def dedup(l):
    x = []
    for v in l:
        if x and x[-1] == v:
            x.pop()
        else:
            x.append(v)
    return x

def blade_to_list(x):
    if not x:
        return []
    if len(x) > 1:
        return NotImplemented
    j = -1
    for i in x:
#        print('btol:',i, x[i])
        j = bin(i)[2:]
        l = [len(j) - i for i in range(len(j)) if j[i] == '1']
        l.reverse()
        return l
        
for iter in range(100):
    x = random_blade()
    y = random_blade()
    l = blade_to_list(x) + blade_to_list(y)
    s = signum(l)
    l = dedup(l)
    z = x*y
    assert blade_to_list(z) == l
    for i in z:
        assert z[i] == s
    
    a = random_blade()
    b = random_blade()
    c = random_blade()
    a1 = random_blade()
    b1 = random_blade()
    c1 = random_blade()

    x = a+2*b-c + 3
    y = 2*a -b - 2*c
    z = x*y
    zz = a*2*a - a*b - a*2*c +2*b*2*a - 2*b*b -2*b*2*c - c*2*a + c*b + c*2*c + 3*2*a - 3*b - 3*2*c
    assert z == zz

    #test MultiVector.__init__
    a = random_rank()
    x = MultiVector(a)
    assert x.dim == a
    assert not x
    assert len(x._data) == 2**a
    for i in x._data:
        assert not i

    #test MultiVector.__getitem__
    x = random_multivector(a)
    y = MultiVector(a)
#    print(x, a)
    for i in x:
        j = bin(i)[2:]
        l = [len(j) - i for i in range(len(j)) if j[i] == '1']
        l.reverse()
#        print(i, j, l, a)
        b = GA(a)[tuple(l)]
#        print(x[i], b)
        y += x[i] * b
#    print( )       
    asserteq(x, y)

    #test MultiVector.__setitem__
    x = random_multivector(a)
    y = MultiVector(a)
    for i in x:
        y[i] = x[i]
    asserteq(x, y)

    #test MultiVector.__delitem__
    x = random_multivector(a)
    y = MultiVector(a)
    j = -1
    for i in x:
        if j == -1:
            j = i
        else:
            y[i] = x[i]
    del x[j]
    asserteq(x, y)

    #test MultiVector.__contains__
    x = random_multivector(a)
    s = set(range(2**a))
    for i in x:
        assert i in x
        s.discard(i)
    for i in s:
        assert i not in x

    #test MultiVector.__iter__
    x = random_multivector(a)
    l = [i for i in x]
    for i in range(2**a):
        if x[i] != 0:
            assert l[0] == i
            l.pop(0)
        else:
            assert i not in l

    #test MultiVector.__len__
    x = random_multivector(a)
    l = [i for i in x]
    assert len(x) == len(l)
    
    #test MultiVector.__eq__
    x = random_multivector(a)    
    y = MultiVector(a)
    for i in x:
        y[i] = x[i]
    asserteq(x, y)
    z = random_blade(a)
    if not z:
        z = 1
    y += z
    assert x != y
    
    #test MultiVector.__pos__
    x = random_multivector(a)        
    assert x == +x
    assert x is not +x
    
    #test MultiVector.__neg__
    x = random_multivector(a)
    assert x == --x
    assert not (x + -x)
    
    #test MultiVector.__float__
    x = random_scalar(a)
    assert float(x) == x[0]

    #test MultiVector.__add__
    x = random_multivector(a)
    y = random_multivector(a)
    z = MultiVector(a)
    for i in range(0, 2**a):
        z[i] = x[i] + y[i]
    asserteq(x+y, z)
    y = random_scalar(a)
    z = +x
    z[0] += y[0]
    asserteq(x+y[0], z)
    asserteq(x+y, z)    
    
    
    #test MultiVector.__radd__
    x = random_multivector(a)
    y = random_multivector(a)
    z = MultiVector(a)
    for i in range(0, 2**a):
        z[i] = x[i] + y[i]
    asserteq(x+y, z)
    y = random_scalar(a)
    z = +x
    z[0] += y[0]
    asserteq(y[0]+x, z)
    asserteq(y+x, z)
    
    #test MultiVector.__sub__
    x = random_multivector(a)
    y = random_multivector(a)
    z = MultiVector(a)
    for i in range(0, 2**a):
        z[i] = x[i] - y[i]
    asserteq(x-y, z)
    y = random_scalar(a)
    z = +x
    z[0] -= y[0]
    asserteq(x-y[0], z)
    asserteq(x-y, z)    
    
    #test MultiVector.__rsub__
    x = random_multivector(a)
    y = random_multivector(a)
    z = MultiVector(a)
    for i in range(0, 2**a):
        z[i] = x[i] - y[i]
    asserteq(x-y, z)
    y = random_scalar(a)
    z = -x
    z[0] = y[0] + z[0]
    asserteq(y[0] - x, z)
    asserteq(y - x, z)
    
    #test MultiVector.__mul__
    x = random_multivector(a)
    y = random_multivector(a)
    z = MultiVector(a)
    for i in x:
        ii = []
        bit=1
        k=1
        while bit <= i:
            if bit & i:
                ii.append(k)
            k += 1
            bit *=2
        for j in y:            
            jj = []
            bit = 1
            k=1
            while bit <= j:
                if bit & j:
                    jj.append(k)
                k += 1
                bit *=2
            l = ii + jj
            s = signum(l)
            l = dedup(l)
            k = 0
            for kk in l:
                k += 2**(kk-1)
            z[k] += s * x[i] * y[j]
    asserteq(x*y, z)        

    y = random_scalar(a)
    z = x*y
    assert(x*y[0] == z)
    #test MultiVector.__rmul__
    x = random_multivector(a)
    y = random_scalar(a)
    z = x*y
    assert(y[0]*x == z)
    
    #test MultiVector.__and__
    x = random_multivector(a)    
    y = random_multivector(a)
    asserteq((x&y) + (x@y), x*y)
    asserteq(x&y, -(y&x))
    
    #test MultiVector.__matmul__
    x = random_multivector(a)    
    y = random_scalar(a)
    asserteq(x*y, x @ y)
    y = random_multivector(a)
    asserteq(x@y, y@x)
    
    #test MultiVector.__or__ TODO
    
    #test MultiVector.__abs__
    x = random_vector(a)
    asserteq(abs(x)*abs(x), x @ x)

    #test MultiVector.__invert__ TODO

    #test MultiVector.rank
    x = random_vector(a)
    assert not x or x.rank() == 1
    assert not x or (x*x).rank() == 2
    assert not (x+x*x) or (x+x*x).rank() == 2
    assert not x or (x*x*x).rank() == 3    
    
    #test MultiVector.cross_product
    x = random_vector(a)
    y = random_vector(a)
    z = x.cross_product(y)
    asserteq(z @ x, 0)
    asserteq(z @ y, 0)
    assert not z or z.rank() == 1
    
    #test MultiVector.left_inv
    x = random_multivector(a)    
    try:
        z = x.left_inv()
        if z is not NotImplemented:
            asserteq(x*z, one)
    except(TypeError, ZeroDivisionError):
        pass
    
    #test MultiVector.right_inv
    x = random_multivector(a)    
    try:
        z = x.right_inv()
        if z is not NotImplemented:
            asserteq(z*x, one)
    except(TypeError, ZeroDivisionError):
        pass
    
    #test MultiVector.__truediv__
    x = random_multivector(a)    
    y = random_multivector(a)
    try:
        z = x/y
        if z is not NotImplemented:
            asserteq(z*y, x)
    except(TypeError, ZeroDivisionError):
        pass
    
    #test MultiVector.dual TODO

    #test MultiVector.I
    x = random_multivector(a)
    assert x.I == GA(a).I
    assert x.I*x.I*x.I*x.I == GA(a).scalar(1)
    
    #test MultiVector.__str__ TODO

    #test MultiVector.__repr__ TODO



    
