#!/usr/bin/env python3

import sys
sys.path.append('./py')

from collections import Iterable
import math
import operator

class Vector:
    """Vectors are immutable.  Coordinates are 1-up, in keeping with
mathematical usage though against common computer programming usage.

    Vector(n_int) -> zero vector in n-dimensional space.

    Vector(n_int, c1_float, c2_float, ..., cn_float) -> vector in
    n-dimensional space specified by coordinates.  Missing coordinates
    get a default value of 0.0.  Extra coordinates are an error.

    Vector(n_int, iterable) -> vector in n-dimensional space whose
    coordinates are specified by an iterable.  Missing coordinates get
    a default value of 0.0.  Extra coordinates are an error.

    >>> z3 = Vector(3)
    >>> bool(z3)
    False
    >>> len(z3)
    3
    >>> z3.as_tuple()
    (0.0, 0.0, 0.0)

    >>> v = Vector(3, 1, 0, 1)
    >>> bool(v)
    True
    >>> len(v)
    3
    >>> v.as_tuple()
    (1.0, 0.0, 1.0)

    >>> w = Vector(6, [1, 0, 1])
    >>> bool(w)
    True
    >>> len(w)
    6
    >>> w.as_tuple()
    (1.0, 0.0, 1.0, 0.0, 0.0, 0.0)
    """

    def __init__(self, n, *args):
        if len(args) == 1 and isinstance(args[0] , Iterable):
            l = [float(x) for x in args[0]]
        else:
            l = [float(x) for x in args]
        if len(l) > n:
            raise IndexError(len(l))
        if len(l) < n:
            l += [0.0]*(n - len(l))
        self._coefs = tuple(l)
            
    
    @property
    def dim(self):
        """x.dim() == len(x) is the dimensionality of the vector space this
        vector lives in, or equivalently, the total number of
        coordinate positions in the vector.
        
        >>> v = Vector(7)
        >>> v.dim
        7
        >>> w = Vector(3, 1, 0, 1)
        >>> w.dim
        3

        """
        return len(self)

    @property
    def x(self):
        """
        v.x -> v[1]

        >>> v = Vector(3, 1, 0, 1)
        >>> v.x
        1.0

        """
        return self[1]

    @property    
    def y(self):
        """
        v.y -> v[2]

        >>> v = Vector(3, 1, 0, 1)
        >>> v.y
        0.0


        """
        return self[2]

    @property
    def z(self):
        """
        v.z -> v[3]

        >>> v = Vector(3, 1, 0, 1)
        >>> v.z
        1.0

        """
        return  self[3]

    def __getitem__(self, i):
        """v[i] is the ith coordinate coefficient of v.  It is an error if i
        is not an integer, i<1, or i>len(v).

        >>> v = Vector(7, [1, 2, 2.5, 3, 1, -1, 0])
        >>> v[3]
        2.5

        """
        return self._coefs[operator.index(i) - 1]

    def __contains__(self, a):
        """
        a_float in v -> True if v[i] == a for some i.

        >>> v = Vector(7, [1, 2, 2.5, 3, 1, -1, 0])
        >>> 2.5 in v
        True
        >>> 5 in v
        False
        """
        return a in self._coefs
        

    def __iter__(self):
        """
        iter(v) -> generates the sequence of coefficients of v in order.

        >>> v = Vector(7, [1, 2, 2.5, 3, 1, -1, 0])
        >>> [x for x in v]
        [1.0, 2.0, 2.5, 3.0, 1.0, -1.0, 0.0]
        
        """
        return iter(self._coefs)

    def __len__(self):
        """x.dim() == len(x) is the dimensionality of the vector space this
        vector lives in, or equivalently, the total number of
        coordinate positions in the vector.
        
        >>> v = Vector(7)
        >>> len(v)
        7
        >>> w = Vector(3, 1, 0, 1)
        >>> len(w)
        3

        """
        
        return len(self._coefs)

    def __eq__(self,  other):
        """v == w if and only if they have the same dimension and
        corresponding coefficients are equal.

        >>> v = Vector(3, 1, 0, 1)
        >>> w = Vector(3, [1.0, 0.0, 1.0])
        >>> v == w
        True
        >>> v != w
        False
        >>> v is w
        False
        >>> x = Vector(4, 1, 0, 1)
        >>> v == x
        False
        >>> v != x
        True
        >>> y = Vector(3, 1, 0, -1)
        >>> v == y
        False
        >>> v != y
        True

        """
        return self._coefs == other._coefs

    def __pos__(self):
        """+v is v

        >>> v = Vector(3, 1, 0, 1)
        >>> +v == v
        True
        >>> +v is v
        True
        """
        return self

    def __neg__(self):
        """-v -> the vector whose coordinates are the negatives of
correwsponding coordinates of v.

        >>> v = Vector(3, 1, 0, 1)
        >>> w = -v
        >>> len(w)
        3
        >>> w.as_tuple() #-0.0 is a valid float value
        (-1.0, -0.0, -1.0)

        """
        return Vector(len(self), [-x for x in self])

    def __add__(self, other):
        """v + w is the vector having the same dimension as both v and w, and
whose coefficients are the sums of the coefficients at corresponding
coordinates in v and w.
        
        >>> v = Vector(3, 1, 0, 1)
        >>> w = Vector(3, 2, 1, -1)
        >>> (v + w).as_tuple()
        (3.0, 1.0, 0.0)

        """
        if len(self) != len(other):
            raise NotImplementedError("dimensions %d != %d" %(len(self), len(other)))
        return Vector(len(self), [x + y for (x,y) in zip(self, other)])

    def __sub__(self, other):
        """v - w is the vector having the same dimension as both v and w, and
whose coefficients are the differences of the coefficients at corresponding
coordinates in v and w.
        
        >>> v = Vector(3, 1, 0, 1)
        >>> w = Vector(3, 2, 1, -1)
        >>> (v - w).as_tuple()
        (-1.0, -1.0, 2.0)
        """
        if len(self) != len(other):
            raise NotImplementedError("dimensions %d != %d" %(len(self), len(other)))
        return Vector(len(self), [x - y for (x,y) in zip(self, other)])

    def __mul__(self, scalar):
        """v*a is the vector having the same dimension as v and whose
coefficients are the product of v's coefficients in the corresponding
coordinates, and the float value a.

        >>> v = Vector(3, 1, 0, 1)       
        >>> (v * 2.5).as_tuple()
        (2.5, 0.0, 2.5)

        """
        return Vector(len(self), [x* scalar for x in self])

    def __truediv__(self, scalar):
        """v/a is the vector having the same dimension as v and whose
coefficients are the product of v's coefficients in the corresponding
coordinates, and the float value 1.0/a.

        >>> v = Vector(3, 1, 0, 1)       
        >>> (v / 2.5).as_tuple()
        (0.4, 0.0, 0.4)

        """
        return Vector(len(self), [x / scalar for x in self])
    
    def __rmul__(self, scalar):
        """a*v is the vector having the same dimension as v and whose
coefficients are the product of v's coefficients in the corresponding
coordinates, and the float value a.
        >>> v = Vector(3, 1, 0, 1)       
        >>> (2.5*v).as_tuple()
        (2.5, 0.0, 2.5)
"""
        return Vector(len(self), [x* scalar for x in self])        

    def __matmul__(self, other):
        """v @ w is the dot product of v and w, the real number that is the
sum of the products of corresponding coefficients of v and w.

        >>> v = Vector(3, 1, 0, 1)       
        >>> w = Vector(3, -1, 2, 3)       
        >>> v @ w
        2.0

        """
        if len(self) != len(other):
            raise NotImplementedError("dimensions %d != %d" %(len(self), len(other)))
        return sum(x*y for (x,y) in zip(self, other))
        
    def __abs__(self):
        """abs(v) is the norm of v, the square root of the some of the squares
        of the coefficients of v.  It is also equal to math.sqrt(v @
        v).

        >>> v = Vector(3, 3, 4, 12)
        >>> abs(v)
        13.0
        """
        x = sum(x*x for x in self)
        return math.sqrt(x)

    def cross_product(self, other):
        """v.cross_product(w) is only defined if len(v)==len(w)==3.  It is
then the geometric cross product of the vectors, producing the unique
vector u such that u @ v == 0, u @ w == 0, u is perpendicular to the
plane of v and w (or 0 if v and w are not linearly independent) and
pointing in the direction that a right-threaded screw would advance if
the head is turned in the direction from v to w, if v and w are
positioned with their origin points at the same place, and such that
abs(u) is equal to the area of the parallelogram defined by v and
w.
        >>> v = Vector(3, 1, 0, 0)
        >>> w = Vector(3, 0, 1, 0)
        >>> v.cross_product(w)
        Vector(3, 0.0, 0.0, 1.0)

        """
        if len(self) != 3:
            raise NotImplementedError("dimension %d != 3" %(len(self), 3))
        if len(other) != 3:
            raise NotImplementedError("dimension %d != 3" %(len(other), 3))
        return Vector(3, self.y*other.z - self.z*other.y, self.z*other.x - self.x*other.z, self.x*other.y - self.y*other.x)

    def to_cylindrical(self):
        """
        Return (r, theta, z), the cylindrical coordinates of this 3d vector.

        >>> v = Vector(3, 1, 1, 0)
        >>> v.to_cylindrical()
        (1.4142135623730951, 0.7853981633974483, 0.0)
        """
        theta = math.atan2(self.y, self.x)
        r = math.sqrt(self.x*self.x + self.y*self.y)
        return (r, theta, self.z)

    def to_spherical(self):
        """Return (r, phi, theta), the spherical coordinates of this 3d
        vector, where the north pole is phi=0, and the south pole is
        phi=math.pi.

        >>> v = Vector(3, 1, 1, 0)
        >>> v.to_spherical()
        (1.4142135623730951, 1.5707963267948966, 0.7853981633974483)

        """        
        rho = math.sqrt(self.x*self.x + self.y*self.y + self.z*self.z)
        theta = math.atan2(self.y, self.x)
        if rho:
            phi = math.acos(self.z/rho)
        else:
            phi = 0
        return (rho, phi, theta)
    
    def __str__(self):
        return str(self.as_tuple())

    def __repr__(self):
        return 'Vector(' + str(len(self)) + ', ' + ', '.join([str(x) for x in self.as_tuple()]) + ')'

    def concat(self, other):
        """v.concat(w) is the vector u such that len(u)==len(v)+len(w) and the
coordinates of u, in order, are the coordinates of v, followed by the
coordinates of w.
        >>> v = Vector(3, 1, 0, 0)
        >>> w = Vector(3, 0, 1, 0)
        >>> str(v.concat(w))
        '(1.0, 0.0, 0.0, 0.0, 1.0, 0.0)'
        """
        return Vector(len(self) + len(other), [x for x in self] + [x for x in other])

    def as_tuple(self):
        """v.as_tuple() is the unique tuple t such that len(v)==len(t) and for
all i where this makes sense, t[i-1] == v[i].  Vector(len(t), *t) == v
if v.as_tuple() == t.
        >>> v = Vector(3, 1, 0, 0)
        >>> v.as_tuple()
        (1.0, 0.0, 0.0)

        """
        return self._coefs
    
    def __bool__(self):
        """bool(v) is True if v has any nonzero coefficients, false if it is a
zero vector.
        >>> v = Vector(3, 1, 0, 0)
        >>> bool(v)
        True
        >>> w = Vector(3)        
        >>> bool(w)
        False
        """
        return any(bool(x) for x in self)

    def normalize(self):
        """v.normalize() is v itself if v is a zero vector (that is,
abs(v)==0).  Otherwise, it is the unique vector u such that abs(u)==1
and there is a positive scalar float a such that a*v==u.
        >>> v = Vector(3, 3, 4, 12)
        >>> v.normalize() == Vector(3, 3/13, 4/13, 12/13)
        True
        """
        a = abs(self)
        if a:
            return self * (1.0/a)
        else:
            return self

    def angle(self, other):
        """v.angle(w) -> the angle in radians between v and w, at least 0 but
less than pi.
        >>> import math
        >>> v = Vector(3, 1, 0, 0)
        >>> w = Vector(3, 1, 1, 0)
        >>> abs(v.angle(w) - math.pi/4) < 1E-10
        True

        """
        return math.acos((self @ other) * (1.0/abs(self)/abs(other)))

    def E(self, i):
        """v.E(i) is the vector having the same dimension as v, but with 0.0
        in every coordinate position except for i (1 <= i <= len(v))
        having a 1.0 in that position.

        >>> v = Vector(7)
        >>> str(v.E(3))
        '(0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0)'

        """
        return Vector(len(self), [0]*(i-1) + [1])

    @classmethod
    def from_spherical(cls, rho, phi, theta):
        """
        Return (r, theta, z), the spherical coordinates of this 3d vector.
        >>> Vector.from_spherical(1.4142135623730951, 1.5707963267948966, 0.7853981633974483)
        Vector(3, 1.0000000000000002, 1.0, 1.0000000000000002)
        """        
        s = math.sin(phi)
        x = rho*s*math.cos(theta)
        y = rho*s*math.sin(theta)
        z = rho*math.cos(theta)
        return Vector(3, x, y, z)

    @classmethod
    def from_cylindrical(cls, r, theta, z):
        """
        Return (r, theta, z), the cylindrical coordinates of this 3d vector.
        >>> Vector.from_cylindrical(1.4142135623730951 , 0.7853981633974483, 0)
        Vector(3, 1.0000000000000002, 1.0, 0.0)
        """                
        x = r*math.cos(theta)
        y = r*math.sin(theta)
        return Vector(3, x, y, z)
    
if __name__ == '__main__':
    import doctest
    doctest.testmod()
