#!/usr/bin/env python3
import sys
sys.path.append('./py')

from numbers import Number

from vector import Vector
from multivector import MultiVector

class GA:
    """GA(n) -> the rank-n geometric algebra with floating point coefficients

    Example:
    >>> ga = GA(3)
    >>> x = 3 + ga[1] + 2*ga[2] - 2*ga[3] + ga[1,2] - ga[2,3] + 2.5*ga[1,3] - ga[1,2,3]
    >>> x
    3.0 + 1.0*GA(3)[1] + 2.0*GA(3)[2] + 1.0*GA(3)[1,2] + -2.0*GA(3)[3] + 2.5*GA(3)[1,3] + -1.0*GA(3)[2,3] + -1.0*GA(3)[1,2,3]
    >>> x*(~x)*ga.I
    -15.0*GA(3)[1,2] + -19.0*GA(3)[1,3] + 2.0*GA(3)[2,3] + 27.25*GA(3)[1,2,3]

    See the multivector.MultiVector class for operations on elements
    of this geometric algebra.  The quantities above produced by GA's
    factory methods are instances of a MultiVector class.

    """
    
    def __init__(self, n):
        self._dim = n

    @property
    def n(self):
        """Return the rank of this geometric algebra"""
        return self._dim

    def scalar(self, s=0.0):
        """Create a scalar in this geometric algebra"""
        x = MultiVector(self.n)
        x[0] = float(s)
        return x

    def coerce_multivector(self, m):
        """Attempt to coerce a multivector from another geometric algebra into
this one.  It will fail if any nonzero terms have rank higher than
this algebra's rank.

        """
        x = MultiVector(self.n)
        for i in range(2**m.dim):
            if m[i]:
                x[i] = m[i]
        return x

    def blade(self, coef, *indices):
        """Create a blade in this algebra.  See the __getitem__ method for a
better interface."""
        x = self.scalar(coef)
        for i in indices:
            y = MultiVector(self.n)
            y[2**(i-1)] = 1.0
            x *= y
        return x

    @property
    def I(self):
        """Return the standard pseudoscalar for this algebra."""
        x = MultiVector(self.n)
        x[-1] = 1.0
        return x
    
    def __getitem__(self, index):
        """Create a blade for this algebra.  The index can be an integer or a
        tuple of integers.

        Example:
        >>> ga = GA(3)
        >>> 3*ga[1,2] #3e_{12} in another notation
        3.0*GA(3)[1,2]
        >>> 2*ga[1,2] - ga[1,2,3] + 5 #a multivector from summing blades and a scalar
        5.0 + 2.0*GA(3)[1,2] + -1.0*GA(3)[1,2,3]

        """
        if isinstance(index, tuple):
            return self.blade(1.0, *index)
        return self.blade(1.0, index)

    @classmethod
    def mv_from_vector(cls, v):
        """
        Create a multivector from a vector.

        >>> v = Vector(3, 1, 0, -1)
        >>> x = GA.mv_from_vector(v)
        >>> x
        1.0*GA(3)[1] + -1.0*GA(3)[3]
        """

        x = MultiVector(len(v))
        bit = 1
        i = 1
        while bit < len(x._data):
            x[bit] = v[i]
            i += 1
            bit *= 2
        return x

    @classmethod
    def vector_from_mv(cls, x):
        """Create a vector from a multivector, raising an exception if it has
        non-rank-one terms.

        >>> x = GA(3)[1] - GA(3)[3]
        >>> v = GA.vector_from_mv(x)
        >>> str(v)
        '(1.0, 0.0, -1.0)'

        """        
        n = x.dim
        l = [0.0]*n
        bit = 1
        i = 0
        while bit < len(x._data):
            l[i] = x[bit]
            i += 1
            bit *= 2
        return Vector(n, l)    
        

if __name__ == '__main__':
    import doctest
    doctest.testmod()    
