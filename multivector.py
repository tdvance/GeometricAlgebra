#!/usr/bin/env python3

import math
from numbers import Number


class MultiVector:
    """
    MultiVector(dim) -> the (mutable) zero multivector in GA(dim), for dim >= 0
    """
    def __init__(self, dim):
        self._data = [0.0]*(2**dim)
        self._dim = dim


    @property
    def dim(self):
        """
        x.dim() -> the order of the geometric algebra containing this multivector.
        """
        return self._dim

    def __getitem__(self, i):
        """Get the ith coefficient, where the bits of i illustrate the
standard basis blade to query.  For example, 11 = 8+2+1 = 2**3 + 2**1
+ 2**0, so x[11] returns the coefficient of [4,2,1] in x.  The index i
must be from 0 to 2**dim - 1.

        This is low-level.  It is better to use the arithmetic
        operations and the GA class to work with multivectors.

        """
        return self._data[i]

    def __setitem__(self, i, coef):
        """Set the ith coefficient, where the bits of i illustrate the
standard basis blade to set.  For example, 11 = 8+2+1 = 2**3 + 2**1 +
2**0, so x[11]=1.2 sets the coefficient of [4,2,1] in x to 1.2,
overwriting any coefficient already there.  The index i must be from 0
to 2**dim - 1.

        This is low-level.  It is better to use the arithmetic
        operations and the GA class to work with multivectors.

        """        
        self._data[i] = float(coef)

    def __delitem__(self, i):
        """Delete (set to zero) the ith coefficient, where the bits of i
illustrate the standard basis blade to set.  For example, 11 = 8+2+1 =
2**3 + 2**1 + 2**0, so del x[11] sets the coefficient of [4,2,1] in x
to 0, overwriting any coefficient already there.  The index i must be
from 0 to 2**dim - 1.

        This is low-level.  It is better to use the arithmetic
        operations and the GA class to work with multivectors.

        """                
        self._data[i] = 0.0

    def __contains__(self, i):
        """Return true if the ith coefficient is nonzero, where the bits of i
illustrate the standard basis blade to query.  For example, 11 = 8+2+1
= 2**3 + 2**1 + 2**0, so 11 in x is true if the coefficient of [4,2,1]
in x is nonzero.  The index i must be from 0 to 2**dim - 1.

        This is low-level.  It is better to use the arithmetic
        operations and the GA class to work with multivectors.

        """                        
        return bool(self._data[i])
        
    def __iter__(self):
        """Generate the indices i having nonzero coefficients, where the bits
of i illustrate the standard basis blade to query.  For example, 11 =
8+2+1 = 2**3 + 2**1 + 2**0, so 11 is among the indices generated if
the coefficient of [4,2,1] in x is nonzero.  

        This is low-level.  It is better to use the arithmetic
        operations and the GA class to work with multivectors.

        """                                
        for i in range(len(self._data)):
            if self._data[i]:
                yield i

    def __len__(self):
        """Return the number of nonzero terms in this multivector."""
        return len(self._data) - self._data.count(0.0)

    def __eq__(self, other):
        """Return True if this multivector is equal to the other multivector,
having equal coefficents and the same dimensionality for the parent
geometric algebra.

        """
        return isinstance(other, MultiVector) and self._data == other._data

    def __pos__(self):
        """
        Return a copy of this multivector.
        """
        x = MultiVector(self.dim)
        x._data = self._data.copy()
        return x

    def __neg__(self):
        """
        Return the negation of this multivector.
        """
        x = MultiVector(self.dim)
        for i in range(len(self._data)):
            x._data[i] = -self._data[i]
        return x

    def __float__(self):
        """Raise a TypeError exception if this multivector is not a scalar.
        Otherwise, return the floating point coercian of the scalar.

        """
        i = len(self)
        if i==0:
            return 0.0
        if i==1:
            return self[0]
        raise TypeError(self)

    def __add__(self, other):
        """
        Add the two multivectors, or a multivector to a number.
        """
        if isinstance(other, Number):
            x = +self
            x[0] += float(other)
            return x
        m = max(self.dim, other.dim)
        x = MultiVector(m)
        for i in range(len(self._data)):
            x._data[i] = self._data[i]
        for i in range(len(other._data)):
            x._data[i] += other._data[i]
        return x

    def __radd__(self, other):
        return self + other
    
    def __sub__(self, other):
        """
        Subtract the two multivectors, or a multivector and a number.
        """
        return self + (-other)

    def __rsub__(self, other):
        return -self + other
    
    def __mul__(self, other):
        """
        Multiply the two multivectors, or a multivector and a number.
        """
        if isinstance(other, Number):
            x = MultiVector(self.dim)
            for i in range(len(self._data)):
                x._data[i] = self._data[i] * float(other)
            return x
        m = max(self.dim, other.dim)
        x = MultiVector(m)
        for i in range(len(self._data)):
            if self._data[i]:
                for j in range(len(other._data)):
                    if other._data[j]:
                        k, s = MultiVector._blade_combine(i, j)
                        x._data[k] += self._data[i]*other._data[j]*s
        return x

    def __rmul__(self, other):
        if isinstance(other, Number):
            return self*other
        if isinstance(other, MultiVector):
            return other*self
        return NotImplemented
    
    def __and__(self, other):
        """Find the outer (wedge) product of the two multivectors, or of a multivector
        and a number.

        """        
        x = self*other - other*self
        for i in range(len(x._data)):
            x._data[i] /= 2
        return x

    def __matmul__(self, other):
        """Find the inner (dot) product of the two multivectors, or of a multivector
        and a number.

        """                
        x = self*other + other*self
        for i in range(len(x._data)):
            x._data[i] /= 2
        return x

    def __or__(self, other):
        """Find the meet (vee) of the two multivectors, or of a multivector
        and a number.

        """                
        return self.dual() * other
    
    def __abs__(self):
        """Find the norm of the multivector"""
        x = 0.0
        for i in self:
            a = self._data[i]
            x += a*a
        return math.sqrt(x)

    def __invert__(self):
        """Return the reversion of the multivector"""
        x = MultiVector(self.dim)
        for i in range(len(self._data)):
            v = self._data[i]
            if v:
                r = MultiVector._rank(i) % 4
                if r==2 or r==3:
                    x[i] = -v
                else:
                    x[i] = v
        return x

    def rank(self):
        """Compute the rank (maximum grade among terms) of the multivector"""
        r = -math.inf
        for i in range(len(self._data)):
            if self._data[i]:
                r = max(r, MultiVector._rank(i))
        return r

    def cross_product(self, other):
        """
        Return the cross product of the two multivectors.
        """
        return (self & other).dual()

    def left_inv(self):
        """Return the left inverse of the multivector if it exists.
        """
        try:
            x = ~self
            s = 1.0 / float(self*x)
        except TypeError:
            return NotImplemented
        for i in range(len(x._data)):
            x._data[i] *= s
        return x

    def right_inv(self):
        """Return the right inverse of the multivector if it exists.
        """        
        try:
            x = ~self
            s = 1.0 / float(x*self)
        except TypeError:
            return NotImplemented
        for i in range(len(x._data)):
            x._data[i] *= s
        return x    
                
    def __truediv__(self, other):
        """Divide the multivectors, if possible.  Or divide the multivector by the scalar.
        """
        if isinstance(other, Number):
            return self * (1.0 / other)
        x = other.left_inv()
        if x == NotImplemented:
            return x
        return self * x
    
    def dual(self):
        """
        Return the dual of the multivector.
        """
        x = self.I
        return self*x*x*x

    @property
    def I(self):
        """
        Return the standard pseudoscalar of the algebra this multivector is in.
        """
        x = MultiVector(self.dim)
        x._data[-1] = 1.0
        return x

    def __str__(self):
        result = ""
        for i in self:
            coef = self[i]
            if result:
                result += ' + '
            result += str(coef)
            bit = 1
            k = 1
            first = True
            while bit <= i:
                if i & bit:
                    if first:
                        result += '*['
                        first = False
                    else:
                        result += ','
                    result += '%d' % k
                k += 1
                bit *= 2
            if not first:
                result += ']'
        if not result:
            return '0'
        return result

    def __repr__(self):
        result = ""
        for i in self:
            coef = self[i]
            if result:
                result += ' + '
            result += str(coef)
            bit = 1
            k = 1
            first = True
            while bit <= i:
                if i & bit:
                    if first:
                        result += '*GA(%d)[' % self._dim
                        first = False
                    else:
                        result += ','
                    result += '%d' % k
                k += 1
                bit *= 2
            if not first:
                result += ']'
        if not result:
            return '0'
        return result

    @staticmethod
    def _rank(a):
        """
        Return the hamming weight of the integer a, viewed as a bitstring.

        >>> MultiVector._rank(0)
        0
        >>> MultiVector._rank(65536)
        1
        >>> MultiVector._rank(65535)
        16
        >>> MultiVector._rank(1+8+16+128+1024+16384 + 1048576)
        7
        """
        return bin(a).count('1')
    
        
    @staticmethod
    def _blade_combine(a, b):
        """Combine standard basis blades (represented by integers) by
        multiplying a and b, returning a tuple (c, s) where c is an
        integer standard basis blade, s is 1 or -1, such that s times
        the blade represented by c is the blade represented by a times
        b.

        >>> MultiVector._blade_combine(1+2+4, 2+4+8)
        (9, -1)
        >>> MultiVector._blade_combine(1, 1)
        (0, 1)
        >>> MultiVector._blade_combine(1, 2)
        (3, 1)
        >>> MultiVector._blade_combine(2,1)
        (3, -1)
        >>> MultiVector._blade_combine(16, 512)
        (528, 1)
        >>> MultiVector._blade_combine(512, 16)
        (528, -1)
        >>> MultiVector._blade_combine(127, 65535-127)
        (65535, 1)
        >>> MultiVector._blade_combine(65535-127, 127)
        (65535, -1)

        """
        if a==0:
            return b, 1
        if b==0:
            return a, 1
        c = a ^ b
        s = 1
        p = max(a,b)
        d = MultiVector._rank(a)
        e = 1
        while e <= p:
            if e & a:
                d -= 1
            if (d&1) and (e & b):
                s = -s
            e *= 2
        return c, s
                
        
if __name__ == '__main__':
    import doctest
    doctest.testmod()
