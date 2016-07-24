from numbers import Number

from multivector import MultiVector

class GA:
    def __init__(self, n):
        self._dim = n

    @property
    def n(self):
        return self._dim

    def scalar(self, s=0.0):
        x = MultiVector(self.n)
        x[0] = float(s)
        return x

    def coerce_multivector(self, m):
        x = MultiVector(self.n)
        for i in range(2**m.dim):
            if m[i]:
                x[i] = m[i]
        return x

    def blade(self, coef, *indices):
        x = self.scalar(coef)
        for i in indices:
            y = MultiVector(self.n)
            y[2**(i-1)] = 1.0
            x *= y
        return x

    @property
    def I(self):
        x = MultiVector(self.n)
        x[-1] = 1.0
        return x
    
    def __getitem__(self, index):
        if isinstance(index, tuple):
            return self.blade(1.0, *index)
        return self.blade(1.0, index)
