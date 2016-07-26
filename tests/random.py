from ga import GA
import random

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

def random_blade(n=3):
    x = GA(n).scalar(0)
    i = random.randint(0, 2**n - 1)
    x._data[i] = 1.0
    return x

def random_rank(n=3):
    return random.randint(0, n)

def random_len(n=3):
    if random.randint(0,100) < 2:
        return 0
    return random.randint(1, 2**n - 1)    

def random_multivector(n=3):
    x = GA(n).scalar(0)
    l = random_len()
    for i  in range(l):
        x += random_scalar(n)*random_blade(n)
    return x
