from math import factorial


def out(l):
    if isinstance(l, dict):
        for i, x in l.iteritems():
            print str(i) + ":", x
    else:
        for i, x in enumerate(l):
            print str(i)+":", x


def mifkad(l):
    d = {}
    for x in l:
        if x not in d:
            d[x] = 1
        else:
            d[x] += 1
    return d


def product(iterable):
    prod = 1
    for n in iterable:
        prod *= n
    return prod


def npr(n, r):
    """
    Calculate the number of ordered permutations of r items taken from a population of size n.
    npr(3, 2) = 6
    """
    assert 0 <= r <= n
    return product(range(n - r + 1, n + 1))


def ncr(n, r):
    """
    Calculate the number of unordered combinations of r items taken from a population of size n.
    ncr(3, 2) = 3
    """
    assert 0 <= r <= n
    if r > n // 2:
        r = n - r
    return npr(n, r) // factorial(r)


def popcount(num): # TODO - improve performance? maybe just use a table...
    return bin(num).count("1")


def indices_and_n_to_num(indices, n):
    """(0,2), 5 ==> 0b10100 which is 20"""
    return sum([2 ** (n - 1 - i) for i in indices])


def num_and_n_to_indices(num, n):
    """20, 5 (0b10100) ==> (0,2)"""
    res = []
    for i in xrange(n):
        if num ^ 2 ** (n - 1 - i) < num: # i'th from left is lit
            res.append(i)
    return res