import itertools
import numpy as np

N_VARS = 3
WANTED_DEG = 2
"""
def mult_comb(vals, idxs):
    #assuming vals is only 0-1 values.
    return all([vals[i] for i in idxs])

def f1(vals):
    return sum(vals) - sum([mult_comb(vals, idxs) for idxs in list(itertools.combinations(range(3),2))])

def f2(vals):
    return sum(vals) - sum([mult_comb(vals, idxs) for idxs in list(itertools.combinations(range(6),2))]) + \
        + sum([mult_comb(vals, idxs) for idxs in [[0,2,3],[0,1,4],[0,3,4],[1,2,3],[1,2,4],[0,1,5],[0,2,5],[1,3,5],[2,4,5],[3,4,5]]])
"""

ZEROS = np.zeros(2**N_VARS).reshape(*[2]*N_VARS)

def monom_sign(vals, idxs):
    #assuming vals is only +-1 values
    res = 1
    for x, idx in zip(vals, idxs):
        if idx != 0:
            res *= x
    return res

def get_poly_deg(poly):
    for deg in xrange(poly.ndim+1):
        poly_not_triv = (poly != 0) 
        if ((poly_not_triv & DEGREE_MASKS[deg]) == poly_not_triv).all():
            return deg

def is_participating(poly, idx):
    poly_not_triv = (poly != 0) 
    return not ((poly_not_triv & PARTICIPATION_MASKS[idx]) == ZEROS).all()

def create_participation_masks(n_vars = N_VARS):
    global PARTICIPATION_MASKS
    PARTICIPATION_MASKS = [np.zeros(2 ** n_vars).reshape(*([2] * n_vars)).astype(int) for i in xrange(n_vars)]
    for idxs in itertools.product(*[[0,1]]*n_vars):
        for idx in xrange(n_vars):
            PARTICIPATION_MASKS[idx][idxs] = idxs[idx]

def create_degree_masks(n_vars = N_VARS):
    global DEGREE_MASKS
    DEGREE_MASKS = [np.zeros(2 ** n_vars).reshape(*([2] * n_vars)).astype(int) for i in xrange(n_vars+1)]
    for idxs in itertools.product(*[[0,1]]*n_vars):
        pop_cnt = idxs.count(1)
        for deg in xrange(pop_cnt, n_vars + 1):
            DEGREE_MASKS[deg][idxs] = 1

def table2polynom(table):
    #this performs table2polynom, {+-1}**n --> +-1
    #what is the difference between rfftn and fftn?
    return (np.fft.rfftn(table) / 2 ** table.ndim).real

def apply_poly(poly, input_vector):
    assert poly.ndim == len(input_vector)
    res = 0
    for idxs in itertools.product(*[[0,1]]*poly.ndim):
        poly_coeff = poly.item(sum([(2**(len(idxs)-i-1)) * idxs[i] for i in xrange(len(idxs))]))
        res += poly_coeff * monom_sign(input_vector, idxs)
    return res

create_degree_masks()
create_participation_masks()

"""
for i,x in enumerate(PARTICIPATION_MASKS):
    print i
    print x
"""
"""
t = np.array([1,1,1,1]*2).reshape(2,2,2)
print get_poly_deg(table2polynom(t))
print [is_participating(table2polynom(t),i) for i in xrange(N_VARS)]
t = np.array([1,-1,1,-1]*2).reshape(2,2,2)
print get_poly_deg(table2polynom(t))
print [is_participating(table2polynom(t),i) for i in xrange(N_VARS)]
t = np.array([1,1,1,-1]*2).reshape(2,2,2)
print get_poly_deg(table2polynom(t))
print [is_participating(table2polynom(t),i) for i in xrange(N_VARS)]
t = np.array([1,1,1,1,1,1,1,-1]).reshape(2,2,2)
print get_poly_deg(table2polynom(t))
print [is_participating(table2polynom(t),i) for i in xrange(N_VARS)]
"""
#print [t[i][j][k] for i,j,k in itertools.product(*[[0,1]]*3)]
#poly = table2polynom(t)
#for inp in itertools.product(*[[1,-1]]*3):
    #print inp, apply_poly(poly, np.array(inp))
    #pass

build_arrays = [[1,-1]] * 2 ** N_VARS
build_arrays[0] = [1]
for i in xrange(N_VARS):
        build_arrays[1 << i] = [-1]
print build_arrays
print 'moo'
count_polys = 0
tqdm_counter = 0
for idxs in itertools.product(*build_arrays):
    tqdm_counter +=1
    if not tqdm_counter & 0xffff:
        print hex(tqdm_counter)
    #print 'idxs', idxs
    t = np.array(idxs).reshape(*[2] * N_VARS)
    poly = table2polynom(t)
    if get_poly_deg(poly) <= WANTED_DEG and all([is_participating(poly, idx) for idx in xrange(N_VARS)]):
        print "found poly!"
        print idxs,
        print poly
        count_polys += 1

print "count_polys", count_polys



