import pickle
import itertools
import numpy
import math
import lagrange_interpolation
from my_utils import *

N_VARS = 12
WANTED_DEG = 5

def apply_uni_poly(poly_coefs, x):
	#poly_coeffs is from high coef to low, so go [::-1]
	return sum([coef * x ** i for i,coef in enumerate(poly_coefs[::-1])])

def poly_change_input(from_format, to_format):
	pass
def poly_change_output(from_format, to_format):
	pass

def almost_eq(x,y, precision=10**-6):
	return abs(x-y) < precision

def compute_num_multivar_polys(poly, n_vars, deg):
	#compute how many polynoms we can create from the given symmetric poly
	#>>> poly = pickle.load(file('./sym_polys_3_6.pickle','rb'))[0]
	#>>> compute_num_multivar_polys(poly, 6, 3)
	# result should be ncr(20,10) because: [(1, 1), (6, 0), (15, 0), (20, 10), (15, 15), (6, 6), (1, 0)]
	inv_quants = [(ncr(n_vars, i)) for i in xrange(n_vars + 1)]
	vals = [int(round(((b+1)/2.) * inv_quants[a])) for a,b in poly]
	num_opts = [ncr(iq,val) for iq, val in zip(inv_quants, vals)]
	print zip(inv_quants, vals)
	print math.log(reduce(lambda x,y: x*y, num_opts),2), "bits enum"

def generate_legit_polys(deg=WANTED_DEG, n_vars=N_VARS):
	number_from_above = (deg - 1) if deg <= 3 else ((deg + 1) / 2)
	number_from_below = deg + 1 - number_from_above
	vals_to_enumerate = range(number_from_below) + [n_vars-i for i in range(number_from_above)[::-1]]
	##print "vals_to_enumerate", vals_to_enumerate
	#TODO - replace the interpolation with scipy.interpolate.lagrange
	#poly = [round(x, 10) for x in numpy.polyfit(xs, ys, deg)]
	inv_quants = [(ncr(n_vars, i)) for i in xrange(n_vars + 1)]
	quants = [1. / inv_quants[i] for i in xrange(len(inv_quants))]

	y_vals_options = [map(lambda x: 1 - 2*x, [x * quants[i] for x in xrange(inv_quants[i] + 1)]) for i in vals_to_enumerate]
	y_vals_options[0] = [1] #to ensure full_sensitivity
	y_vals_options[1] = [-1] #to ensure full_sensitivity
	#print y_vals_options
	
	polys = []
	counter = 0
	print "searching for symmetric polys...", n_vars, deg
	for y_vals in itertools.product(*y_vals_options):
		counter += 1
		#print "vals", vals_to_enumerate, y_vals
		poly_func = lagrange_interpolation.lagrange(zip(vals_to_enumerate, y_vals))
		#######poly = numpy.polyfit(vals_to_enumerate, y_vals, len(vals_to_enumerate) - 1)
		#######poly_func = lambda x: apply_uni_poly(poly, x)

		full_poly = [(x, poly_func(x)) for x in xrange(n_vars + 1)]
		#print "full_poly", full_poly
		if all([abs(y) <= 1.01 and almost_eq(inv_quants[i]*y, round(inv_quants[i]*y)) for i,y in full_poly]):
			#print "legit symmetric poly!"
			polys.append(full_poly)
	print "enumerated over %d polys" % counter
	return polys

def generate_legit_polys_sym_or_antisym(deg=WANTED_DEG, n_vars=N_VARS):
	number_from_above = (deg - 1) if deg <= 3 else ((deg + 1) / 2)
	number_from_below = max(deg + 1 - number_from_above, number_from_above)
	vals_to_enumerate = range(number_from_below)
	vals_to_mirror = [n_vars-i for i in range(number_from_above)[::-1]]
	print vals_to_enumerate
        
	inv_quants = [(ncr(n_vars, i)) for i in xrange(n_vars + 1)]
	quants = [1. / inv_quants[i] for i in xrange(len(inv_quants))]

	y_vals_options = [map(lambda x: 1 - 2*x, [x * quants[i] for x in xrange(inv_quants[i] + 1)]) for i in vals_to_enumerate]
	y_vals_options[0] = [1] #to ensure full_sensitivity
	y_vals_options[1] = [-1] #to ensure full_sensitivity
	#print y_vals_options
	
	polys = []
	counter = 0
	print "searching for symmetric polys...", n_vars, deg
	for y_vals in itertools.product(*y_vals_options):
		for sign in xrange(2): #determine if sym or anti sim
			counter += 1
			#print "vals", vals_to_enumerate, y_vals
			inp = zip(vals_to_enumerate, y_vals) + zip(vals_to_mirror, [x * (-1)**sign for x in y_vals[::-1]])
			poly_func = lagrange_interpolation.lagrange(inp)

			full_poly = [(x, poly_func(x)) for x in xrange(n_vars + 1)]
			#print "full_poly", full_poly
			if all([abs(y) <= 1.01 and almost_eq(inv_quants[i]*y, round(inv_quants[i]*y)) for i,y in full_poly]):
				print "legit symmetric poly!"
				polys.append(full_poly)
	print "enumerated over %d polys" % counter
	return polys

possible_parameters = {}
best_parameters = (5,15) #we know by now of [many] (28,7) that is sym, and there is no (29,7) sym/antisym.
#weird that I didn't find 21/6 along the way... look for precision bugs, maybe really work  over Q (rationals)
for n_vars in xrange(30):
	print n_vars
	for deg in xrange(2, n_vars):
		if math.log(deg,n_vars) > math.log(best_parameters[0],best_parameters[1]):
			continue
		#legit_polys = generate_legit_polys(deg, n_vars)
		legit_polys = generate_legit_polys_sym_or_antisym(deg, n_vars)
		if len(legit_polys) > 0:
			possible_parameters[n_vars] = (deg, len(legit_polys), legit_polys[0])
			best_parameters = (deg, n_vars)
			print "best_parameters = (%d, %d)" % (deg, n_vars)
			print possible_parameters[n_vars]
			pickle.dump(legit_polys, file("sym_polys_%d_%d.pickle" % (deg, n_vars),'wb'))
	#out(legit_polys)

print best_parameters
