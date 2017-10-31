"""
This code supports functions that are of the form {0,1} ** n ==> {+1, -1}
This has two reasons:
* {0, 1} ** n allows us to use the binary representation which is convenient
* {+1, -1} is better when we want to work with fourier. Note that "+1 < -1" for that matter!
"""

import numpy
import itertools
import random
from utils import popcount, indices_and_n_to_num, num_and_n_to_indices, out


class Point(object):
    # TESTED
    def __init__(self, n, num):
        """for example: n=5, indices=[0,2] represents 10100"""
        self.n = n
        self.num = num
        self.indices = num_and_n_to_indices(num, n)
        self.binary_tuple = tuple(str(bin(self.num))[2:])
        self.popcount = popcount(num)

    def get_popcount(self):
        return self.popcount

    # TESTED
    def __neg__(self):
        return Point(self.n, 2 ** self.n - 1 - self.num)

    def __str__(self):
        return bin(self.num)[2:].zfill(self.n)

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        return self.n == other.n and self.num == other.num

    def __ne__(self, other):
        return not self.__eq__(other)


def get_upper_neighbours(point, ctx):
    """gets a Point, returns all neighbors from one level above
    all_points is a dictionary from num to point"""
    all_points = ctx.NUM_TO_POINT
    assert point in all_points.values()
    res = []
    for i in xrange(point.n):
        if point.num ^ 2 ** (point.n - 1 - i) > point.num: # i'th from left is off
            neigh = all_points[point.num ^ 2 ** (point.n - 1 - i)]
            res.append(neigh)
    return res


def get_lower_neighbours(point, ctx):
    """gets a Point, returns all neighbors from one level below
    all_points is a dictionary from num to point"""
    all_points = ctx.NUM_TO_POINT
    assert point in all_points.values()
    res = []
    for i in xrange(point.n):
        if point.num ^ 2 ** (point.n - 1 - i) < point.num: # i'th from left is on
            neigh = all_points[point.num ^ 2 ** (point.n - 1 - i)]
            res.append(neigh)
    return res


def get_all_upper_neighbours(point, ctx):
    """gets a Point, returns all points that are above it
    all_points is a dictionary from num to point"""
    return get_all_directional_neighbours(point, ctx, direction="UPPER")


def get_all_lower_neighbours(point, ctx):
    """gets a Point, returns all points that are below it
    all_points is a dictionary from num to point"""
    return get_all_directional_neighbours(point, ctx, direction="LOWER")


def get_all_directional_neighbours(point, ctx, direction):
    # If we have k zeros in the point, we expect 2**k points (because we include the point itself).
    all_points = ctx.NUM_TO_POINT
    one_indices = num_and_n_to_indices(point.num, point.n)
    zero_indices = [x for x in xrange(point.n) if x not in one_indices]
    indices_to_flip = []
    if direction == "UPPER":
        indices_to_flip = zero_indices
    elif direction == "LOWER":
        indices_to_flip = one_indices
    else:
        print "no direction!!!!!"
        exit(1)

    res = set()
    for num_to_flip in xrange(len(indices_to_flip) + 1):
        for indices in itertools.combinations(indices_to_flip, num_to_flip):
            num_diff = indices_and_n_to_num(indices, point.n)
            new_num = point.num ^ num_diff
            res.add(all_points[new_num])
    return res


class BooleanFunction(object):
    def __init__(self, ctx, points_to_values):
        self.ctx = ctx
        self.n = ctx.N
        self.points_to_values = points_to_values
        self.val_array = numpy.array([points_to_values[ctx.NUM_TO_POINT[num]] for num in xrange(2 ** self.n)])
        self.val_array = numpy.resize(self.val_array, [2] * self.n)
        self.fft_arr = self.calc_fft()

    def __str__(self, by_levels=True):
        res = []
        if by_levels:
            for level_num, level in sorted(ctx.LEVEL_TO_POINTS.items()):
                res.append("*" * 10 + (" %d " % level_num) + "*" * 10 + "\n")
                for point in level:
                    res.append(str(point) + "\t" + str(self[point]) + "\n")
        else:
            for num in xrange(2 ** self.n):
                point = ctx.NUM_TO_POINT[num]
                res.append(str(point) + "\t" + str(self[point]) + "\n")
        return "".join(res)

    def __setitem__(self, point, val):
        self.points_to_values[point] = val

    def __getitem__(self, point):
        return self.points_to_values[point]

    # TESTED
    def calc_fft(self):
        fft_arr = numpy.fft.fftn(self.val_array) / 2 ** self.n
        return map(float, fft_arr.reshape([2 ** self.n]))

    def get_fourier_coeff(self, indices):
        num = indices_and_n_to_num(indices, self.n)
        return self.fft_arr[num]

    def is_monotone(self):
        for point, val in self.points_to_values.iteritems():
            upper_neighs = get_upper_neighbours(point, ctx)
            for un in upper_neighs:
                if self[un] > val: #reverse because +1,-1
                    print point, val, un, self[un]
                    return False
        return True


# TESTED
def generate_dual_function(bool_func):
    new_points_to_values = {p: - bool_func[-p] for p in bool_func.points_to_values.keys()}
    return BooleanFunction(bool_func.ctx, new_points_to_values)


# TESTED
def generate_not_of_function(bool_func):
    new_points_to_values = {p: - bool_func[p] for p in bool_func.points_to_values.keys()}
    return BooleanFunction(bool_func.ctx, new_points_to_values)


# TESTED
def multiply_functions(bool_func_1, bool_func_2):
    new_points_to_values = {p: bool_func_1[p] * bool_func_2[p] for p in bool_func_1.points_to_values.keys()}
    return BooleanFunction(bool_func_1.ctx, new_points_to_values)


def sample_monotone_function(ctx):
    # This raises a natural question - how to randomly sample a monotone function?
    # for now, we will start with the zero function, randomize point, and flip them with some probability
    points_to_values = {p: 1 for p in ctx.NUM_TO_POINT.values()}
    n = ctx.N
    for num_flip in xrange(max(2 ** int(n / 2), 100)): # for now, magic numbers...
        num = random.randint(0, 2 ** n - 1)
        point = ctx.NUM_TO_POINT[num]
        level = popcount(num)
        dist_from_mid_level = 1 + abs(level - float(n) / 2)
        # TODO - maybe flip with respect to level (extreme levels don't flip
        to_flip = random.random() < 0.5 * (1 / dist_from_mid_level) ** 2
        if to_flip:
            if points_to_values[point] == 1:
                all_to_flip = get_all_upper_neighbours(point, ctx)
                for p_to_flip in all_to_flip:
                    points_to_values[p_to_flip] = -1
            else:
                all_to_flip = get_all_lower_neighbours(point, ctx)
                for p_to_flip in all_to_flip:
                    points_to_values[p_to_flip] = 1

    return BooleanFunction(ctx, points_to_values)

def sample_random_function(ctx):
    points_to_values = {p: 1 - 2 * random.randint(0, 1) for p in ctx.NUM_TO_POINT.values()}
    return BooleanFunction(ctx, points_to_values)


def generate_majority_function(ctx):
    assert ctx.N % 2 == 1
    new_points_to_values = {p: -1 if (float(p.popcount) / p.n > 0.5) else 1 for p in ctx.NUM_TO_POINT.values()}
    return BooleanFunction(ctx, new_points_to_values)


class Level(object):
    def __init__(self, ctx, points_to_values):
        self.ctx = ctx
        self.k = points_to_values.keys()[0].popcount
        self.rep_as_num = "".join([str(b) for p, b in sorted(points_to_values.items(), key=lambda (pt, b): pt.num)])
        self.points_to_values = points_to_values

    def __str__(self):
        res_str = ["N=" + str(ctx.N) + ", " + "K=" + str(self.k) + ":"]
        res_str.extend([str(p) + " " + str(b) for p, b in sorted(self.points_to_values.items(), key=lambda (pt, b): pt.num)])
        return '\n'.join(res_str)

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        return self.rep_as_num == other.rep_as_num

    def __ne__(self, other):
        return not self.__eq__(other)


# TESTED
def generate_all_levels(ctx):
    n = ctx.N
    levels = ctx.LEVEL_TO_POINTS
    res = [None] * (n + 1)
    for level_num, l in levels.iteritems():
        all_level_options = []
        for num in xrange(2 ** len(l)):
            indices = set(num_and_n_to_indices(num, len(l)))
            level_option = {}
            for idx_p, p in enumerate(l):
                level_option[p] = -1 if (idx_p in indices) else 1
            all_level_options.append(Level(ctx, level_option))
        res[level_num] = all_level_options
    return res


def levels_are_consistent(ctx, lower_level, upper_level):
    for point, val in lower_level.points_to_values.iteritems():
        if val == 1:
            continue #no harm can be done
        upper_neighbors = get_upper_neighbours(point, ctx)
        for upper_neigh in upper_neighbors:
            if upper_level.points_to_values[upper_neigh] == 1:
                return False #found a counterexample to monotonicity
    return True


def generate_all_legal_upgrades(ctx, all_levels):
    n = ctx.N
    level_option_to_legal_upgrades = {}
    for level_num in xrange(n):
        levels_lower = all_levels[level_num]
        levels_upper = all_levels[level_num + 1]
        for lower_level in levels_lower:
            lower_level_legal_upgrades = []
            for upper_level in levels_upper:
                if levels_are_consistent(ctx, lower_level, upper_level):
                    lower_level_legal_upgrades.append(upper_level)
            level_option_to_legal_upgrades[lower_level] = lower_level_legal_upgrades

    return level_option_to_legal_upgrades

def generate_all_monotones(ctx, all_levels, level_option_to_legal_upgrades):
    #build iteratively, level after level...
    n = ctx.N
    all_monotones = [[x] for x in all_levels[0]]
    for level in xrange(1, n + 1):
        new_all_monotones = []
        for partial_monotone in all_monotones:
            highest_level = partial_monotone[-1]
            legal_upgrades = level_option_to_legal_upgrades[highest_level]
            for legal_upgrade in legal_upgrades:
                new_all_monotones.append(partial_monotone[:] + [legal_upgrade])
        all_monotones = new_all_monotones
    # TODO - create them as BooleanFunction
    res = []
    for monotone in all_monotones:
        points_to_values = {}
        for level in monotone:
            points_to_values.update(level.points_to_values)
        res.append(BooleanFunction(ctx, points_to_values))
    return res






class GlobalContext(object):
    def __init__(self, N):
        self.N = N
        self.NUM_TO_POINT = {}
        self.LEVEL_TO_POINTS = {i: [] for i in xrange(N + 1)}
        for i in xrange(2 ** N):
            self.NUM_TO_POINT[i] = Point(N, i)
            curr_level = popcount(i)
            self.LEVEL_TO_POINTS[curr_level].append(self.NUM_TO_POINT[i])

        for l in self.LEVEL_TO_POINTS.values():
            l.sort(key=lambda p: p.num)

    def __getitem__(self, item):
        return self.NUM_TO_POINT[item]


if __name__ == "__main__":
    print "moo main"
    N = 3
    ctx = GlobalContext(N=N)
    all_levels = generate_all_levels(ctx)
    legal_upgrades = generate_all_legal_upgrades(ctx, all_levels)
    all_monotones = generate_all_monotones(ctx, all_levels, legal_upgrades)
    print len(all_monotones)
    print all_monotones[14]
    assert all([x.is_monotone for x in all_monotones])

    


    """
    a = get_all_upper_neighbours(point, ctx)
    b = get_all_lower_neighbours(point, ctx)
    c = get_upper_neighbours(point, ctx)
    d = get_lower_neighbours(point, ctx)
    mon_funcs = [sample_monotone_function(ctx) for i in xrange(100)]
    rand_funcs = [sample_random_function(ctx) for i in xrange(1000)]
    print sum([bf.is_monotone() for bf in mon_funcs])
    for x in [bf for bf in rand_funcs if bf.is_monotone()]:
        print x
    print "what"
    mon_func = sample_monotone_function(ctx)
    #print mon_func

    maj_bf = BooleanFunction(ctx, {Point(N, num) : val for num, val in zip(range(2**N), 2 * [1, 1, 1, -1, 1, -1, -1, -1])})
    print maj_bf
    print map(float, maj_bf.fft_arr)
    not_maj_bf = generate_not_of_function(maj_bf)
    print not_maj_bf
    print map(float, not_maj_bf.fft_arr)

    mult_func = multiply_functions(maj_bf, not_maj_bf)
    print mult_func
    print mult_func.fft_arr

    print '*'*40
    random_func = sample_random_function(ctx)
    print random_func
    print random_func.fft_arr
    dual_func = generate_dual_function(random_func)
    print dual_func
    print dual_func.fft_arr
    """

