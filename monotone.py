"""
This code supports functions that are of the form {0,1} ** n ==> {+1, -1}
This has two reasons:
* {0, 1} ** n allows us to use the binary representation which is convenient
* {+1, -1} is better when we want to work with fourier. Note that "+1 < -1" for that matter!
"""

import numpy
import itertools
import random
from utils import popcount, indices_and_n_to_num, num_and_n_to_indices


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


def get_upper_neighbors(point, all_points):
    """gets a Point, returns all neighbors from one level above
    all_points is a dictionary from num to point"""
    assert point in all_points
    res = []
    for i in xrange(point.n):
        if point.num ^ 2 ** (point.n - 1 - i) > point.num: # i'th from left is off
            neigh = all_points[point.num ^ 2 ** (point.n - 1 - i)]
            res.append(neigh)
    return res


def get_lower_neighbors(point, all_points):
    """gets a Point, returns all neighbors from one level below
    all_points is a dictionary from num to point"""
    assert point in all_points
    res = []
    for i in xrange(point.n):
        if point.num ^ 2 ** (point.n - 1 - i) < point.num: # i'th from left is on
            neigh = all_points[point.num ^ 2 ** (point.n - 1 - i)]
            res.append(neigh)
    return res


def get_all_upper_neighbours(point, all_points):
    """gets a Point, returns all points that are above it
    all_points is a dictionary from num to point"""
    return get_all_directional_neighbours(point, all_points, direction="UPPER")


def get_all_lower_neighbours(point, all_points):
    """gets a Point, returns all points that are below it
    all_points is a dictionary from num to point"""
    return get_all_directional_neighbours(point, all_points, direction="LOWER")


def get_all_directional_neighbours(point, all_points, direction):
    # If we have k zeros in the point, we expect 2**k - 1 points (because we exclude the point itself).
    one_indices = num_and_n_to_indices(point.num, point.n)
    zero_indices = [x for x in xrange(n) if x not in one_indices]
    indices_to_flip = []
    if direction == "UPPER":
        indices_to_flip = zero_indices
    elif direction == "LOWER":
        indices_to_flip = one_indices

    res = set()
    for num_to_flip in xrange(1, len(indices_to_flip) + 1):
        for indices in itertools.combinations(indices_to_flip, num_to_flip):
            num_diff = indices_and_n_to_num(indices, point.n)
            new_num = point.num ^ num_diff
            return res.add(all_points[new_num])

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
            upper_neighs = get_upper_neighbors(point, ctx.NUM_TO_POINT)
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


def sample_monotone_function(ctx, max_pq = 0.1):
    # This raises a natural question - how to randomly sample a monotone function?
    pass


def sample_random_function(ctx):
    points_to_values = {p: 1 - 2 * random.randint(0, 1) for p in ctx.NUM_TO_POINT.values()}
    return BooleanFunction(ctx, points_to_values)


def generate_majority_function(ctx):
    assert ctx.N % 2 == 1
    new_points_to_values = {p: popcount(p) for num, p in ctx.NUM_TO_POINT.iteritems()}
    return BooleanFunction(ctx, new_points_to_values)


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
        if isinstance(item, Point):
            return self.NUM_TO_POINT[num]
        elif isinstance(item, int):
            return self.LEVEL_TO_POINTS[item]
        else:
            raise ValueError()


if __name__ == "__main__":
    print "moo main"
    N = 4
    ctx = GlobalContext(N=N)
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

