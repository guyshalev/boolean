
class Level(object):
    """
    represents a level of a boolean function. (points and values)
    This only helps us in the generation of monotone function.
    """
    def __init__(self, n, k, points_to_values):
        self._n = n
        self._k = k
        self._points_to_values = points_to_values

    def get_points(self):
        return self._points_to_values.keys()

    def set_value_for_point(self, p, val):
        self._points_to_values[p] = val

    def get_value_for_point(self, p):
        return self._points_to_values[p]

    def __str__(self):
        res_str = "N=" + str(n) + ", " + "K=" + str(k) + ":" + "\n"
        res_str.append(''.join(str(point) for point in self.points))
        return res_str