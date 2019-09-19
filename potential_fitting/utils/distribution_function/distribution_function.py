import math


class DistributionFunction(object):

    def get_value(self, x):
        raise NotImplementedError


class LinearDistributionFunction(DistributionFunction):

    def __init__(self, intercept, slope):
        self.intercept = intercept
        self.slope = slope

    def get_value(self, x):
        return self.intercept + x * self.slope

    @staticmethod
    def get_function_from_2_points(x1, y1, x2, y2):

        slope = (y2 - y1) / (x2 - x1)
        intercept = y1 - slope * x1

        return LinearDistributionFunction(intercept, slope)


class LogarithmicDistributionFunction(DistributionFunction):

    def __init__(self, min_val, max_val, min_x, max_x):
        self.min_val = min_val
        self.max_val = max_val
        self.min_x = min_x
        self.max_x = max_x

        self.dx = (math.log(max_val) - math.log(min_val)) / (max_x - min_x)

    def get_value(self, x):
        return math.e ** (math.log(self.min_val) + (x - self.min_x) * self.dx)

class ConstantDistributionFunction(DistributionFunction):

    def __init__(self, val):
        self.val = val

    def get_value(self, x):
        return self.val

class PiecewiseDistributionFunction(DistributionFunction):

    def __init__(self, functions, cutoffs):
        self.functions = functions
        self.cutoffs = cutoffs

        if not len(functions) + 1 != len(cutoffs):
            raise Error

    def get_value(self, x):
        for function, cutoff in zip(self.functions, self.cutoffs):
            if x < cutoff:
                return function.get_value(x)

        return self.functions[-1].get_value(x)