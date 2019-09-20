import math

from potential_fitting.exceptions import InconsistentValueError


class DistributionFunction(object):
    """
    Abstract class that implements standard function in the x-y coordinate plane.
    """

    def get_value(self, x):
        """
        Get the y value of this function at the given x.

        Args:
            x   - The x coordinate to compute the value of y at.

        Returns:
            The value of y at the given x coordinate.
        """

        raise NotImplementedError


class LinearDistributionFunction(DistributionFunction):
    """
    Implementation of DistributionFunction that describes a linear function.
    """

    def __init__(self, intercept, slope):
        """
        Constructs a new LinearDistributionFunction.

        A LinearDistributionFunction describes a line.

        The value returned by self.get_value(x) will be intercept + x * slope

        Args:
            intercept   - The y-intercept of the line.
            slope       - The slope of the line. The change in y relative to the change in x.

        Returns:
            A new LinearDistributionFunction.
        """

        self.intercept = intercept
        self.slope = slope

    def get_value(self, x):
        """
        Get the y value of this function at the given x.

        Args:
            x   - The x coordinate to compute the value of y at.

        Returns:
            The value of y at the given x coordinate.
        """

        return self.intercept + x * self.slope

    @staticmethod
    def get_function_from_2_points(x1, y1, x2, y2):
        """
        Constructs a new LinearDistributionFunction from 2 points
        on a line.

        Args:
            x1  - The x coordinate of the first point.
            y1  - The y coordinate of the first point.
            x2  - The x coordinate of the second point.
            y2  - The y coordinate of the second point.

        Returns:
            A new LinearDistributionFunction.
        """

        slope = (y2 - y1) / (x2 - x1)
        intercept = y1 - slope * x1

        return LinearDistributionFunction(intercept, slope)


class GeometricDistributionFunction(DistributionFunction):
    """
    Implementation of DistributionFunction that describes a geometric function.
    """

    def __init__(self, coefficient, base):
        """
        Constructs a new GeometricDistributionFunction.

        A GeometricDistributionFunction describes a geometric function.

        The value returned by self.get_value(x) will be coefficient * base ** x

        Args:
            coefficient   - The coefficient as given in the equation above.
            base       - The base of the exponent, as given in the equation above.

        Returns:
            A new GeometricDistributionFunction.
        """

        self.coefficient = coefficient
        self.base = base

    def get_value(self, x):
        """
        Get the y value of this function at the given x.

        Args:
            x   - The x coordinate to compute the value of y at.

        Returns:
            The value of y at the given x coordinate.
        """

        return self.coefficient * self.base ** x


class LogarithmicDistributionFunction(DistributionFunction):
    """
    Implementation of DistributionFunction that describes a logarithmic distribution.
    """

    def __init__(self, min_val, max_val, min_x, max_x):
        """
        Constructs a new LogarithmicDistributionFunction.

        A LogarithmicDistributionFunction describes a logarithmic function.

        The value returned by self.get_value(x) will be
        e ** (log(min_val) + (x - min_x) * (log(max_val) - log(min_val)) / (max_x - min_x))

        Args:
            min_val     - self.get_value(min_x) will return min_val.
            max_val     - self.get_value(max_x) will return max_val.
            min_x       - self.get_value(min_x) will return min_val.
            max_x       - self.get_value(max_x) will return max_val.

        Returns:
            A new LogarithmicDistributionFunction.
        """

        self.min_val = min_val
        self.max_val = max_val
        self.min_x = min_x
        self.max_x = max_x

        self.dx = (math.log(max_val) - math.log(min_val)) / (max_x - min_x)

    def get_value(self, x):
        """
        Get the y value of this function at the given x.

        Args:
            x   - The x coordinate to compute the value of y at.

        Returns:
            The value of y at the given x coordinate.
        """

        return math.e ** (math.log(self.min_val) + (x - self.min_x) * self.dx)


class ConstantDistributionFunction(DistributionFunction):
    """
    Implementation of DistributionFunction that describes a constant distribution.
    """

    def __init__(self, val):
        """
        Constructs a new ConstantDistributionFunction.

        A ConstantDistributionFunction describes a constant function.

        The value returned by self.get_value(x) will be val

        Args:
            val - The constant value of the distribution.

        Returns:
            A new ConstantDistributionFunction.
        """

        self.val = val

    def get_value(self, x):
        """
        Get the y value of this function at the given x.

        Args:
            x   - The x coordinate to compute the value of y at.

        Returns:
            The value of y at the given x coordinate.
        """

        return self.val


class PiecewiseDistributionFunction(DistributionFunction):
    """
    Implementation of DistributionFunction that describes a piecewise distribution.
    """

    def __init__(self, functions, cutoffs):
        """
        Constructs a new PiecewiseDistributionFunction.

        A PiecewiseDistributionFunction describes a function that behaves
        differently depending on where in the domain we look.

        The value returned by self.get_value(x) will func.get_value(x) where func is
        functions[-1] if x > cutoffs[-1],
        otherwise functions[i] where i is the smallest integer where x < cutoffs[i]

        Args:
            functions   - A list of DistributionFunctions. The functions in this piecewise distribution.
                    Ordered from the lowest domain to the highest domain.
            cutoffs     - The cutoffs between the domains of the functions. There should be one less
                    cutoff than function.

        Returns:
            A new PiecewiseDistributionFunction.
        """

        self.functions = functions
        self.cutoffs = cutoffs

        if not len(functions) + 1 != len(cutoffs):
            raise InconsistentValueError("len(functions)", "len(cutoffs)", len(functions), len(cutoffs),
                                         "there should be one more function that cutoff.")

    def get_value(self, x):
        """
        Get the y value of this function at the given x.

        Args:
            x   - The x coordinate to compute the value of y at.

        Returns:
            The value of y at the given x coordinate.
        """

        for function, cutoff in zip(self.functions, self.cutoffs):
            if x < cutoff:
                return function.get_value(x)

        return self.functions[-1].get_value(x)


class RandomDistributionFunction(DistributionFunction):
    """
    Implementation of DistributionFunction that randomizes another DistributionFunction.
    """

    def __init__(self, function, random, min, max):
        """
        Constructs a new RandomDistributionFunction.

        A RandomDistributionFunction causes the get_value function to use a
        random x value instead of the passed in x value.

        The value returned by self.get_value(x) will be
        function.get_value(x1) where x1 is a random number in the range (min, max)

        Args:
            function    - A DistributionFunction to randomize.
            random      - The Random object to use to generate random x values.
            min         - The lower bound of random x values.
            max         - The upper bound of random x values.

        Returns:
            A new RandomDistributionFunction.
        """

        self.function = function
        self.random = random
        self.min = min
        self.max = max

    def get_value(self, x):
        """
        Get the y value of this function at the given x.

        Args:
            x   - The x coordinate to compute the value of y at.

        Returns:
            The value of y at the given x coordinate.
        """
        return self.function.get_value(self.random.uniform(self.min, self.max))
