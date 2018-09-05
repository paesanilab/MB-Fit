class Filter()

    def keep(monomial, variables)
        raise NotImplementedError

class MasterFilter()

    def __init__():
        self.subfilters = []

    def keep(monomial, variables):
        for filter in self.subfilters:
            if filter.keep():
                return True
        return False

class MaxDegreeFilter():
    def __init__(max_degree):
        self.max_degree = max_degree

    def(
