class Filter(object):

    def keep(self, monomial, variables):
        raise NotImplementedError

class DegreeFilter(Filter):

    def __init__(self, variable_string, degree_string, term_string):
        self.variables = variable_string.split("/")

        self.degrees = degree_string.split("/")
        self.terms = term_string.split("/")

        print(self.variables, self.degrees, self.terms)

    def keep(self, monomial, variables):
        term = sum(monomial)

        applicable_term = False

        for term_string in self.terms:
            if term_string == "*":
                applicable_term = True
                break
            elif term_string.endswith("-"):
                if term <= int(term_string[:-1]):
                    applicable_term = True
                    break
            elif term_string.endswith("+"):
                if term >= int(term_string[:-1]):
                    applicable_term = True
                    break
            elif "-" in term_string:
                if term >= int(term_string[:term_string.index("-")]) and term <= int(term_string[term_string.index("-") + 1:]):
                    applicable_term = True
                    break
            else:
                if term == int(term_string):
                    applicable_term = True
                    break

        if not applicable_term:
            return True

        for degree, variable in zip(monomial, variables):

            applicable_variable = False

            for variable_string in self.variables:
                if variable_string == "*":
                    applicable_variable = True
                    break
                if "*" in variable_string:
                    fits_type = False
                    fits_atoms = False

                    if variable_string[-2] == "*" and variable_string[-1] == "*":
                        fits_atoms = True
                    elif variable_string[-2] == "*":
                        fits_atoms = variable.category[-1] == variable_string[-1] or variable.category[-2] == variable_string[-1]
                    elif variable_string[-1] == "*":
                        fits_atoms = variable.category[-1] == variable_string[-2] or variable.category[-2] == variable_string[-2]
                    else:
                        fits_atoms = variable.category[-2:] == variable_string[-2:]

                    if variable_string[2] == "*" and variable_string[3] != "*":
                        fits_type = True
                    else:
                        fits_type = variable_string[:-2] == variable.category[:-2]

                    if fits_type and fits_atoms:
                        applicable_variable = True
                        break 
                elif variable_string == variable.category:
                    applicable_variable = True
                    break
            

            if not applicable_variable:
                continue

            for degree_string in self.degrees:

                if degree_string == "*":
                    return False
                elif degree_string.endswith("-"):
                    if degree <= int(degree_string[:-1]):
                        return False
                elif degree_string.endswith("+"):
                    if degree >= int(degree_string[:-1]):
                        return False
                elif "-" in degree_string:
                    if degree >= int(degree_string[:degree_string.index("-")]) and degree <= int(degree_string[degree_string.index("-") + 1:]):
                        return False
                else:
                    if degree == int(degree_string):
                        return False

        return True


