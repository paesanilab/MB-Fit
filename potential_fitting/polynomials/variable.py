class Variable(object):
    """
    Holds all information relevant to a single variable.
    """

    def __init__(self, atom1_name, atom1_index, atom1_fragment, atom2_name, atom2_index, atom2_fragment, category):
        """
        Creates a new Variable from all required information.

        Args:
            atom1_name              - The type of the first atom in this variable.
            atom1_index             - The index of the first atom in this variable.
                    Should be 1 if it is the first atom with its name, 2 if it is the second, etc.
            atom1_fragment          - The fragment of the first atom in this variable.
            atom2_name              - The type of the second atom in this variable.
            atom2_index             - The index of the second atom in this variable.
                    Should be 1 if it is the first atom with its name, 2 if it is the second, etc.
            atom2_fragment          - The fragment of the second atom in this variable.
            category                - The category of this variable, either 'inter' or 'intra'.


        Returns:
            A new Variable.
        """

        self.atom1_name = atom1_name + str(atom1_index)
        self.atom1_fragment = atom1_fragment
        self.atom2_name = atom2_name + str(atom2_index)
        self.atom2_fragment = atom2_fragment
        self.category = category
