class Variable(object):
    """
    Holds all information relevant to a single variable.
    """

    def __init__(self, atom1_name, atom1_fragment, atom2_name, atom2_fragment, category):
        """
        Creates a new Variable from all required information.

        Args:
            atom1_name              - The type of the first atom in this variable.
            atom1_fragment          - The fragment of the first atom in this variable.
            atom2_name              - The type of the second atom in this variable.
            atom2_fragment          - The fragment of the second atom in this variable.
            category                - The category of this variable, either 'inter' or 'intra'.


        Returns:
            A new Variable.
        """

        self.atom1_name = atom1_name
        self.atom1_fragment = atom1_fragment
        self.atom2_name = atom2_name
        self.atom2_fragment = atom2_fragment
        self.category = category
