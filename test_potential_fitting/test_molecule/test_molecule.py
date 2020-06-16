import unittest, random, os

from test_potential_fitting.test_case_with_id import TestCaseWithId
from potential_fitting.utils import Quaternion
from potential_fitting.molecule import Atom
from potential_fitting.molecule import Fragment
from potential_fitting.molecule import Molecule
from potential_fitting.exceptions import InconsistentValueError

"""
Test cases for molecule class
"""
class TestMolecule(TestCaseWithId):
    def __init__(self, *args, **kwargs):
        super(TestMolecule, self).__init__(*args, **kwargs)
        self.test_folder = os.path.dirname(os.path.abspath(__file__))

    def test_constructor(self):

        with self.assertRaises(InconsistentValueError):
            Molecule([Fragment([Atom("H", "A", 0, 0, 0)], "frag1", 0, 1, "H"),
                      Fragment([Atom("O", "B", 0, 0, 0)], "frag1", 0, 1, "O")])

        with self.assertRaises(InconsistentValueError):
            Molecule([Fragment([Atom("H", "A", 0, 0, 0)], "frag1", 0, 1, "H"),
                      Fragment([Atom("H", "B", 0, 0, 0)], "frag1", 0, 1, "H")])

        with self.assertRaises(InconsistentValueError):
            Molecule([Fragment([Atom("H", "A", 0, 0, 0)], "frag1", 0, 1, "H"),
                      Fragment([Atom("H", "A", 0, 0, 0)], "frag2", 0, 1, "H")])

        self.test_passed = True

    def test_get_name(self):
        mol = Molecule([Fragment([], "frag1", 0, 1, "")])

        self.assertEqual(mol.get_name(), "frag1")

        mol = Molecule([Fragment([], "frag1", 0, 1, ""),
                        Fragment([], "frag1", 0, 1, "")])

        self.assertEqual(mol.get_name(), "frag1-frag1")

        mol = Molecule([Fragment([], "frag1", 0, 1, ""),
                        Fragment([], "frag2", 0, 1, "")])

        self.assertEqual(mol.get_name(), "frag1-frag2")

        mol = Molecule([Fragment([], "frag1", 0, 1, ""),
                        Fragment([], "frag2", 0, 1, ""),
                        Fragment([], "frag1", 0, 1, "")])

        self.assertEqual(mol.get_name(), "frag1-frag2-frag1")

        self.test_passed = True

    def test_get_symmetry(self):
        mol = Molecule([])

        self.assertEqual(mol.get_symmetry(), "")

        mol = Molecule([Fragment([Atom("H", "A", 0, 0, 0)], "frag1", 0, 1, "H")])

        self.assertEqual(mol.get_symmetry(), "A1")

        mol = Molecule([Fragment([Atom("H", "A", 0, 0, 0),
                                  Atom("H", "A", 1, 1, 1)], "frag1", 0, 1, "HH")])

        self.assertEqual(mol.get_symmetry(), "A2")

        mol = Molecule([Fragment([Atom("H", "A", 0, 0, 0),
                                  Atom("H", "B", 1, 1, 1)], "frag1", 0, 1, "HH")])

        self.assertEqual(mol.get_symmetry(), "A1B1")

        mol = Molecule([Fragment([Atom("H", "A", 0, 0, 0),
                                  Atom("H", "A", 1, 1, 1)], "frag1", 0, 1, "HH"),
                        Fragment([Atom("H", "A", 0, 0, 0),
                                  Atom("H", "A", 1, 1, 1)], "frag1", 0, 1, "HH")
                        ])

        self.assertEqual(mol.get_symmetry(), "A2_A2")

        mol = Molecule([Fragment([Atom("O", "A", 0, 0, 0),
                                  Atom("H", "B", 1, 1, 1),
                                  Atom("H", "B", 1, 1, 1)], "frag1", 0, 1, "O(H)H"),
                        Fragment([Atom("H", "C", 0, 0, 0),
                                  Atom("H", "D", 1, 1, 1)], "frag2", 0, 1, "HH")
                        ])

        self.assertEqual(mol.get_symmetry(), "A1B2_C1D1")

        self.test_passed = True

    def test_get_fragments(self):
        molecule = Molecule([])

        # get_fragments() should return list of len 0 before any fragments added to Molecule
        self.assertEqual(len(molecule.get_fragments()), 0)

        fragment0 = Fragment([Atom("H", "A", 0, 0, 0)], "H", -3, 2, "H")

        molecule = Molecule([fragment0])

        self.assertEqual(molecule.get_fragments()[0], fragment0)
        # get_fragments() should return list of len 1 after 1 fragment added to Molecule
        self.assertEqual(len(molecule.get_fragments()), 1)
        
        fragment1 = Fragment([Atom("H", "B", 0, 0, 0), Atom("He", "C", 234, -0.1, 32), Atom("He", "D", 34, 43.53, 0)], "HHe2", 1, 1, "H[He][He]")

        molecule = Molecule([fragment0, fragment1])

        self.assertEqual(molecule.get_fragments()[0], fragment0)
        self.assertEqual(molecule.get_fragments()[1], fragment1)

        # get_fragments() should return list of len 2 after 2 fragments added to Molecule
        self.assertEqual(len(molecule.get_fragments()), 2)

        self.test_passed = True

    def test_get_atoms(self):
        molecule = Molecule([])

        # get_atoms() should return list of len 0 before any atoms added to Molecule
        self.assertEqual(len(molecule.get_atoms()), 0)

        atom0 = Atom("H", "A", 0, 0, 0)

        fragment0 = Fragment([atom0], "H", -3, 2, "H")

        molecule = Molecule([fragment0])

        self.assertEqual(molecule.get_atoms()[0], atom0)
        # get_atoms() should return list of len 1 after 1 atom added to Molecule
        self.assertEqual(len(molecule.get_atoms()), 1)

        atom1 = Atom("H", "B", 0, 0, 0)
        atom2 = Atom("He", "C", 234, -0.1, 32)
        atom3 = Atom("He", "D", 34, 43.53, 0)

        fragment1 = Fragment([atom1, atom2, atom3], "HHe2", 0, 1, "H[He][He]")

        molecule = Molecule([fragment0, fragment1])

        self.assertEqual(molecule.get_atoms()[0], atom0)
        self.assertEqual(molecule.get_atoms()[1], atom1)
        self.assertEqual(molecule.get_atoms()[2], atom2)
        self.assertEqual(molecule.get_atoms()[3], atom3)
        # get_atoms() should return list of len 4 after 4 atoms added to molecule
        self.assertEqual(len(molecule.get_atoms()), 4)

        self.test_passed = True

    def test_get_charge(self):
        molecule = Molecule([])

        # get_charge() should return 0 before any fragments added to Molecule
        self.assertEqual(molecule.get_charge(), 0)

        molecule = Molecule([Fragment([], "Name", -1, 3, "")])

        # get_charge() should return -1 after first fragment added to Molecule
        self.assertEquals(molecule.get_charge(), -1)

        molecule = Molecule([Fragment([], "Name", -1, 3, ""), Fragment([], "Name", 3, 1, "")])

        # get_charge() should return 2 after second fragment added to Molecule
        self.assertEquals(molecule.get_charge(), 2)

        self.test_passed = True

    def test_get_spin_multiplicity(self):
        molecule = Molecule([])

        # get_spin_multiplicity() should return 1 before any fragments added to Molecule
        self.assertEqual(molecule.get_spin_multiplicity(), 1)

        molecule = Molecule([Fragment([], "Name", -1, 3, "")])

        # get_spin_multiplicity() should return 3 after first fragment added to Molecule
        self.assertEquals(molecule.get_spin_multiplicity(), 3)

        molecule = Molecule([Fragment([], "Name", -1, 3, ""), Fragment([], "Name", 3, 1, "")])

        # get_spin_multiplicity() should return 3 after second fragment added to Molecule
        self.assertEquals(molecule.get_spin_multiplicity(), 3)

        self.test_passed = True

    def test_get_num_fragments(self):
        molecule = Molecule([])
        
        # get_num_fragments() should return 0 before any fragments added to molecule
        self.assertEqual(molecule.get_num_fragments(), 0)

        molecule = Molecule([Fragment([Atom("F", "A", 0, 0, 0)], "F", 1, 2, "F")])

        # get_num_fragments() should return 1 after 1 fragment added to molecule
        self.assertEqual(molecule.get_num_fragments(), 1)

        molecule = Molecule([Fragment([Atom("F", "A", 0, 0, 0)], "F", 1, 2, "F"), Fragment([Atom("F", "B", 0, 0, 3), Atom("Cl", "C", 0, 0, 0)], "FCl", 3, 2, "F[Cl]")])

        # get_num_fragments() should return 2 after 2 fragment added to molecule
        self.assertEqual(molecule.get_num_fragments(), 2)

        self.test_passed = True

    def test_get_num_atoms(self):
        molecule = Molecule([])
        
        # get_num_atoms() should return 0 before any atoms added to molecule
        self.assertEqual(molecule.get_num_atoms(), 0)

        molecule = Molecule([Fragment([Atom("F", "A", 0, 0, 0)], "F", 1, 2, "F")])

        # get_num_fragments() should return 1 after 1 atom added to molecule
        self.assertEqual(molecule.get_num_atoms(), 1)

        molecule = Molecule([Fragment([Atom("F", "A", 0, 0, 0)], "F", 1, 2, "F"), Fragment([Atom("F", "B", 0, 0, 3), Atom("Cl", "C", 0, 0, 0)], "FCl", 3, 2, "F[Cl]")])

        # get_num_fragments() should return 3 after 3 atoms added to molecule
        self.assertEqual(molecule.get_num_atoms(), 3)

        self.test_passed = True

    def test_translate(self):

        for i in range(100):
            x1, y1, z1 = random.random() * 100, random.random() * 100, random.random() * 100
            x2, y2, z2 = random.random() * 100, random.random() * 100, random.random() * 100
            x3, y3, z3 = random.random() * 100, random.random() * 100, random.random() * 100
            x4, y4, z4 = random.random() * 100, random.random() * 100, random.random() * 100

            dx, dy, dz = random.random() * 100, random.random() * 100, random.random() * 100

            mol = Molecule([Fragment([Atom("O", "A", x1, y1, z1),
                                      Atom("H", "B", x2, y2, z2)], "OH-", -1, 1, "OH"),
                            Fragment([Atom("O", "A", x3, y3, z3),
                                      Atom("H", "B", x4, y4, z4)], "OH-", -1, 1, "OH")
                            ])

            ref_mol = Molecule([Fragment([Atom("O", "A", x1 + dx, y1 + dy, z1 + dz),
                                          Atom("H", "B", x2 + dx, y2 + dy, z2 + dz)], "OH-", -1, 1, "OH"),
                                Fragment([Atom("O", "A", x3 + dx, y3 + dy, z3 + dz),
                                          Atom("H", "B", x4 + dx, y4 + dy, z4 + dz)], "OH-", -1, 1, "OH")
                                ])

            mol.translate(dx, dy, dz)

            self.assertEqual(mol, ref_mol)

        self.test_passed = True

    def test_rotate(self):

        for i in range(100):

            ref_x, ref_y, ref_z = random.random() * 100, random.random() * 100, random.random() * 100

            ref_atom = Atom("O", "A", ref_x, ref_y, ref_z)

            x1, y1, z1 = random.random() * 100, random.random() * 100, random.random() * 100
            x2, y2, z2 = random.random() * 100, random.random() * 100, random.random() * 100
            x3, y3, z3 = random.random() * 100, random.random() * 100, random.random() * 100
            x4, y4, z4 = random.random() * 100, random.random() * 100, random.random() * 100

            mol = Molecule([Fragment([Atom("O", "A", x1, y1, z1),
                                      Atom("H", "B", x2, y2, z2)], "OH-", -1, 1, "OH"),
                            Fragment([Atom("O", "A", x3, y3, z3),
                                      Atom("H", "B", x4, y4, z4)], "OH-", -1, 1, "OH")
                            ])

            pre_dist0 = ref_atom.distance(mol.get_atoms()[0])
            pre_dist1 = ref_atom.distance(mol.get_atoms()[1])
            pre_dist2 = ref_atom.distance(mol.get_atoms()[2])
            pre_dist3 = ref_atom.distance(mol.get_atoms()[3])

            mol.rotate(Quaternion.get_random_rotation_quaternion(), origin_x=ref_x, origin_y=ref_y, origin_z=ref_z)

            post_dist0 = ref_atom.distance(mol.get_atoms()[0])
            post_dist1 = ref_atom.distance(mol.get_atoms()[1])
            post_dist2 = ref_atom.distance(mol.get_atoms()[2])
            post_dist3 = ref_atom.distance(mol.get_atoms()[3])

            self.assertAlmostEqual(pre_dist0, post_dist0)
            self.assertAlmostEqual(pre_dist1, post_dist1)
            self.assertAlmostEqual(pre_dist2, post_dist2)
            self.assertAlmostEqual(pre_dist3, post_dist3)

        self.test_passed = True

    def test_move_to_center_of_mass(self):

        mol = Molecule([Fragment([Atom("O", "A", 2, 0, 0),
                                  Atom("H", "B", 3, 0, 0)], "OH-", -1, 1, "OH"),
                        Fragment([Atom("O", "A", 5, 0, 0),
                                  Atom("H", "B", 4, 0, 0)], "OH-", -1, 1, "OH")
                        ])

        ref_mol = Molecule([Fragment([Atom("O", "A", -1.5, 0, 0),
                                      Atom("H", "B", -0.5, 0, 0)], "OH-", -1, 1, "OH"),
                            Fragment([Atom("O", "A", 1.5, 0, 0),
                                      Atom("H", "B", 0.5, 0, 0)], "OH-", -1, 1, "OH")
                            ])

        mol.move_to_center_of_mass()

        self.assertEqual(mol, ref_mol)

        mol = Molecule([Fragment([Atom("O", "A", 2, -2, 0),
                                  Atom("H", "B", 3, -3, 0)], "OH-", -1, 1, "OH"),
                        Fragment([Atom("O", "A", 5, -5, 0),
                                  Atom("H", "B", 4, -4, 0)], "OH-", -1, 1, "OH")
                        ])

        ref_mol = Molecule([Fragment([Atom("O", "A", -1.5, 1.5, 0),
                                      Atom("H", "B", -0.5, 0.5, 0)], "OH-", -1, 1, "OH"),
                            Fragment([Atom("O", "A", 1.5, -1.5, 0),
                                      Atom("H", "B", 0.5, -0.5, 0)], "OH-", -1, 1, "OH")
                            ])

        mol.move_to_center_of_mass()

        self.assertEqual(mol, ref_mol)

        self.test_passed = True

    def rotate_on_principal_axes(self):

        ref_mols = [Molecule([Fragment([Atom("O", "A", 0, 0, -2),
                                      Atom("O", "A", 0, 0, 2)], "OO", 0, 1, "OO")]),
                    Molecule([Fragment([Atom("O", "A", 0, 0, 2),
                                        Atom("O", "A", 0, 0, -2)], "OO", 0, 1, "OO")])
                    ]

        mol = Molecule([Fragment([Atom("O", "A", 2, 0, 0),
                                  Atom("O", "A", -2, 0, 0)], "OO", 0, 1, "OO")])

        mol.rotate_on_principal_axes()

        self.assertIn(mol, ref_mols)

        mol = Molecule([Fragment([Atom("O", "A", -2, 0, 0),
                                  Atom("O", "A", 2, 0, 0)], "OO", 0, 1, "OO")])

        mol.rotate_on_principal_axes()

        self.assertIn(mol, ref_mols)

        mol = Molecule([Fragment([Atom("O", "A", 1.41421356238, 1.41421356238, 0),
                                  Atom("O", "A", -1.41421356238, -1.41421356238, 0)], "OO", 0, 1, "OO")])

        mol.rotate_on_principal_axes()

        self.assertIn(mol, ref_mols)

        mol = Molecule([Fragment([Atom("O", "A", 1.15470053838, -1.15470053838, 1.15470053838),
                                  Atom("O", "A", -1.15470053838, 1.15470053838, -1.15470053838)], "OO", 0, 1, "OO")])

        mol.rotate_on_principal_axes()

        self.assertIn(mol, ref_mols)

        self.test_passed = True

    def test_get_excluded_pairs(self):
        mol = Molecule([Fragment([Atom("O", "A", 1.41421356238, 1.41421356238, 0),
                                  Atom("O", "A", -1.41421356238, -1.41421356238, 0)], "OO", 0, 1, "OO")])

        self.assertEqual(mol.get_excluded_pairs(), [[[[0, 1]]], [[]], [[]]])

        mol = Molecule([Fragment([Atom("O", "A", -1.5, 1.5, 0),
                                  Atom("H", "B", -0.5, 0.5, 0)], "OH-", -1, 1, "OH"),
                        Fragment([Atom("O", "A", 1.5, -1.5, 0),
                                  Atom("H", "B", 0.5, -0.5, 0)], "OH-", -1, 1, "OH")
                        ])

        self.assertEqual(mol.get_excluded_pairs(), [[[[0, 1]], [[0, 1]]], [[], []], [[], []]])

        mol = Molecule([Fragment([Atom("O", "A", -1.5, 1.5, 0),
                                  Atom("H", "B", -0.5, 0.5, 0),
                                  Atom("H", "B", -0.5, -0.5, 0)], "OH2", 0, 1, "O(H)H"),
                        Fragment([Atom("O", "A", 1.5, -1.5, 0),
                                  Atom("H", "B", 0.5, -0.5, 0),
                                  Atom("H", "B", 0.5, 0.5, 0)], "OH2", 0, 1, "O(H)H")
                        ])

        self.assertEqual(mol.get_excluded_pairs(), [[[[0, 1], [0, 2]], [[0, 1], [0, 2]]], [[[1, 2]], [[1, 2]]], [[], []]])

        mol = Molecule([Fragment([Atom("O", "A", -1.5, 1.5, 0),
                                  Atom("H", "B", -0.5, 0.5, 0),
                                  Atom("H", "B", -0.5, -0.5, 0)], "OH2", 0, 1, "O(H)H"),
                        Fragment([Atom("O", "A", 1.5, -1.5, 0),
                                  Atom("H", "B", 0.5, -0.5, 0),
                                  Atom("H", "B", 0.5, 0.5, 0)], "OH2", 0, 1, "O(H)H")
                        ])

        self.assertEqual(mol.get_excluded_pairs(max_exclusion=2), [[[[0, 1], [0, 2]], [[0, 1], [0, 2]]], [[[1, 2]], [[1, 2]]]])

        self.test_passed = True

    def test_to_xyz(self):

        fragment0 = Fragment([Atom("H", "A", 5, 3, 0.00343), Atom("Cl", "B", 2, 0, -13), Atom("He", "C", 6, 2, 0.343)], "HClHe", -1, 1, "H[Cl][He]")

        molecule = Molecule([fragment0])

        self.assertEqual(molecule.to_xyz(), fragment0.to_xyz()[:-1])

        fragment1 = Fragment([Atom("Ar", "D", 0.23430523424, -34, -234.5235)], "Ar", -2, 1, "[Ar]")

        molecule = Molecule([fragment0, fragment1])
        
        self.assertEqual(molecule.to_xyz(), fragment0.to_xyz() + fragment1.to_xyz()[:-1])

        fragment2 = Fragment([Atom("Xe", "E", 0, 0, 0), Atom("Br", "F", 62, 5, 0.001)], "XeBr", 0, 2, "[Xe][Br]")

        molecule = Molecule([fragment0, fragment1, fragment2])

        self.assertEqual(molecule.to_xyz(), fragment0.to_xyz() + fragment1.to_xyz() + fragment2.to_xyz()[:-1])

        self.test_passed = True

suite = unittest.TestLoader().loadTestsFromTestCase(TestMolecule)
