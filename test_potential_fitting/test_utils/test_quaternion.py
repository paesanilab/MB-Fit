import unittest, random, math

from potential_fitting.utils import quaternion

class TestQuaternion(unittest.TestCase):

    def test_add(self):
        for i in range(100):
            r1 = random.random()
            i1 = random.random()
            j1 = random.random()
            k1 = random.random()

            r2 = random.random()
            i2 = random.random()
            j2 = random.random()
            k2 = random.random()

            q1 = quaternion.Quaternion(r1, i1, j1, k1)
            q2 = quaternion.Quaternion(r2, i2, j2, k2)

            q3 = q1 + q2

            self.assertEqual(q3.get_r(), r1 + r2)
            self.assertEqual(q3.get_i(), i1 + i2)
            self.assertEqual(q3.get_j(), j1 + j2)
            self.assertEqual(q3.get_k(), k1 + k2)

    def test_sub(self):
        for i in range(100):
            r1 = random.random()
            i1 = random.random()
            j1 = random.random()
            k1 = random.random()

            r2 = random.random()
            i2 = random.random()
            j2 = random.random()
            k2 = random.random()

            q1 = quaternion.Quaternion(r1, i1, j1, k1)
            q2 = quaternion.Quaternion(r2, i2, j2, k2)

            q3 = q1 - q2

            self.assertEqual(q3.get_r(), r1 - r2)
            self.assertEqual(q3.get_i(), i1 - i2)
            self.assertEqual(q3.get_j(), j1 - j2)
            self.assertEqual(q3.get_k(), k1 - k2)

    def test_mul(self):
        for i in range(100):
            r1 = random.random()
            i1 = random.random()
            j1 = random.random()
            k1 = random.random()

            r2 = random.random()
            i2 = random.random()
            j2 = random.random()
            k2 = random.random()

            q1 = quaternion.Quaternion(r1, i1, j1, k1)
            q2 = quaternion.Quaternion(r2, i2, j2, k2)

            q3 = q1 * q2

            self.assertEqual(q3.get_r(), r1 * r2 - i1 * i2 - j1 * j2 - k1 * k2)
            self.assertEqual(q3.get_i(), r1 * i2 + i1 * r2 + j1 * k2 - k1 * j2)
            self.assertEqual(q3.get_j(), r1 * j2 - i1 * k2 + j1 * r2 + k1 * i2)
            self.assertEqual(q3.get_k(), r1 * k2 + i1 * j2 - j1 * i2 + k1 * r2)

    def test_abs(self):
        for i in range(100):
            r1 = random.random()
            i1 = random.random()
            j1 = random.random()
            k1 = random.random()

            q1 = quaternion.Quaternion(r1, i1, j1, k1)

            self.assertEqual(abs(q1), math.sqrt(r1 ** 2 + i1 ** 2 + j1 ** 2 + k1 ** 2))

    def test_neg(self):
        for i in range(100):
            r1 = random.random()
            i1 = random.random()
            j1 = random.random()
            k1 = random.random()

            q1 = quaternion.Quaternion(r1, i1, j1, k1)

            q2 = -q1

            self.assertEqual(q1.get_r(), -q2.get_r())
            self.assertEqual(q1.get_i(), -q2.get_i())
            self.assertEqual(q1.get_j(), -q2.get_j())
            self.assertEqual(q1.get_k(), -q2.get_k())
            self.assertEqual(abs(q1), abs(q2))

    def test_pos(self):
        for i in range(100):
            r1 = random.random()
            i1 = random.random()
            j1 = random.random()
            k1 = random.random()

            q1 = quaternion.Quaternion(r1, i1, j1, k1)

            q2 = +q1

            self.assertEqual(q1.get_r(), q2.get_r())
            self.assertEqual(q1.get_i(), q2.get_i())
            self.assertEqual(q1.get_j(), q2.get_j())
            self.assertEqual(q1.get_k(), q2.get_k())
            self.assertEqual(abs(q1), abs(q2))

    def test_eq(self):
        self.assertTrue(quaternion.Quaternion(1, 2, 3, 4) == quaternion.Quaternion(1, 2, 3, 4))
        self.assertFalse(quaternion.Quaternion(2, 2, 3, 4) == quaternion.Quaternion(1, 2, 3, 4))
        self.assertFalse(quaternion.Quaternion(1, 3, 3, 4) == quaternion.Quaternion(1, 2, 3, 4))
        self.assertFalse(quaternion.Quaternion(1, 2, 4, 4) == quaternion.Quaternion(1, 2, 3, 4))
        self.assertFalse(quaternion.Quaternion(1, 2, 3, 5) == quaternion.Quaternion(1, 2, 3, 4))
        self.assertFalse(quaternion.Quaternion(3, 4, 1, 2) == quaternion.Quaternion(1, 2, 3, 4))

    def test_ne(self):
        self.assertFalse(quaternion.Quaternion(1, 2, 3, 4) != quaternion.Quaternion(1, 2, 3, 4))
        self.assertTrue(quaternion.Quaternion(2, 2, 3, 4) != quaternion.Quaternion(1, 2, 3, 4))
        self.assertTrue(quaternion.Quaternion(1, 3, 3, 4) != quaternion.Quaternion(1, 2, 3, 4))
        self.assertTrue(quaternion.Quaternion(1, 2, 4, 4) != quaternion.Quaternion(1, 2, 3, 4))
        self.assertTrue(quaternion.Quaternion(1, 2, 3, 5) != quaternion.Quaternion(1, 2, 3, 4))
        self.assertTrue(quaternion.Quaternion(3, 4, 1, 2) != quaternion.Quaternion(1, 2, 3, 4))

    def test_normalize(self):
        for i in range(100):
            r1 = random.random()
            i1 = random.random()
            j1 = random.random()
            k1 = random.random()

            q1 = quaternion.Quaternion(r1, i1, j1, k1)

            q2 = q1.normalize()

            self.assertAlmostEqual(abs(q2), 1)
            self.assertAlmostEqual(q2.r, q1.r / abs(q1))
            self.assertAlmostEqual(q2.i, q1.i / abs(q1))
            self.assertAlmostEqual(q2.j, q1.j / abs(q1))
            self.assertAlmostEqual(q2.k, q1.k / abs(q1))

    def test_conjugate(self):
        for i in range(100):
            r1 = random.random()
            i1 = random.random()
            j1 = random.random()
            k1 = random.random()

            q1 = quaternion.Quaternion(r1, i1, j1, k1)

            q2 = q1.conjugate()

            self.assertEqual(q1.r, q2.r)
            self.assertEqual(q1.i, -q2.i)
            self.assertEqual(q1.j, -q2.j)
            self.assertEqual(q1.k, -q2.k)
            self.assertEqual(abs(q1), abs(q2))

    def test_rotate(self):

        self.assertEqual(quaternion.Quaternion(0, 0, 1, 0).rotate(-1, 1, 1), (1, 1, -1))
        self.assertEqual(quaternion.Quaternion(0, 1, 0, 0).rotate(-1, 1, 1), (-1, -1, -1))
        self.assertEqual(quaternion.Quaternion(0, 0, 0, 1).rotate(-1, 1, 1), (1, -1, 1))
        x1, y1, z1 = -1, 1, 1
        x2, y2, z2 = quaternion.Quaternion( 0.70710678, 0,  0.70710678, 0).rotate(x1, y1, z1)
        self.assertAlmostEqual(x2, 1)
        self.assertAlmostEqual(y2, 1)
        self.assertAlmostEqual(z2, 1)


suite = unittest.TestLoader().loadTestsFromTestCase(TestQuaternion)