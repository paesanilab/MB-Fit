import math

class Quaternion(object):
    
    def __init__(self, r, i, j, k):
        self.r = r
        self.i = i
        self.j = j
        self.k = k

    def __add__(self, other):
        return Quaternion(self.r + other.r, self.i + other.i, self.j + other.j, self.k + other.k)

    def __sub__(self, other):
        return Quaternion(self.r - other.r, self.i - other.i, self.j - other.j, self.k - other.k)

    def __mul__(self, other):
        return Quaternion(
            self.r * other.r - self.i * other.i - self.j * other.j - self.k * other.k,
            self.r * other.i + self.i * other.r + self.j * other.k - self.k * other.j,
            self.r * other.j - self.i * other.k + self.j * other.r + self.k * other.r,
            self.r * other.k + self.i * other.j - self.j * other.i + self.k * other.r
        )

    def __abs__(self):
        return math.sqrt(self.r ** 2 + self.i ** 2 + self.j ** 2 + self.k ** 2)

    def __neg__(self):
        return Quaternion(-self.r, -self.i, -self.j, -self.k)

    def __pos__(self):
        return Quaternion(+self.r, +self.i, +self.j, +self.k)

    def __eq__(self, other):
        return self.r == other.r and self.i == other.i and self.j == other.j and self.k == other.k

    def __ne__(self, other):
        return not self == other

    def conjugate(self):
        return Quaternion(self.r, -self.i, -self.j, -self.k)

    def rotate(self, x, y, z):
        vector_quaternion = Quaternion(0, x, y, z)

        rotated_quaternion = self * vector_quaternion * self.conjugage()

        x = rotated_quaternion.i
        y = rotated_quaternion.j
        z = rotated_quaternion.k

        return x, y, z
