class Model:
	def __init__(self, method, basis, cp = False):
		"""
		Constructor for the Model class.

		Args:
			method			- The method of this model.
			basis			- The basis set of this model.
			cp				- Is counterpoise correction used in this model?
				Default: False

		Returns:
			None.
		"""

		self.method, self.basis, self.cp = method, basis, cp

	def get_method(self):
		"""
		Retrieves the method of this model.

		Args:
			None.

		Returns:
			None.
		"""

		return self.method

	def get_basis(self):
		"""
		Retrieves the basis set of this model.

		Args:
			None.

		Returns:
			None.
		"""

		return self.basis

	def get_cp(self):
		"""
		Retrieves whether this model uses counterpoise correction.

		Args:
			None.

		Returns:
			None.
		"""

		return self.cp

	def __equ__(self, other):
		return self.method == other.method and self.basis == other.basis and self.cp == other.cp