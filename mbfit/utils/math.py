def test_difference_under_threshold(val1, val2, threshold):
	"""
	Tests if the difference between val1 and val2 is under the threshold

	Args:
		val1				- The first value to test.
		val2				- The second value to test.
		threshold 			- The cutoff threshold.

	Returns:
		True if the difference is strictly less than the threhsold, false otherwise
	"""

	return abs(val1 - val2) < threshold