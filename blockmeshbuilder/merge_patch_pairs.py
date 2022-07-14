from .boundary_tags import BoundaryTag


class _MergePatchPair:
	def __init__(self, large_boundary, small_boundary):

		boundary_pair = (large_boundary, small_boundary)
		for boundary in boundary_pair:
			if not isinstance(boundary, BoundaryTag):
				raise TypeError('both large_boundary and small_boundary must be a BoundaryTag.')

		self.boundary_pair = boundary_pair

	def format(self):
		return f'({self.boundary_pair[0].name} {self.boundary_pair[1].name})'
