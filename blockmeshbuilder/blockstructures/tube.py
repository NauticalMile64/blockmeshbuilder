import numpy as np
import warnings
from ..blockelements import cyl_conv_pair, Point, ProjectionEdge, Face
from ..geometry import Cylinder, Cone
from .cartesian import BaseBlockStruct


def are_coterminal(a, b, tol):
	modulus = 2 * np.pi
	diff = (a - b) % modulus
	return min(diff, modulus - diff) <= tol


class TubeBlockStruct(BaseBlockStruct):

	def __new__(cls, rs, ts, zs, nr, nt, nz, is_complete=False, angle_r_tol=1e-8, **kwargs):

		if np.any(np.asarray(rs) < 0):
			raise ValueError('Negative values supplied to array of radial values in TubeBlockStruct.')

		are_first_and_last_close = are_coterminal(ts[0], ts[-1], angle_r_tol)
		if not are_first_and_last_close and is_complete:
			warnings.warn(f'TubeBlockStruct had the is_complete flag raised, while the '
						  f'first and last angles are unequal; make sure these are separated by 2*pi. '
						  f'Alternatively, increase angle_tol.'
						  f'Setting is_complete=False')
			is_complete = False

		if not is_complete and are_first_and_last_close:
			warnings.warn(f'TubeBlockStruct had the is_complete flag down, while the '
						  f'first and last angles are equal or nearly equal; separated by a value less than 2*pi. '
						  f'The resulting tube struct may visually appear closed, but a circumferential \'wall\' '
						  f'will be present at theta={ts[0]}.')

		if is_complete and len(ts) < 3:
			raise ValueError(f'Too few elements in ts array supplied to TubeBlockStruct. At least 3 elements are '
							 f'required when the is_complete flag is raised.')

		if np.any(np.diff(ts) >= np.pi):
			warnings.warn(f'The ts array passed to TubeBlockStruct contains neighbouring angles '
						  f'separated by approximately 180 degrees or more. When these edges are rendered in blockmesh, '
						  f'either an error will be triggered, or it may constitute a degenerate case.')

		block_structure = super(TubeBlockStruct, cls).__new__(cls, rs, ts, zs, nr, nt, nz, cyl_conv_pair, **kwargs)

		b_vts = block_structure.baked_vertices
		edges = block_structure.edges
		faces = block_structure.faces
		rshape = np.array(block_structure.shape) - 1

		if is_complete:
			b_vts[:, -1] = b_vts[:, 0]

		block_structure.is_full = np.isclose(rs[0], 0.)

		# Re-assign vertices, edges, and faces at the axis of the tube
		if block_structure.is_full:
			b_vts[0] = b_vts[0, 0]
			block_structure.face_mask[0, ..., 0] = True
			block_structure.edge_mask[0, ..., 1] = True

			for j in range(rshape[1] + 1):
				for k in range(rshape[2] + 1):
					edges[0, j, k, 0] = ProjectionEdge(b_vts[0:2, j, k])

			for j in range(rshape[1] + 1):
				for k in range(rshape[2]):
					faces[0, j, k, 1] = Face(b_vts[0:2, j, k:k + 2])

		block_structure.is_complete = is_complete

		return block_structure

	def __array_finalize__(self, other_block_structure):
		if (other_block_structure is None) or (self.dtype is not other_block_structure.dtype):
			return

		# Right now this function gets called 3 times after slicing by recarray internals,
		# but it shouldn't impact the user

		# Assume index ordering hasn't been changed using a function such as np.moveaxis
		self.is_complete = other_block_structure.is_complete and \
						   np.all(self.baked_vertices[:, 0] == self.baked_vertices[:, -1])

		self.is_full = other_block_structure.is_full and \
					   np.all(np.isin(self.baked_vertices[0], other_block_structure.baked_vertices[0]))

		BaseBlockStruct.__array_finalize__(self, other_block_structure)

	def write(self, block_mesh_dict):

		vts = self.vertices

		if self.is_full and not np.all(np.isclose(vts[0, ..., 0], 0.)):
			warnings.warn(f'When initialized, the inner radii in TubeBlockStruct in zone_tag {self.zone_tag} were '
						  f'set to 0, but some were changed before writing. The nodes along the centerline of '
						  f'the tube may not be positioned as expected.')

		# Define cylinder axis in canonical (untransformed) coordinate system
		axis_length = 1e5  # Large value to ensure it spans the domain

		# Start and end points along the cylindrical axis (Z-axis) in curvilinear coords
		canonical_start = np.array([0., 0., -axis_length])
		canonical_end = np.array([0., 0., axis_length])

		# Convert to Cartesian, apply transform, create Points
		start_cart = self.conv_funcs[0](canonical_start)  # cyl -> cart
		end_cart = self.conv_funcs[0](canonical_end)

		# Apply the transformation (rotation + translation + scale)
		transformed_start = self.transform.apply(start_cart.reshape(1, 3)).ravel()
		transformed_end = self.transform.apply(end_cart.reshape(1, 3)).ravel()

		s_pt = Point(transformed_start)
		e_pt = Point(transformed_end)

		# Create cylinders for each unique radius
		cyls = {}
		for i, r in np.ndenumerate(np.unique(vts[..., 0])):
			if not np.isclose(r, 0.):
				# Radius also needs to be scaled if non-uniform scaling is used
				# For uniform scaling, all radii scale the same
				# For non-uniform, we use the average of X and Y scale factors
				# (assuming cylindrical coords are in X-Y plane)
				scaled_r = r * np.mean(self.transform.scale[:2])
				cyls[r] = Cylinder(s_pt, e_pt, scaled_r, f'blockcyl-{i}')

		# Mask axial and radial edges as well as circumferential faces so
		# no redundant edges or faces are written to file
		if self.is_complete:
			self.edge_mask[:, -1, :, [0, 2]] = True
			self.face_mask[:, -1, :, 1] = True

		b_vts = self.baked_vertices
		vertex_mask = self.vertex_mask
		proj_rcrds = vts[..., 0]

		shape = self.shape
		shp = (shape[0], shape[1] - 1 if self.is_complete else shape[1], shape[2])
		for ind in np.ndindex(shp):
			if not (vertex_mask[ind] or np.isclose(proj_rcrds[ind], 0.)):
				b_vts[ind].proj_geom(cyls[proj_rcrds[ind]])

		# Test and see if the circumferential edges or axial edges need to be projected
		a_edges = self.edges[..., 2]
		a_edge_mask = self.edge_mask[..., 2]
		c_edges = self.edges[..., 1]
		c_edge_mask = self.edge_mask[..., 1]
		acrds = self.vertices[..., 1]

		cones = {}
		are_cones = Cone in block_mesh_dict.of_distribution_features['geometries']

		for ind in np.ndindex(shape):

			a_edge = a_edges[ind]
			if (not a_edge_mask[ind]) and isinstance(a_edge, ProjectionEdge):
				nind = (ind[0], ind[1], ind[2] + 1)
				if a_edge.is_relevant() or (not np.allclose(acrds[ind], acrds[nind])):
					if not np.allclose(proj_rcrds[ind], proj_rcrds[nind]):
						if are_cones:
							cone_tuple = (proj_rcrds[ind], proj_rcrds[nind], acrds[ind], acrds[nind])
							if cone_tuple not in cones:

								# Cone endpoints at specific Z-coordinates in cylindrical coords
								cone_start_cyl = np.array([0., 0., acrds[ind]])
								cone_end_cyl = np.array([0., 0., acrds[nind]])

								# Convert to Cartesian
								cone_start_cart = self.conv_funcs[0](cone_start_cyl)
								cone_end_cart = self.conv_funcs[0](cone_end_cyl)

								# Apply transformation
								cone_start_transformed = self.transform.apply(cone_start_cart.reshape(1, 3)).ravel()
								cone_end_transformed = self.transform.apply(cone_end_cart.reshape(1, 3)).ravel()

								# Create Points
								cone_s_pt = Point(cone_start_transformed)
								cone_e_pt = Point(cone_end_transformed)

								# Scale radii
								scaled_r1 = proj_rcrds[ind] * np.mean(self.transform.scale[:2])
								scaled_r2 = proj_rcrds[nind] * np.mean(self.transform.scale[:2])

								cones[cone_tuple] = Cone(
									cone_s_pt, cone_e_pt,
									scaled_r1, scaled_r2,
									'block_cone'
								)
							a_edge.proj_geom(cones[cone_tuple])

					elif not np.allclose(proj_rcrds[ind], 0.):
						a_edge.proj_geom(cyls[proj_rcrds[ind]])

			c_edge = c_edges[ind]
			if (not c_edge_mask[ind]) and isinstance(c_edge, ProjectionEdge):
				c_edge.proj_geom(cyls[proj_rcrds[ind]])

		BaseBlockStruct.write(self, block_mesh_dict)
