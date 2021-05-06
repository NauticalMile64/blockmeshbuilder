import numpy as np
import warnings
from ..blockelements import cyl_conv_pair, cart_to_cyl, cyl_to_cart, Point, ProjectionEdge, Face
from ..geometry import Cylinder, Cone
from ..zone_tags import DEFAULT_ZONE_TAG
from .cartesian import BaseBlockStruct


def wrap_radians(values):
	return values % (2 * np.pi)


class TubeBlockStruct(BaseBlockStruct):

	def __new__(cls, rs, ts, zs, nr, nt, nz, zone_tag=DEFAULT_ZONE_TAG, is_complete=False, offset=Point((0., 0., 0.))):

		if np.any(np.asarray(rs) < 0):
			raise ValueError('Negative values supplied to array of radial values in TubeBlockStruct.')

		are_first_and_last_close = np.isclose(wrap_radians(ts[0]), wrap_radians(ts[-1]))
		if is_complete and ~are_first_and_last_close:
			warnings.warn(f'TubeBlockStruct in zone_tag {zone_tag} had the is_complete flag raised, while the '
						  f'first and last angles are unequal; make sure these are separated by 2*pi. '
						  f'Setting is_complete=False')
			is_complete = False

		if not is_complete and are_first_and_last_close:
			warnings.warn(f'TubeBlockStruct in zone_tag {zone_tag} had the is_complete flag down, while the '
						  f'first and last angles are equal or nearly equal; separated by a value less than 2*pi. '
						  f'The resulting tube struct may visually appear closed, but a circumferential \'wall\' '
						  f'will be present at theta={ts[0]}.')

		if is_complete and len(ts) < 3:
			raise ValueError(f'Too few elements in ts array supplied to TubeBlockStruct. At least 3 elements are '
							 f'required when the is_complete flag is raised.')

		if np.any(np.diff(ts) >= np.pi):
			warnings.warn(f'The ts array passed to TubeBlockStruct in zone_tag {zone_tag} contains neighbouring angles '
						  f'separated by approximately 180 degrees or more. When these edges are rendered in blockmesh, '
						  f'either an error will be triggered, or it may constitute a degenerate case.')

		block_structure = super(TubeBlockStruct, cls).__new__(cls, rs, ts, zs, nr, nt, nz, cyl_conv_pair, zone_tag)

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
		block_structure.offset = offset

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

		self.offset = other_block_structure.offset

	def write(self, block_mesh_dict):

		vts = self.vertices

		if self.is_full and not np.all(np.isclose(vts[0, ..., 0], 0.)):
			warnings.warn(f'When initialized, the inner radii in TubeBlockStruct in zone_tag {self.zone_tag} were '
						  f'set to 0, but some were changed before writing. The nodes along the centerline of '
						  f'the tube may not be positioned as expected.')

		'''
		if np.any(np.diff(self.vertices[..., 1], axis=1) >= np.pi):
			raise ValueError(f'One or more angles assigned to a tube block struct are separated from their neighbour by 180 '
					   f'degrees or more. When these edges are rendered in blockmesh, either an error will be '
					   f'triggered, or it may constitute a degenerate case.')
		'''

		cyls = {}
		s_pt = Point([0, 0, -1e5]) + self.offset
		e_pt = Point([0, 0, 1e5]) + self.offset
		for i, r in np.ndenumerate(np.unique(vts[..., 0])):
			if not np.isclose(r, 0.):
				cyls[r] = Cylinder(s_pt, e_pt, r, f'blockcyl-{i}')

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
								cones[cone_tuple] = Cone(Point((0, 0, acrds[ind])) + self.offset,
														 Point((0, 0, acrds[nind])) + self.offset,
														 proj_rcrds[ind], proj_rcrds[nind], 'block_cone')
							a_edge.proj_geom(cones[cone_tuple])
					else:
						a_edge.proj_geom(cyls[proj_rcrds[ind]])

			c_edge = c_edges[ind]
			if (not c_edge_mask[ind]) and isinstance(c_edge, ProjectionEdge):
				c_edge.proj_geom(cyls[proj_rcrds[ind]])

		vts[:] = cart_to_cyl(cyl_to_cart(vts) + self.offset.get_cart_crds())

		BaseBlockStruct.write(self, block_mesh_dict)
