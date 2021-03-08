# -*- coding: future_fstrings -*-
from .core import *
import numpy as np
import warnings

from sys import version_info

if version_info[0] == 2:
	from itertools import izip as zip

headers = ['vertices', 'num_divisions', 'grading', 'baked_vertices', 'edges', 'faces',
		   'block_mask', 'vertex_mask', 'edge_mask', 'face_mask', 'zone_tags', 'boundary_tags']
formats = ['3f8', '3u4', '3O', 'O', '3O', '3O', '?', '?', '3?', '3?', 'O', '3O']
block_struct_dtype = np.dtype({'names': headers, 'formats': formats})
init_pos = np.arange(3)
init_pos.setflags(write=False)


def wrap_radians(values):
	return values % (2 * np.pi)


class BaseBlockStruct(np.recarray):

	def __new__(cls, x0, x1, x2, nd0, nd1, nd2, conv_func=cart_to_cart, zone_tag=DEFAULT_ZONE_TAG):
		x0 = np.asarray(x0)
		x1 = np.asarray(x1)
		x2 = np.asarray(x2)

		nd0 = np.asarray(nd0, dtype=np.int32)
		nd1 = np.asarray(nd1, dtype=np.int32)
		nd2 = np.asarray(nd2, dtype=np.int32)

		# Error Checking
		dim_array_error_msg = '\nDimension arrays must be 1D arrays of floating points/integers, ' \
							  '\nlisted in ascending order, with a minimum of 2 elements.'
		nd_array_error_msg = '\nThe number of cell divisions parameters must be positive integers or 1D arrays thereof, ' \
							 '\nhaving lengths equal to or one less than the corresponding dimension array.'
		for x, nd, arr_order in zip((x0, x1, x2), (nd0, nd1, nd2), ('first', 'second', 'third')):
			if not (np.issubdtype(x.dtype, np.floating) or np.issubdtype(x.dtype, np.integer)):
				raise TypeError(f'The values contained in the {arr_order} dimension array '
								f'are not floating point or integer numeric values. {dim_array_error_msg}')

			if x.ndim != 1:
				raise IndexError(f'The {arr_order} dimension array is not 1D. {dim_array_error_msg} '
								 f'\nSee examples for moving vertices after block structure creation.')

			if x.size < 2:
				raise IndexError(f'The {arr_order} dimension array has fewer than 2 elements. '
								 f'\nAt least two coordinates are needed for each dimension array to define the end '
								 f'points of a block. {dim_array_error_msg}')

			if np.any(np.diff(x) < 0.0):
				raise ValueError(f'The {arr_order} dimension array is in non-ascending order. {dim_array_error_msg}')

			if nd.ndim > 1:
				raise IndexError(f'The {arr_order} number of cell divisions array has two or more dimensions.'
								 f'{nd_array_error_msg}')

			if nd.ndim == 1 and not (x.size - 1 <= nd.size <= x.size):
				raise IndexError(f'The length of the number of cell divisions array is invalid. It must be equal to or '
								 f'\none element shorter than the corresponding dimension array. {nd_array_error_msg}')

			if np.any(nd <= 0):
				raise ValueError(f'The {arr_order} number of cell divisions array contains '
								 f'a non-positive integer. {nd_array_error_msg}')

		# Adjust number of divisions arrays if necessary
		nd0 = nd0 if (nd0.ndim == 0 or nd0.size == x0.size) else np.append(nd0, 0)
		nd1 = nd1 if (nd1.ndim == 0 or nd1.size == x1.size) else np.append(nd1, 0)
		nd2 = nd2 if (nd2.ndim == 0 or nd2.size == x2.size) else np.append(nd2, 0)

		if isinstance(zone_tag, str):
			zone_tag = ZoneTag(zone_tag)
		elif not isinstance(zone_tag, ZoneTag):
			raise TypeError(f'The zone_tag parameter was neither a string nor a ZoneTag.')

		shape = (x0.size, x1.size, x2.size)
		block_structure = super(BaseBlockStruct, cls).__new__(cls, shape, dtype=block_struct_dtype)
		rshape = (shape[0] - 1, shape[1] - 1, shape[2] - 1)

		# Initialize vertices
		X0, X1, X2 = np.meshgrid(x0, x1, x2, indexing='ij')

		vts = block_structure.vertices
		vts[..., 0] = X0
		vts[..., 1] = X1
		vts[..., 2] = X2

		for ind in np.ndindex(shape):
			block_structure.baked_vertices[ind] = Vertex(vts[ind], conv_func)

		# Initialize number of divisions
		ND0, ND1, ND2 = np.meshgrid(nd0, nd1, nd2, indexing='ij')

		nds = block_structure.num_divisions
		nds[..., 0] = ND0
		nds[..., 1] = ND1
		nds[..., 2] = ND2

		# Initialize grading
		block_structure.grading[:] = uniformGradingElement

		# Initialize edges and faces
		for s in range(3):
			roll_pos = np.roll(init_pos, s)

			d_edges = np.moveaxis(block_structure.edges[..., s], init_pos, roll_pos)
			d_faces = np.moveaxis(block_structure.faces[..., s], init_pos, roll_pos)
			d_vts = np.moveaxis(block_structure.baked_vertices, init_pos, roll_pos)

			for i in range(rshape[s]):
				for j in range(shape[(s + 1) % 3]):
					for k in range(shape[(s + 2) % 3]):
						d_edges[i, j, k] = ProjectionEdge(d_vts[i:i + 2, j, k])

			for i in range(shape[s]):
				for j in range(rshape[(s + 1) % 3]):
					for k in range(rshape[(s + 2) % 3]):
						d_faces[i, j, k] = Face(d_vts[i, j:j + 2, k:k + 2])

		block_structure.zone_tags[:] = zone_tag

		return block_structure

	def project_structure(self, direction, geometry, face_ind):

		# Get the subarray relevant to the face being projected
		roll_pos = np.roll(init_pos, direction)
		struct = np.moveaxis(self, init_pos, roll_pos)[face_ind]
		shape = struct.shape
		rshape = np.roll(np.array(self.shape) - 1, -direction)[1:]

		# Project vertices
		b_vts = struct.baked_vertices
		vt_mask = struct.vertex_mask
		for ind in np.ndindex(shape):
			if not vt_mask[ind]:
				b_vts[ind].proj_geom(geometry)

		# Project edges
		edges = np.roll(struct.edges, -direction, axis=-1)[..., 1:]
		edge_mask = np.roll(struct.edge_mask, -direction, axis=-1)[..., 1:]
		for j in range(shape[0]):
			for k in range(shape[1]):
				for s in range(2):
					edge = edges[j, k, s]
					if (not edge_mask[j, k, s]) and isinstance(edge, ProjectionEdge):
						edge.proj_geom(geometry)

		# Project faces
		faces = struct.faces[..., direction]
		face_mask = struct.face_mask[..., direction]
		for j in range(rshape[0]):
			for k in range(rshape[1]):
				if not face_mask[j, k]:
					faces[j, k].proj_geom(geometry)

	@staticmethod
	def _get_grading(gt):

		# Get relevant edges
		grd_arr = np.array([
			np.moveaxis(gt[..., s], s, 0)[0].T for s in range(3)
		])

		# Get simplest grading type
		if np.all(grd_arr == uniformGradingElement):
			return uniformGrading

		elif np.all(np.array([grd_arr[s] == grd_arr[s, 0, 0] for s in range(3)])):
			return SimpleGrading(grd_arr[:, 0, 0])

		else:
			for gc in grd_arr:
				gc[1, 0], gc[1, 1] = gc[1, 1], gc[1, 0]

			return EdgeGrading(grd_arr.flatten())

	@staticmethod
	def _get_block_vertices(c_vs):
		return tuple((c_vs[0, 0, 0], c_vs[1, 0, 0],
					  c_vs[1, 1, 0], c_vs[0, 1, 0],

					  c_vs[0, 0, 1], c_vs[1, 0, 1],
					  c_vs[1, 1, 1], c_vs[0, 1, 1]))

	def write(self, block_mesh_dict):

		# Reset block_mask edge
		bmask = self.block_mask
		bmask[-1] = True
		bmask[:, -1] = True
		bmask[:, :, -1] = True

		shape = self.shape
		rshape = (shape[0] - 1, shape[1] - 1, shape[2] - 1)

		# TO DO: See if numpy moving_window function will help here
		for i in range(rshape[0]):
			for j in range(rshape[1]):
				for k in range(rshape[2]):

					if not bmask[i, j, k]:
						# Get sub-array
						block_data = self[i:i + 2, j:j + 2, k:k + 2]

						gt = block_data.grading.copy()
						grading = self._get_grading(gt)

						nd = block_data.num_divisions[0, 0, 0]

						vts = self._get_block_vertices(block_data.baked_vertices)

						block_zone_tag = block_data.zone_tags[0, 0, 0]

						block = HexBlock(vts, nd, block_zone_tag, grading)
						block_mesh_dict.add_hexblock(block)

		b_vts = self.baked_vertices
		for index in np.ndindex(shape):
			block_mesh_dict.add_geometries(b_vts[index].proj_g)

		# write relevant edges and faces
		for s in range(3):
			roll_pos = np.roll(init_pos, s)

			d_bmask = np.moveaxis(bmask, init_pos, roll_pos)
			d_edges = np.moveaxis(self.edges[..., s], init_pos, roll_pos)
			d_faces = np.moveaxis(self.faces[..., s], init_pos, roll_pos)
			d_boundary_tags = np.moveaxis(self.boundary_tags[..., s], init_pos, roll_pos)

			for i in range(rshape[s]):
				for j in range(shape[(s + 1) % 3]):
					for k in range(shape[(s + 2) % 3]):
						edge = d_edges[i, j, k]

						if not edge.is_relevant():
							continue

						lmsk = d_bmask[i, max(j - 1, 0):j + 1, max(k - 1, 0):k + 1]
						if not np.all(lmsk):
							block_mesh_dict.add_edge(edge)
							if isinstance(edge, ProjectionEdge):
								block_mesh_dict.add_geometries(edge.proj_g)

			for i in range(shape[s]):
				for j in range(rshape[(s + 1) % 3]):
					for k in range(rshape[(s + 2) % 3]):
						face = d_faces[i, j, k]

						lmsk = d_bmask[max(i - 1, 0):i + 1, j, k]
						if not np.all(lmsk):

							if d_boundary_tags[i, j, k]:
								if i > 0 and lmsk.sum() == 0:
									warnings.warn(
										f'Attempting to specify a boundary at an interior face. '
										f'Please check assignment to boundary_tags.')
								else:
									block_mesh_dict.add_boundary_face(d_boundary_tags[i, j, k], face)

							if face.proj_g:
								block_mesh_dict.add_face(face)
								block_mesh_dict.add_geometries((face.proj_g,))


# There is truly no need to subclass the BaseBlockStruct into CartBlockStruct
CartBlockStruct = BaseBlockStruct


class TubeBlockStruct(BaseBlockStruct):

	def __new__(cls, rs, ts, zs, nr, nt, nz, zone_tag=DEFAULT_ZONE_TAG, is_complete=False):

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

		block_structure = super(TubeBlockStruct, cls).__new__(cls, rs, ts, zs, nr, nt, nz, cyl_to_cart, zone_tag)

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

	def write(self, block_mesh_dict):

		shape = self.shape
		shp = tuple((shape[0], shape[1] - 1, shape[2]))

		vts = self.vertices

		if self.is_full and not np.all(np.isclose(vts[0, ..., 0], 0.)):
			warnings.warn(f'When initialized, the inner radii in TubeBlockStruct in zone_tag {self.zone_tag} were '
						  f'set to 0, but some were changed before writing. The nodes along the centerline of '
						  f'the tube may not be positioned as expected.')

		cyls = {}
		s_pt = Point([0, 0, -1e5])
		e_pt = Point([0, 0, 1e5])
		for i, r in np.ndenumerate(np.unique(vts[..., 0])):
			if np.isclose(r, 0.):
				continue

			cyls[r] = Cylinder(s_pt, e_pt, r, f'blockcyl-{i}')

		# Mask axial and radial edges as well as circumferential faces so
		# no redundant edges or faces are written to file
		if self.is_complete:
			self.edge_mask[:, -1, :, [0, 2]] = True
			self.face_mask[:, -1, :, 1] = True

		b_vts = self.baked_vertices
		vertex_mask = self.vertex_mask
		proj_rcrds = self.vertices[..., 0]

		for ind in np.ndindex(shp):
			if not (vertex_mask[ind] or np.isclose(proj_rcrds[ind], 0.)):
				b_vts[ind].proj_geom(cyls[proj_rcrds[ind]])

		edges = self.edges
		edge_mask = self.edge_mask
		acrds = self.vertices[..., 1]

		# Test and see if the circumferential edges or axial edges need to be projected
		a_edges = edges[..., 2]
		a_edge_mask = edge_mask[..., 2]
		c_edges = edges[..., 1]
		c_edge_mask = edge_mask[..., 1]

		cones = {}
		are_cones = Cone in block_mesh_dict.of_available_geometries

		for ind in np.ndindex(shape):

			a_edge = a_edges[ind]
			if (not a_edge_mask[ind]) and isinstance(a_edge, ProjectionEdge):
				nind = (ind[0], ind[1], ind[2] + 1)
				if a_edge.is_relevant() or (not np.allclose(acrds[ind], acrds[nind])):
					if not np.allclose(proj_rcrds[ind], proj_rcrds[nind]):
						if are_cones:
							cone_tuple = (proj_rcrds[ind], proj_rcrds[nind], acrds[ind], acrds[nind])
							if cone_tuple not in cones:
								cones[cone_tuple] = Cone(Point((0, 0, acrds[ind])), Point((0, 0, acrds[nind])),
														 proj_rcrds[ind], proj_rcrds[nind], 'block_cone')
							a_edge.proj_geom(cones[cone_tuple])
					else:
						a_edge.proj_geom(cyls[proj_rcrds[ind]])

			c_edge = c_edges[ind]
			if (not c_edge_mask[ind]) and isinstance(c_edge, ProjectionEdge):
				c_edge.proj_geom(cyls[proj_rcrds[ind]])

		BaseBlockStruct.write(self, block_mesh_dict)


dummy_vertex = Vertex((0, 0, 0))
dummy_edge = Edge([dummy_vertex] * 2, name='dummy')
_drt2 = 1. / np.sqrt(2)


class CylBlockStructContainer(object):
	# O-grid curvature offsets for core-oriented cylinders, and tube-oriented cylinders
	_og_core_vectors = np.array([Point([0, -1, 0]), Point([1, 0, 0]), Point([0, 1, 0]), Point([-1, 0, 0])])
	_og_core_vectors.setflags(write=False)
	_og_tube_vectors = np.array([Point([_drt2, _drt2, 0]), Point([-_drt2, _drt2, 0]),
								 Point([-_drt2, -_drt2, 0]), Point([_drt2, -_drt2, 0])])
	_og_tube_vectors.setflags(write=False)

	def __init__(self, rs, ts, zs, nr, nt, nz, zone_tag=DEFAULT_ZONE_TAG, inner_arc_curve=0.25, is_core_aligned=True):

		self.tube_struct = TubeBlockStruct(rs, ts, zs, nr, nt, nz, zone_tag=zone_tag, is_complete=True)

		if np.isclose(rs[0], 0.):
			raise ValueError(f'CylBlockStructContainer in zone_tag {zone_tag} has an inner tube radius of {rs[0]}, '
							 f'which is too close to 0, such that an o-grid cannot be accomodated. '
							 f'Consider using TubeBlockStruct instead.')

		if (len(ts) - 1) % 4 > 0:
			raise ValueError(f'The ts array length provided to CylBlockStructContainer is invalid. In order to make '
							 f'sure the o-grid is properly formed, the number of divisions must be `4*n+1`, where `n` '
							 f'is any positive integer.')

		self.inner_arc_curve = inner_arc_curve
		self.is_core_aligned = is_core_aligned

		num_side_blocks = (len(ts) - 1) // 4
		xs = ys = np.linspace(-rs[0], rs[0], num_side_blocks + 1) * _drt2

		if isinstance(nt, int):
			nx = ny = nt
		else:
			nx = nt[:num_side_blocks]
			ny = nt[num_side_blocks:2 * num_side_blocks]

		self.core_struct = CartBlockStruct(xs, ys, zs, nx, ny, nz, zone_tag=zone_tag)

		if is_core_aligned:  # Rotate the tube to match the core
			self.tube_struct.vertices[..., 1] -= 3 / 4 * np.pi
		else:  # Rotate the core to match the tube
			cyl_vts = np_cart_to_cyl(self.core_struct.vertices)
			cyl_vts[..., 1] -= 5 / 4 * np.pi
			self.core_struct.vertices[:] = np_cyl_to_cart(cyl_vts)

		core_b_vts = self.core_struct.baked_vertices
		tube_b_vts = self.tube_struct.baked_vertices

		# Connect the outer tube structure to the core
		tube_indices = np.arange(ts.size - 1).reshape(4, num_side_blocks)

		for s in range(4):
			tube_b_vts[0, tube_indices[s], :] = np.rot90(core_b_vts, k=-s)[:-1, 0, :]

		tube_b_vts[0, -1, :] = tube_b_vts[0, 0, :]

	def write(self, block_mesh_dict):

		self.tube_struct.face_mask[0, :, :, 0] = True
		self.tube_struct.vertex_mask[0, :, :] = True

		iac = self.inner_arc_curve
		og_vectors = self._og_core_vectors if self.is_core_aligned else self._og_tube_vectors

		if not np.isclose(iac, 0.0):
			tube = self.tube_struct
			shape = tube.shape

			tube_vts = tube.vertices[0]
			tube.edges[0, ..., 1:] = dummy_edge

			core = self.core_struct

			# Create a dictionary of cylinder geometries of the innermost vertices on the tube struct
			cyl_dict = {}
			s_pt = Point([0, 0, -1e5])
			e_pt = Point([0, 0, 1e5])
			cyl_arr = np.empty((shape[2], 4), dtype=Cylinder)
			for k, r in np.ndenumerate(tube_vts[0, :, 0]):
				if r not in cyl_dict:
					u = np.arctan(iac)
					a = r * _drt2
					r_cyl = a / np.sin(u)
					offset = a * (1 / iac - 1)
					offset_axes = og_vectors * offset
					local_cyls = np.array([Cylinder(s_pt - offset_axes[i], e_pt - offset_axes[i],
													r_cyl, f'o-grid-cyl-{r_cyl}-{i}') for i in range(4)])
					cyl_dict[r] = local_cyls

				cyl_arr[k] = cyl_dict[r]

			# For each edge of the O-grid square
			for s in range(4):
				core_side_b_vts = np.rot90(core.baked_vertices, k=-s)[:-1, 0, :]
				core_side_edges = np.rot90(core.edges, k=-s)[::np.sign(3 - (2 * s)), 0, :, s % 2][:-1]

				for index, core_vertex in np.ndenumerate(core_side_b_vts):
					core_vertex.proj_geom(cyl_arr[index[1], s])

				for index, edge in np.ndenumerate(core_side_edges):
					edge.proj_geom(cyl_arr[index[1], s])

		else:
			self.tube_struct.edge_mask[0, :, :, [1, 2]] = True

		self.tube_struct.write(block_mesh_dict)
		self.core_struct.write(block_mesh_dict)
