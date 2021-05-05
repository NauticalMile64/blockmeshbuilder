import numpy as np
import warnings
from ..blockelements import cart_conv_pair, Vertex, ProjectionEdge, Face, HexBlock
from ..grading import SimpleGrading, EdgeGrading, uniformGrading, uniformGradingElement
from ..tags import ZoneTag, DEFAULT_ZONE_TAG


headers = ['vertices', 'num_divisions', 'grading', 'baked_vertices', 'edges', 'faces',
		   'block_mask', 'vertex_mask', 'edge_mask', 'face_mask', 'zone_tags', 'boundary_tags']
formats = ['3f8', '3u4', '3O', 'O', '3O', '3O', '?', '?', '3?', '3?', 'O', '3O']
block_struct_dtype = np.dtype({'names': headers, 'formats': formats})
init_pos = np.arange(3)
init_pos.setflags(write=False)


class BaseBlockStruct(np.recarray):

	def __new__(cls, x0, x1, x2, nd0, nd1, nd2, conv_funcs=cart_conv_pair, zone_tag=DEFAULT_ZONE_TAG):
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
				raise IndexError(f'The length of the number of cell divisions array of the {arr_order} array is invalid.'
								 f'\nIt must be equal to or one element shorter than the corresponding dimension array.'
								 f'\n{nd_array_error_msg}')

			if np.any(nd <= 0):
				raise ValueError(f'The {arr_order} number of cell divisions array contains '
								 f'a non-positive integer. {nd_array_error_msg}')

		# Adjust number of divisions arrays if necessary
		nd0 = nd0 if (nd0.ndim == 0 or nd0.size == x0.size) else np.append(nd0, 0)
		nd1 = nd1 if (nd1.ndim == 0 or nd1.size == x1.size) else np.append(nd1, 0)
		nd2 = nd2 if (nd2.ndim == 0 or nd2.size == x2.size) else np.append(nd2, 0)

		if not callable(conv_funcs[0]) or not callable(conv_funcs[1]):
			raise TypeError(f'The conversion function object passed is not callable.')

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
			block_structure.baked_vertices[ind] = Vertex(vts[ind], conv_funcs)

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
