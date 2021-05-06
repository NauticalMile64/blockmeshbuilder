from ..blockelements import cart_to_cyl, cyl_to_cart, Point, Vertex, Edge
from ..geometry import Cylinder
from ..zone_tags import DEFAULT_ZONE_TAG
from .cartesian import CartBlockStruct
from .tube import TubeBlockStruct
import numpy as np

dummy_vertex = Vertex((0, 0, 0))
dummy_edge = Edge([dummy_vertex] * 2, name='dummy')
_drt2 = 1. / np.sqrt(2)


class CylBlockStructContainer:
	# O-grid curvature offsets for core-oriented cylinders, and tube-oriented cylinders
	_og_core_vectors = np.array([Point([0, -1, 0]), Point([1, 0, 0]), Point([0, 1, 0]), Point([-1, 0, 0])])
	_og_core_vectors.setflags(write=False)
	_og_tube_vectors = np.array([Point([_drt2, _drt2, 0]), Point([-_drt2, _drt2, 0]),
								 Point([-_drt2, -_drt2, 0]), Point([_drt2, -_drt2, 0])])
	_og_tube_vectors.setflags(write=False)

	def __init__(self, rs, ts, zs, nr, nt, nz, zone_tag=DEFAULT_ZONE_TAG, offset=Point((0., 0., 0.)),
				 inner_arc_curve=0.25, is_core_aligned=True):

		self.tube_struct = TubeBlockStruct(rs, ts, zs, nr, nt, nz, zone_tag=zone_tag, is_complete=True, offset=offset)

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
		self.offset = offset

		num_side_blocks = (len(ts) - 1) // 4
		xs = ys = np.linspace(-rs[0], rs[0], num_side_blocks + 1) * _drt2

		nt_arr = np.asarray(nt, dtype=np.int32)
		if nt_arr.ndim == 0:
			nx = ny = nt_arr
		else:
			nx = nt_arr[:num_side_blocks]
			ny = nt_arr[num_side_blocks:2 * num_side_blocks]

		self.core_struct = CartBlockStruct(xs, ys, zs, nx, ny, nz, zone_tag=zone_tag)

		if is_core_aligned:  # Rotate the tube to match the core
			self.tube_struct.vertices[..., 1] -= 3 / 4 * np.pi
		else:  # Rotate the core to match the tube
			cyl_vts = cart_to_cyl(self.core_struct.vertices)
			cyl_vts[..., 1] -= 5 / 4 * np.pi
			self.core_struct.vertices[:] = cyl_to_cart(cyl_vts)

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
			s_pt = Point([0, 0, -1e5]) + self.offset
			e_pt = Point([0, 0, 1e5]) + self.offset
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

		self.core_struct.vertices += self.offset.get_cart_crds()

		self.tube_struct.write(block_mesh_dict)
		self.core_struct.write(block_mesh_dict)
