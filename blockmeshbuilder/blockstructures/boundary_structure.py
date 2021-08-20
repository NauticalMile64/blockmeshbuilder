import numpy as np
from .cartesian import CartBlockStruct
from ..blockelements import Point, PolyLineCurvedEdge
import shapely.geometry as shp
# from shapely.ops import nearest_points
from numpy.linalg import norm


p = 1e-8


def create_boundary_block_structure(boundary_points, boundary_indices, ns, zs, nt, nn, nz, closed_boundary=False, **kwargs):

	# boundary_points = np.asarray(boundary_points)
	# if not (np.issubdtype(boundary_points.dtype, np.floating) or np.issubdtype(boundary_points.dtype, np.integer)):
	#	raise TypeError(f'The surface coordinates contained in the boundary_points array array '
	#					f'are not floating point or integer numeric values.')

	# For now assume boundary_points is a 1D sequence of ndarrays of length num_z_divs x num_poly_pts x 3
	# Assume boundary_indices has dimensions num_z_divs x num_block_divs

	boundary_points = list(boundary_points)  # Make a copy of the input points list
	ns = np.asarray(ns)
	zs = np.asarray(zs)
	nt = np.asarray(nt, dtype=np.uint32)
	nn = np.asarray(nn, dtype=np.uint32)
	nz = np.asarray(nz, dtype=np.uint32)

	boundary_indices = np.asarray(boundary_indices, dtype=np.integer)

	if closed_boundary:
		for k in range(zs.size):
			bnd_pts = boundary_points[k]
			bnd_idx = boundary_indices[k]

			if not np.allclose(bnd_pts[-1], bnd_pts[0]):

				bnd_idx[bnd_idx < 0] -= 1
				bnd_idx[bnd_idx >= 0] += 1

				boundary_points[k] = np.pad(bnd_pts, ((1, 1), (0, 0)), mode='wrap')

		boundary_indices = np.pad(boundary_indices, ((0, 0), (0, 1)), mode='constant', constant_values=((0, 0), (0, -1)))
		print(boundary_indices)

	else:
		for k in range(zs.size):
			bnd_pts = boundary_points[k]
			bnd_idx = boundary_indices[k]

			bnd_idx[bnd_idx < 0] -= 1
			bnd_idx[bnd_idx >= 0] += 1

			start_pt = (1 + p) * bnd_pts[0] - p * bnd_pts[1]
			end_pt = (1 + p) * bnd_pts[-1] - p * bnd_pts[-2]
			boundary_points[k] = np.vstack(([start_pt], bnd_pts, [end_pt]))

	num_t_blocks = boundary_indices.shape[1]
	tangent_vecs = np.zeros((num_t_blocks, 3))
	offset_indices = np.zeros(num_t_blocks, dtype=np.integer)

	boundary_struct = CartBlockStruct(np.arange(num_t_blocks), ns, zs, nt, nn, nz, **kwargs)
	bl_edges = boundary_struct.edges
	bl_b_vts = boundary_struct.baked_vertices

	for k in range(zs.size):
		boundary_pts_z = np.asarray(boundary_points[k])
		num_b_pts = boundary_pts_z.shape[0]
		boundary_indices_z = boundary_indices[k] % num_b_pts

		disp_vecs = np.diff(boundary_pts_z, axis=0, append=((1 + p) * boundary_pts_z[-1] - p * boundary_pts_z[-2],))
		disp_vecs /= norm(disp_vecs, axis=1)[:, np.newaxis]

		for i in range(num_t_blocks):
			b_idx = boundary_indices_z[i]
			tangent_vecs[i, :2] = np.sum(disp_vecs[b_idx - 1:b_idx + 1], axis=0)

		unit_norm_vecs = np.cross(tangent_vecs, (0, 0, 1))[:, :2]
		unit_norm_vecs /= norm(unit_norm_vecs, axis=1)[:, np.newaxis]

		if closed_boundary:
			base_geometry = shp.LinearRing(boundary_pts_z)
		else:
			base_geometry = shp.LineString(boundary_pts_z)

		# Create offset curves
		for n in range(ns.size):
			if np.isclose(ns[n], 0.):
				bl_pts = boundary_pts_z
				offset_indices[:] = boundary_indices_z
				boundary_struct.vertices[:, n, k, :2] = bl_pts[offset_indices]
			else:

				if closed_boundary:
					# TODO: add a try statement to catch cases where multiple linear rings created, and pick exterior
					x, y = base_geometry.buffer(ns[n]).exterior.coords.xy
					bl_pts = np.array((x, y)).T[:-1]  # Exclude the last point, which is co-incident with the first

				else:
					geometry_string = base_geometry.parallel_offset(ns[n], 'left')
					try:
						bl_line_string = geometry_string[1]
					except TypeError:
						bl_line_string = geometry_string

					bl_pts = np.array(bl_line_string)

				c_pt = boundary_pts_z[1] - unit_norm_vecs[0] * ns[n]
				zero_index = (norm(bl_pts - c_pt, axis=1)).argmin()
				bl_pts = np.roll(bl_pts, -zero_index, axis=0)
				plt.plot(*c_pt, color=f'C{n}', marker='v')
				bl_pts = np.pad(bl_pts, ((0, 1), (0, 0)), mode='wrap')

				for i in range(boundary_indices_z.size):
					if i == 0 and boundary_indices_z[0] == 1:
						nearest_index = 0
						print(bl_pts[nearest_index])
					elif closed_boundary and i == num_t_blocks - 1:
						nearest_index = bl_pts.shape[0] - 1
						print(bl_pts[nearest_index] - bl_pts[0])
						print(f'last {bl_pts[-1]}')
					else:
						close_pt = boundary_pts_z[boundary_indices_z[i]] - unit_norm_vecs[i] * ns[n]
						nearest_index = (norm(bl_pts - close_pt, axis=1)).argmin()

					boundary_struct.vertices[i, n, k, :2] = bl_pts[nearest_index]
					offset_indices[i] = nearest_index

			for i in range(boundary_indices_z.size - 1):
				if offset_indices[i + 1] - offset_indices[i] > 1:
					pts = [Point(np.append(s_pt, zs[k])) for s_pt in bl_pts[offset_indices[i] + 1:offset_indices[i + 1]]]
					bl_edges[i, n, k, 0] = PolyLineCurvedEdge(bl_b_vts[i:i + 2, n, k], pts)

	if closed_boundary:
		boundary_struct.baked_vertices[-1] = boundary_struct.baked_vertices[0]

	return boundary_struct
