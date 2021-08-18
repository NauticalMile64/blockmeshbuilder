import numpy as np
from .cartesian import CartBlockStruct
from ..blockelements import Point, PolyLineCurvedEdge
import shapely.geometry as shp
# from shapely.ops import nearest_points
from numpy.linalg import norm

p = 1e-7


def create_boundary_block_structure(boundary_points, boundary_indices, ns, zs, nt, nn, nz, closed_boundary=False, **kwargs):

	# For now assume boundary_points is a 1D sequence of ndarrays of length num_z_divs x num_poly_pts x 3
	# Assume boundary_indices has dimensions num_z_divs x num_block_divs

	boundary_indices = np.asarray(boundary_indices, dtype=np.integer)
	num_t_blocks = boundary_indices.shape[1]

	if closed_boundary:
		for k in range(zs.size):
			bnd_pts = boundary_points[k]
			bnd_idx = boundary_indices[k]
			if not np.allclose(bnd_pts[-1], bnd_pts[0]):
				bnd_idx[bnd_idx < 0] -= 4
				tiny_vec1 = p * (bnd_pts[-1] - bnd_pts[-2])
				tiny_vec2 = p * (bnd_pts[1] - bnd_pts[0])
				app_arr = np.array([bnd_pts[-1] + tiny_vec1 - bnd_pts[0], -tiny_vec2, [0, 0], tiny_vec2]) + bnd_pts[0]
				boundary_points[k] = np.append(bnd_pts, app_arr, axis=0)

			bnd_idx[:] = np.asarray(bnd_idx) % boundary_points[k].shape[0]

		boundary_indices = np.append(boundary_indices, [[-2]] * ns.size, axis=1)
		num_t_blocks = boundary_indices.shape[1]

	tangent_vecs = np.zeros((num_t_blocks, 3))
	offset_indices = np.zeros(num_t_blocks, dtype=np.integer)

	boundary_struct = CartBlockStruct(np.arange(num_t_blocks), ns, zs, nt, nn, nz, **kwargs)
	bl_edges = boundary_struct.edges
	bl_b_vts = boundary_struct.baked_vertices

	for k in range(zs.size):
		boundary_pts_z = np.asarray(boundary_points[k])
		num_b_pts = boundary_pts_z.shape[0]
		boundary_indices_z = np.asarray(boundary_indices[k]) % num_b_pts

		for i in range(num_t_blocks):
			b_idx = boundary_indices_z[i]
			tangent_vecs[i, :2] = boundary_pts_z[min(b_idx + 1, num_b_pts - 1)] - boundary_pts_z[max(b_idx - 1, 0)]

		unit_norm_vecs = np.cross(tangent_vecs, [0, 0, 1])[:, :2]
		unit_norm_vecs /= norm(unit_norm_vecs, axis=1)[:, np.newaxis]

		base_line_string = shp.LineString(boundary_pts_z)

		# Create offset curves
		for n in range(ns.size):
			if np.isclose(ns[n], 0.):
				bl_pts = boundary_pts_z
				offset_indices = boundary_indices_z
				boundary_struct.vertices[:, n, k, :2] = bl_pts[offset_indices]
			else:
				geometry_string = base_line_string.parallel_offset(ns[n], 'left')
				try:
					bl_line_string = geometry_string[1]
				except TypeError:
					bl_line_string = geometry_string

				bl_pts = np.array(bl_line_string)

				for i in range(boundary_indices_z.size):
					close_pt = boundary_pts_z[boundary_indices_z[i]] - unit_norm_vecs[i] * ns[n]

					# TEMPORARY CODE TO WORK AROUND BUG
					if i == 0 and boundary_indices_z[0] == 0:
						nearest_index = 0
					elif closed_boundary and i == num_t_blocks - 1:
						nearest_index = bl_pts.shape[0] - 2
					else:
						nearest_index = (norm(bl_pts - close_pt, axis=1)).argmin()

					boundary_struct.vertices[i, n, k, :2] = bl_pts[nearest_index]
					offset_indices[i] = nearest_index
					# close_pts[i] = close_pt

			for i in range(boundary_indices_z.size - 1):
				if offset_indices[i + 1] - offset_indices[i] > 1:
					pts = [Point(np.append(s_pt, zs[k])) for s_pt in bl_pts[offset_indices[i] + 1:offset_indices[i + 1]]]
					bl_edges[i, n, k, 0] = PolyLineCurvedEdge(bl_b_vts[i:i + 2, n, k], pts)

	if closed_boundary:
		boundary_struct.baked_vertices[-1] = boundary_struct.baked_vertices[0]

	return boundary_struct
