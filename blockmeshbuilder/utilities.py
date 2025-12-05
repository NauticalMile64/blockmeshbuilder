import numpy as np


def get_index(arr, entry):
	arr = np.asarray(arr)
	return np.where(np.isclose(arr, entry))[0][0]


def get_indices(arr, entries):
	arr = np.asarray(arr)
	return (get_index(arr, e) for e in entries)


# Experimental Interpolate Coordinates Function
def interp_coordinates(vertices, axis, C):
	C = np.asarray(C)
	init_pos = np.arange(vertices.ndim)
	roll_pos = np.r_[np.roll(init_pos[:-1], axis), vertices.ndim - 1]
	t_vts = np.moveaxis(vertices, init_pos, roll_pos)

	x_diff = t_vts[..., axis] - t_vts[0, np.newaxis, ..., axis]
	x_range_inv = 1. / (t_vts[-1, ..., axis] - t_vts[0, ..., axis])
	y_range = t_vts[-1, ..., C] - t_vts[0, ..., C]
	y_start = t_vts[0, np.newaxis, ..., C]
	t_vts[..., C] = np.einsum('i...,...,...->i...', x_diff, x_range_inv, y_range) + y_start


def number_of_divisions(arr_s, scale, weights=None):
	return np.round(np.maximum(np.multiply(np.diff(arr_s) * scale, weights or 1.0), 1)).astype(np.uint32)


def safe_int(i, default=0):
	try:
		return int(i)
	except ValueError:
		try:
			return int(''.join(filter(str.isdigit, s)))
		except:
			return default
