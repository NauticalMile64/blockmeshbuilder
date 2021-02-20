import numpy as np
from airfoil_points import NACA4
from blockmeshbuilder import BlockMeshDict, CartBlockStruct, Point, \
	BSplineCurvedEdge, PolyLineCurvedEdge, SimpleGradingElement, MultiGradingElement, get_grading_info, BoundaryTag
import shapely.geometry as shp

# Airfoil shape and dimensions
close_trailing_edge = False
airfoil = NACA4((4, 4, 24), close_trailing_edge=close_trailing_edge)
chord_length = 1.0
num_chords_upstream = 1.5
num_chords_downstream = 2.0
num_chords_normal = 1.5
num_chords_bl_thickness = 0.05
num_chords_span = 0.1

# Node counts
nx_upstream = np.round(60*num_chords_upstream)
nx_chord = 100
nx_downstream = np.round(75*num_chords_downstream)
ny_wall_normal = np.round(75*num_chords_downstream)
nz_span = np.round(10*num_chords_span)
n_bl = np.round(350*num_chords_bl_thickness)
n_nose = 22

# Spline Resolution
spine_points = 100
spline_res = 1

# Create the blockMeshDict object
block_mesh_dict = BlockMeshDict()
block_mesh_dict.set_metric('mm')

# Create boundary layer struct
bl_ys_positive = np.array([0.05, 1.0, num_chords_downstream + 1])
bl_ys = np.concatenate((-bl_ys_positive[::-1], bl_ys_positive))
bl_xs = np.array([0., num_chords_bl_thickness]) - 0.5
zs = np.array([-0.5, 0.5]) * num_chords_span

bl_ndy_positive = [nx_chord, nx_downstream]
bl_ndy = np.array(bl_ndy_positive[::-1] + [n_nose] + bl_ndy_positive + [0])
bl_ndx = np.array([n_bl, 0])
ndz = np.array([nz_span, 0])

bl_struct = CartBlockStruct(bl_xs, bl_ys, zs, bl_ndx, bl_ndy, ndz, zone='boundary_layer')

# Get airfoil surface points using cosine spacing: http://airfoiltools.com/airfoil/naca4digit
x = (1 - np.cos(np.linspace(0, np.pi, spine_points))) / 2
xu, yu, xl, yl = airfoil.calculate_surface_points(x)
u_pts = np.dstack((xu, yu))
l_pts = np.dstack((xl, yl))
s_pts = np.hstack((l_pts[:, ::-1], u_pts))[0]

# Get the boundary layer outer profile
af_line_string = shp.LineString(s_pts)
bl_pts = np.array(af_line_string.parallel_offset(num_chords_bl_thickness, 'left'))

# Need to separate the lists of points into y, x arrays for gradient function
with np.errstate(divide='ignore', invalid='ignore'):
	sf_slope = np.nan_to_num(np.gradient(*s_pts[:, ::-1].T), copy=False, nan=np.inf)
	bl_slope = np.nan_to_num(np.gradient(*bl_pts[:, ::-1].T), copy=False, nan=np.inf)

# Find projection co-ordinate indices
leading_edge_slope = np.gradient(airfoil.calculate_camber_line([0., 0.001]), [0., 0.001])[0]
p_slope = np.tan(np.pi / 4 + np.arctan(leading_edge_slope))
m_slope = -np.tan(np.pi / 4 - np.arctan(leading_edge_slope))

u_nose_index = (np.abs(sf_slope - p_slope)).argmin()
l_nose_index = (np.abs(sf_slope - m_slope)).argmin()
ub_nose_index = (np.abs(bl_slope - p_slope)).argmin()
lb_nose_index = (np.abs(bl_slope - m_slope)).argmin()

sb_idx = [0, l_nose_index, u_nose_index, -1]
bb_idx = [0, lb_nose_index, ub_nose_index, -1]

bl_struct['vertices'][-1, 1:5, :, :2] = s_pts[sb_idx][:, np.newaxis]
bl_struct['vertices'][0, 1:5, :, :2] = bl_pts[bb_idx][:, np.newaxis]

bl_edges = bl_struct['edges']
bl_b_vts = bl_struct['baked_vertices']

for j in range(1, 4):
	for k in range(zs.size):
		pts = [Point(np.append(b_pt, zs[k])) for b_pt in bl_pts[bb_idx[j - 1] + 1:bb_idx[j]:spline_res]]
		bl_edges[0, j, k, 1] = PolyLineCurvedEdge(bl_b_vts[0, j:j + 2, k], pts)

		pts = [Point(np.append(s_pt, zs[k])) for s_pt in s_pts[sb_idx[j - 1] + 1:sb_idx[j]:spline_res]]
		bl_edges[-1, j, k, 1] = PolyLineCurvedEdge(bl_b_vts[-1, j:j + 2, k], pts)

# Align downstream boundary layer block structure
bl_struct['vertices'][:, [0, -1], :, 0] = 1 + num_chords_downstream
bl_struct['vertices'][0, 0, :, 1] = -num_chords_bl_thickness
bl_struct['vertices'][-1, [0, -1], :, 1] = bl_struct['vertices'][-1, [1, -2], :, 1]
bl_struct['vertices'][0, -1, :, 1] = num_chords_bl_thickness

te_struct = None
if close_trailing_edge:
	bl_struct['baked_vertices'][-1, [-1, -2]] = bl_struct['baked_vertices'][-1, [0, 1]]
else:
	te_struct = CartBlockStruct([0, 1], [0, 1], zs, [3, 0], [nx_downstream, 0], ndz)
	te_struct['baked_vertices'][0] = bl_struct['baked_vertices'][-1, :2]
	te_struct['baked_vertices'][-1] = bl_struct['baked_vertices'][-1, -2:][::-1]

# Create the far-field structure
ff_xs = np.array([-num_chords_upstream, 0., 1., 1 + num_chords_downstream])
ff_ys_positive = np.array([num_chords_bl_thickness, num_chords_normal])
ff_ys = np.concatenate((-ff_ys_positive[::-1], ff_ys_positive))

ff_ndx = np.array([nx_upstream, nx_chord, nx_downstream, 0])
ff_ndy = np.array([ny_wall_normal, n_nose, ny_wall_normal, 0])
ff_struct = CartBlockStruct(ff_xs, ff_ys, zs, ff_ndx, ff_ndy, ndz, zone='far_field')

# Cut out the center of the far-field struct for the boundary layer struct to fit
ff_struct['block_mask'][1:, 1, :] = True

# Mate far-field structure to boundary layer structure
ff_struct['baked_vertices'][1:, 1] = bl_struct['baked_vertices'][0, 2::-1]
ff_struct['baked_vertices'][1:, 2] = bl_struct['baked_vertices'][0, 3:]

# Mesh tuning
bl_struct['grading'][0, ..., 0] = SimpleGradingElement(1./2)

len_pcts = np.array([0.8, 0.2])
dens = np.array([1., 1., 3.])
ff_struct['grading'][0, ..., 0] = MultiGradingElement(*get_grading_info(len_pcts, dens))
ff_struct['grading'][:, 0, :, 1] = MultiGradingElement(*get_grading_info(len_pcts, dens))
ff_struct['grading'][:, 2, :, 1] = MultiGradingElement(*get_grading_info(len_pcts[::-1], dens[::-1]))

# Define boundary_tags
ff_struct['boundary_tags'][0, ..., 0] = BoundaryTag('inlet')

outlet_tag = BoundaryTag('outlet')
ff_struct['boundary_tags'][-1, [0, -2], :, 0] = outlet_tag
bl_struct['boundary_tags'][0, [0, -1], :, 1] = outlet_tag
if te_struct:
	te_struct['boundary_tags'][:, 0, :, 1] = outlet_tag

bl_struct['boundary_tags'][-1, 1:-2, :, 0] = BoundaryTag('wall-airfoil')

# Write the blocks to the blockMeshDict
bl_struct.write(block_mesh_dict)
ff_struct.write(block_mesh_dict)
if te_struct:
	te_struct.write(block_mesh_dict)

# Write to file
with open(r'OF_case/system/blockMeshDict', 'w') as infile:
	infile.write(block_mesh_dict.format())

try:
	import matplotlib.pyplot as plt
	plt.plot(*s_pts.T)
	plt.plot(*bl_pts.T)
	plt.scatter(*s_pts[[u_nose_index, l_nose_index]].T)
	plt.scatter(*bl_pts[[ub_nose_index, lb_nose_index]].T)
	ax = plt.gca()
	ax.set_aspect('equal')
	plt.show()
except ImportError:
	pass
