import numpy as np
from airfoil_points import NACA4
from blockmeshbuilder import BlockMeshDict, create_boundary_block_structure, BoundaryTag


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

# Create boundary layer struct
bl_ys_positive = np.array([0.05, 1.0, num_chords_downstream + 1])
bl_ys = np.concatenate((-bl_ys_positive[::-1], bl_ys_positive))
bl_xs = np.array([0., num_chords_bl_thickness]) - 0.5
zs = np.array([-0.5, 0.5]) * num_chords_span

bl_ndy_positive = [nx_chord, nx_downstream]
bl_ndy = np.array(bl_ndy_positive[::-1] + [n_nose] + bl_ndy_positive)
bl_ndx = n_bl
ndz = nz_span

# Get airfoil surface points using cosine spacing: http://airfoiltools.com/airfoil/naca4digit
x = (1 - np.cos(np.linspace(0, np.pi, spine_points))) / 2
xu, yu, xl, yl = airfoil.calculate_surface_points(x)
u_pts = np.dstack((xu, yu))
l_pts = np.dstack((xl, yl))
s_pts = np.hstack((l_pts[:, ::-1][:, :-1], u_pts))[0]

# Need to separate the lists of points into y, x arrays for gradient function
with np.errstate(divide='ignore', invalid='ignore'):
	sf_slope = np.nan_to_num(np.gradient(*s_pts[:, ::-1].T), copy=False, nan=np.inf)

# Find projection co-ordinate indices
leading_edge_slope = np.gradient(airfoil.calculate_camber_line([0., 0.001]), [0., 0.001])[0]
p_slope = np.tan(np.pi / 4 + np.arctan(leading_edge_slope))
m_slope = -np.tan(np.pi / 4 - np.arctan(leading_edge_slope))

u_nose_index = (np.abs(sf_slope - p_slope)).argmin()
l_nose_index = (np.abs(sf_slope - m_slope)).argmin()

sb_idx = [0, l_nose_index, u_nose_index, -1]

s_pts = [s_pts, s_pts]
sb_idx = np.tile(sb_idx, (2, 1))

nt = np.append(bl_ndy[1:-1], 20)

bl_struct = create_boundary_block_structure(s_pts, sb_idx, np.array((0.0, 0.02, 0.04)), zs, nt, 1, ndz, closed_boundary=True)

bl_struct.boundary_tags[:, 0, :, 1] = BoundaryTag('airfoilWall')

# Write the blocks to the blockMeshDict
block_mesh_dict = BlockMeshDict(metric='mm')
bl_struct.write(block_mesh_dict)
block_mesh_dict.write_file('OF_case', run_blockMesh=True)
