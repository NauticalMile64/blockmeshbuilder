# Creates a sample structured mesh in cartesian co-ordinates
import numpy as np
from blockmeshbuilder import BlockMeshDict, CartBlockStruct, SimpleGradingElement, MultiGradingElement, \
	get_grading_info, Cylinder, Point, PlanePointAndNormal, BSplineCurvedEdge

# Create arrays of points for the co-ordinates
xs = np.linspace(0., 1., 4)
ys = np.linspace(0., 1., 3)
zs = np.linspace(0., 0.3, 3)

ndx = np.array([8, 12, 8])
ndy = 8
ndz = 6

# Create the block structure
rect_struct = CartBlockStruct(xs, ys, zs, ndx, ndy, ndz, zone_tag='test_zone')

# The basic form of the structure is given by the definitions of xs,ys,zs: a rectilinear block-structured mesh.
# To make the final structure more interesting, we can use the power of Numpy's slicing and indexing to make
# precise adjustments to the structure before it is written to the blockMeshDict.

# First we'll change the y-coordinates of some of the nodes on the top face of the block structure.
# The first three indices index the nodes, and the 4th index specifies which component (i.e. x, y, or z)
# is to be changed.

# Move the top nodes of the first block upward
rect_struct.vertices[:2, -1, :, 1] += 0.25

# The above statement applies to the first two nodes in the x-direction (:2), the last nodes in the y-direction (-1),
# all the nodes in the z-direction (:), and the final index selects the y-component of those nodes (1).
# The right hand expression adds the offset value to the sliced array.

# We can also apply this to other features of the mesh we want local control over, such as the grading:
GD = rect_struct.grading

# All edges on the first two rows of blocks graded in the y-direction
GD[:, 0, :, 1] = SimpleGradingElement(1.0 / 3)
GD[:, 1, :, 1] = SimpleGradingElement(3.)

# All edges in the first and third columns of blocks graded in the x-direction
GD[0, :, :, 0] = SimpleGradingElement(1.0 / 3)
GD[2, :, :, 0] = SimpleGradingElement(3.)

# We can also add more complicated grading to the x-direction of the blocks in the second column

# Divide the x-edges of the block into 3 chunks: the first chunk is the first 20% of the length,
# the second is the next 60%, and the last chunk is the final 20%.
len_pcts = np.array([0.2, 0.6, 0.2])

# Now assign the grid densities at the boundaries of each of the chunks:
# at the beginning, the grid is densest with a relative value of 2.5.
# At the end of chunk 1 (20% of the block edge length), the grid density is reduced to 1.
# At the end of chunk 2 (80% of the length), the grid density is still 1,
# so the grid density is uniform through the second chunk.
# Finally we increase the density again at the right edge of the block so the grid is refined towards the boundary.
dens = np.array([2.5, 1., 1., 2.])

# We create a MulitGradingElement from this information, using the getGradingInfo helper function to
# translate this information into the length percent, cell percent, and expansion ratios expected by blockMesh
grd_elm = MultiGradingElement(*get_grading_info(len_pcts, dens))

# Assign the grading element to the blocks in the second column and the x-direction
GD[1, :, :, 0] = grd_elm

# Remove the block at the (1,0,0) index. Notice 3 indices are needed
# this time since blocks don't have an explicit position.
rect_struct.block_mask[1, 0, :] = True

# Project the left side to a plane
plane_point = Point([0., 0, 0.])
plane_normal = Point([-1., -1., 0.])
plane = PlanePointAndNormal(plane_point, plane_normal, 'plane')

rect_struct.project_structure(0, plane)

# Create a spline edge 'awning' to top the 'doorway'
vertices = rect_struct.vertices
door_x_coordinates = vertices[1:3, 1, -1, 0]
door_x_values = np.linspace(door_x_coordinates[0], door_x_coordinates[1], 5, endpoint=True)
door_center = np.average(door_x_coordinates)

door_top_spline_coords = np.tile(vertices[1, 1, -1], (5, 1))
door_top_spline_coords[:, 0] = door_x_values
norm_x = (door_x_values - door_center) * 2 / np.diff(door_x_coordinates)[0]
door_top_spline_coords[:, 2] += 0.05 * (1.1 - norm_x**4)
door_top_spline_pts = [Point(c) for c in door_top_spline_coords]

door_top_baked_vts = rect_struct.edges[1, 1, -1, 0].vertices
rect_struct.edges[1, 1, -1, 0] = BSplineCurvedEdge(door_top_baked_vts, door_top_spline_pts)

# Create a cylinder geometry along the right hand side of the block structure
vts = rect_struct.vertices[-1]
rad = (zs[-1] - zs[0]) / 2
pt1 = Point([xs[-1], -2, zs[1]])
pt2 = Point([xs[-1], 2, zs[1]])

cyl = Cylinder(pt1, pt2, rad, 'cyl')

vts[:, 1, 0] += 0.1

# project the right side of the block structure onto the cylinder.
rect_struct.project_structure(0, cyl, -1)

# Write the blocks to the blockMeshDict
block_mesh_dict = BlockMeshDict(metric='mm')
rect_struct.write(block_mesh_dict)
block_mesh_dict.write_file('OF_case/system/')
