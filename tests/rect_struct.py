# Creates a sample structured mesh in cartesian co-ordinates

import numpy as np
from ofblockmeshdicthelper import BlockMeshDict, CartBlockStruct, SimpleGradingElement, MultiGradingElement, \
	get_grading_info, Cylinder, Point

bmd = BlockMeshDict()  # Create a container to hold the objects

bmd.set_metric('mm')

# Create arrays of points for the co-ordinates
xs = np.linspace(0., 1., 4)
ys = np.linspace(0., 1., 3)
zs = np.linspace(0., 0.3, 3)

ndx = np.array([8, 12, 8, 0])
ndy = np.array([8, 8, 0])
ndz = np.array([6, 6, 0])

# Create the block structure
test_struct = CartBlockStruct(xs, ys, zs, ndx, ndy, ndz, zone='ts')

# The basic form of the structure is given by the definitions of xs,ys,zs: a rectilinear block-strucutred mesh. To make the final structure more interesting, we can use the power of numpy's slicing and indexing to make precise adjustments to the structure before it is written to the blockMeshDict.

# First we'll change the y-coordinates of some of the nodes on the top face of the block structure. The first three indices index the nodes, and the 4th index specifies which component (i.e. x, y, or z) is to be changed.

# Move the top nodes of the first block upward
test_struct['vertices'][:2, -1, :, 1] += 0.25

# The above statement applies to the first two nodes in the x-direction (:2), the last nodes in the y-direction (-1), all the nodes in the z-direction (:), and the final index selects the y-component of those nodes (1). The right hand expression adds the offset value to the sliced array.

# We can also apply this to other features of the mesh we want local control over, such as the grading
GD = test_struct['grading']

# All edges on the first two rows of blocks graded in the y-direction
GD[:, 0, :, 1] = SimpleGradingElement(1.0 / 3)
GD[:, 1, :, 1] = SimpleGradingElement(3.)

# All edges in the first and third columns of blocks graded in the x-direction
GD[0, :, :, 0] = SimpleGradingElement(1.0 / 3)
GD[2, :, :, 0] = SimpleGradingElement(3.)

# We can also add more complicated grading to the x-direction of the blocks in the second column

# Divide the x-edges of the block into 3 chunks: the first chunk is the first 20% of the length, the second is the next 60%, and the last chunk is the final 20%.
len_pcts = np.array([0.2, 0.6, 0.2])

# Now assign the grid densities at the boundaries of each of the chunks: at the beginning, the grid is densest with a relative value of 2.5. At the end of chunk 1 (20% of the block edge length), the grid density is reduced to 1. At the end of chunk 2 (80% of the length), the grid density is still 1, so the grid density is uniform through the second chunk. Finally we increase the density again at the right edge of the block so the grid is refined towards the boundary.
dens = np.array([2.5, 1., 1., 2.])

# We create a MulitGradingElement from this information, using the getGradingInfo helper function to translate this information into the length percent, cell percent, and expansion ratios expected by blockMesh
grd_elm = MultiGradingElement(*get_grading_info(len_pcts, dens))

# Assign the grading element to the blocks in the second column and the x-direction
GD[1, :, :, 0] = grd_elm

# Remove the block at the (1,0,0) index. Notice 3 indices are needed, this time since blocks don't have an explicit position.
test_struct['block_mask'][1, 0, :] = True

# Create a cylinder geometry along the right hand side of the block structure
vts = test_struct['vertices'][-1]
rad = (zs[-1] - zs[0]) / 2
pt1 = Point([xs[-1], -2, zs[1]])
pt2 = Point([xs[-1], 2, zs[1]])

cyl = Cylinder(pt1, pt2, rad, 'cyl')
bmd.add_geometry(cyl)

vts[:, 1, 0] += 0.1

# project the right side of the block structure onto the cylinder.
test_struct.project_structure(0, -1, cyl)

# Write the blocks to the blockMeshDict
test_struct.write(bmd)

# Write to file
with open(r'OF_case/system/blockMeshDict', 'w') as infile:
	infile.write(bmd.format())
