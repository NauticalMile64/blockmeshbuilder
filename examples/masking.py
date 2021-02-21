# Creates a sample structured mesh in cartesian co-ordinates

import numpy as np
from blockmeshbuilder import BlockMeshDict, CartBlockStruct, ZoneTag

bmd = BlockMeshDict()  # Create a container to hold the objects

bmd.set_metric('mm')

# Create arrays of points for the co-ordinates
xs = np.linspace(0., 1., 4)
ys = np.linspace(0., 1 + 0.5, 6)
zs = np.linspace(0., 0.3, 3)

nElms = 8
ndx = np.array([8, 8, 8, 0])
ndy = np.array([8, 8, 8, 8, 8, 0])
ndz = np.array([8, 8, 0])

# Create the block structure
f_struct = CartBlockStruct(xs, ys, zs, ndx, ndy, ndz, zone_tag=ZoneTag('ts'))

# Mask structure
f_struct['block_mask'][:-1, :-1, :] = np.array(
	[[0, 0, 0],
	 [0, 1, 1],
	 [0, 0, 0],
	 [0, 1, 1],
	 [0, 1, 1]]).T[:, ::-1, np.newaxis]

ixs = np.linspace(1., 1.5, 2) + 0.1
ind_x = np.array([8, 0])
i_struct = CartBlockStruct(ixs, ys + 0.5, zs, ind_x, ndy, ndz, zone_tag=ZoneTag('ts'))

# Align the x and y co-ordinates of the relevant i_struct blocks with the f_struct prior to mating
i_struct['vertices'][..., 0] -= 0.1  # x co-ordinates
i_struct['vertices'][:, :-2, :, 1] = f_struct['vertices'][-1, 2:, :, 1][np.newaxis, ...]  # y co-ordinates

# Mate the structures by setting the relevant vertex objects (or 'baked_vertices') equal to one another
i_struct['baked_vertices'][0, :-2, :] = f_struct['baked_vertices'][-1, 2:, :]

# When blockmeshbuilder writes this blockMeshDict, it will use the same vertices, and therefore will identify
# common edges and faces of the mated blocks, and treat them as a contiguous structure.

# After the blocks are mated, the mated vertices are only modifiable through the f_struct.
f_struct['vertices'][-1, 2:, :, 1] += 0.1  # This line moves the mated vertices up
i_struct['vertices'][0, :-2, :, 1] -= 0.2  # This line does nothing

# Write the blocks to the blockMeshDict
f_struct.write(bmd)
i_struct.write(bmd)

# Write to file
with open(r'OF_case/system/blockMeshDict', 'w') as infile:
	infile.write(bmd.format())
