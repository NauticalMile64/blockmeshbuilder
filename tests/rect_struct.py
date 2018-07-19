#Creates a sample structured mesh in cartesian co-ordinates

import numpy as np
from ofblockmeshdicthelper import BlockMeshDict, CartBlockStruct, SimpleGradingElement

bmd = BlockMeshDict()		#Create a container to hold the objects

#Create arrays of points for the co-ordinates
xs = np.linspace(0.,1.,4)
ys = np.linspace(0.,1.,3)
zs = np.linspace(0.,0.1,2)

ndx = np.array([8,6,8,0])
ndy = np.array([8,8,0])
ndz = np.array([1,0])

#Create the block structure
test_struct = CartBlockStruct(xs,ys,zs,ndx,ndy,ndz,'ts')

#The basic form of the structure is given by the definitions of xs,ys,zs: a rectilinear block-strucutred mesh. To make the final structure more interesting, we can use the power of numpy's slicing and indexing to make precise adjustments to the structure before it is written to the blockMeshDict.

#First we'll change the y-coordinates of some of the nodes on the top face of the block structure. The first three indices index the nodes, and the 4th index specifies which component (i.e. x, y, or z) is to be changed.

#Move the top nodes of the first two blocks upward
test_struct['vertices'][:3,-1,:,1] += 0.25

#The above statement applies to the first three nodes in the x-direction (:3), the last nodes in the y-direction (-1), all the nodes in the z-direction (:), and the final index selects the y-component of those nodes (1). The right hand expression adds the offset value to the sliced array.

#We can also apply this to other features of the mesh we want local control over, such as the grading
GD = test_struct['grading']

#All edges on the first two rows of blocks graded in the y-direction
GD[:,0,:,1] = SimpleGradingElement(1.0/3)
GD[:,1,:,1] = SimpleGradingElement(3.)

#All edges in the first and third columns of blocks graded in the x-direction
GD[0,:,:,0] = SimpleGradingElement(1.0/3)
GD[2,:,:,0] = SimpleGradingElement(3.)

#Remove the block at the (1,0,0) index. Notice 3 indices are needed, this time since blocks don't have an explicit position.
test_struct['block_mask'][1,0,0] = True

#The modifications are finished, so now the objects can be created.

test_struct.bake_vertices() #Create the Vertex objects
test_struct.create_blocks() #Create the HexBlock objects

#Note that in more complicated cases we may also want to define blocks separately from the block_structure, and / or combine multiple block structures, so the blockmeshdict is needed to assemble everything.

test_struct.write_blocks(bmd) #Write the blocks to the blockMeshDict

#Write to file
with open(r'OF_case/system/blockMeshDict','w') as infile:
	infile.write(bmd.format())