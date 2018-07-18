import numpy as np
from ofblockmeshdicthelper import BlockMeshDict, CartBlockStruct, SimpleGradingElement

xs = np.linspace(0.,1.,4)
ys = np.linspace(0.,1.,2)
zs = np.linspace(0.,0.1,2)

ndx = np.array([8,4,4,0])
ndy = np.array([8,0])
ndz = np.array([1,0])

test_struct = CartBlockStruct(xs,ys,zs,ndx,ndy,ndz,'ts')

#The first three indices index the nodes / blocks of interest, and the 4th index specifies which component (i.e. x,y, or z) is to be changed.

#Move the top nodes of the first two blocks upward
test_struct['vertices'][:3,-1,:,1] += 0.25

#The above statement applies to the first three nodes in the x-direction (:3), the last nodes in the y-direction (-1), all the nodes in the z-direction (:), and final index selects the y-component of those nodes (1). The right hand expression adds the offset value to the sliced array.

#We can also apply this to other features of the mesh we want local control over, such as the grading

#All blocks graded in the y-direction
test_struct['grading'][:,:,:,1] = SimpleGradingElement(1.0/3)

test_struct.bake_vertices() #Create the Vertex objects
test_struct.create_blocks() #Create the HexBlock objects

bmd = BlockMeshDict()		#Create a container to hold the objects
#Note that in more complicated cases we may also want to define blocks separately from the block_structure, and / or combine multiple block structures, so the blockmeshdict is needed to assemble everything.

test_struct.write_blocks(bmd) #Write the blocks to the blockMeshDict

#Write to file
with open('blockMeshDict','w') as infile:
	infile.write(bmd.format())