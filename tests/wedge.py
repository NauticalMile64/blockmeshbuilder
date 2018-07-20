""" example of ofblockmeshdicthelper
Reproduce the wedge shown at
https://openfoamwiki.net/index.php/Main_ContribExamples/AxiSymmetric
"""

from math import radians
from ofblockmeshdicthelper import BlockMeshDict, Vertex, HexBlock, Boundary, cyl_to_cart

# wedge dimensions
wd = radians(10.0)
r = 0.19
l = 1.1

# prepare ofblockmeshdicthelper.BlockMeshDict instance to
# gather vertices, blocks, faces and boundaries.
bmd = BlockMeshDict()

# set metrics
bmd.set_metric('mm')

# Base points
b_vs = [Vertex([0.,0.,0.],cyl_to_cart,'p0'),
		Vertex([r,wd/2,0.],cyl_to_cart,name = 'p1'),
		Vertex([r,-wd/2,0.],cyl_to_cart,name = 'p2'),
		Vertex([0.,0.,l],cyl_to_cart,name = 'p3'),
		Vertex([r,wd/2,l],cyl_to_cart,name = 'p4'),
		Vertex([r,-wd/2,l],cyl_to_cart,name = 'p5')]

#Number of divisions along each direction
nds = (10, 1, 10)

block_name = 'wedge_block'
block = HexBlock((b_vs[0],b_vs[2],b_vs[1],b_vs[0],
					b_vs[3],b_vs[5],b_vs[4],b_vs[3]),
					nds,block_name)
bmd.add_hexblock(block,block_name)

#Front and back boundaries
front_bnd = Boundary('patch', 'front', faces = [block.face('s')])
bmd.add_boundary(front_bnd)

back_bnd = Boundary('patch', 'back', faces = [block.face('n')])
bmd.add_boundary(back_bnd)

# output
with open(r'OF_case/system/blockMeshDict','w') as infile:
	infile.write(bmd.format())