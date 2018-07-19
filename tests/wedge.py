""" example of ofblockmeshdicthelper
Reproduce the wedge shown at
https://openfoamwiki.net/index.php/Main_ContribExamples/AxiSymmetric
"""

from math import radians,sin,cos
from ofblockmeshdicthelper import BlockMeshDict, Vertex, HexBlock, Boundary

# wedge dimensions
wd = radians(10.0)
r = 0.19
l = 1.1

# prepare ofblockmeshdicthelper.BlockMeshDict instance to
# gather vertices, blocks, faces and boundaries.
bmd = BlockMeshDict()

# set metrics
bmd.set_metric('m')

# Base points
b_vs = [Vertex(0.,0.,0.,'p0'),Vertex(r*cos(wd/2),r*sin(wd/2),0.,'p1'),
						Vertex(r*cos(-wd/2),r*sin(-wd/2),0.,'p2'),
		Vertex(0.,0.,l,'p3'),Vertex(r*cos(wd/2),r*sin(wd/2),l,'p4'),
						Vertex(r*cos(-wd/2),r*sin(-wd/2),l,'p5')]

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