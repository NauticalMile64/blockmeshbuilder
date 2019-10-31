""" example of ofblockmeshdicthelper
Reproduce the wedge shown at
https://openfoamwiki.net/index.php/Main_ContribExamples/AxiSymmetric
"""
import numpy as np
from math import radians
from ofblockmeshdicthelper import BlockMeshDict, TubeBlockStruct, Boundary

# wedge dimensions
wd = radians(10.0)
r = 0.19
l = 1.1

# prepare ofblockmeshdicthelper.BlockMeshDict instance to gather vertices, blocks, faces and boundaries.
bmd = BlockMeshDict()

# set metrics
bmd.set_metric('mm')

rs = np.array([0,r])
ts = np.array([-wd/2,wd/2])
zs = np.array([0,1])

nr = np.array([10])
nt = np.array([1])
nz = np.array([10])

wedge = TubeBlockStruct(rs,ts,zs,nr,nt,nz,'wedge')
wFaces = wedge['faces']

#Front and back boundaries
front_bnd = Boundary('patch', 'front', faces = [wFaces[0,0,0,2]])
bmd.add_boundary(front_bnd)

back_bnd = Boundary('patch', 'back', faces = [wFaces[0,0,1,2]])
bmd.add_boundary(back_bnd)

wedge.write(bmd)

# output
with open(r'OF_case/system/blockMeshDict','w') as infile:
	infile.write(bmd.format())