"""
Reproduce the wedge shown at
https://openfoamwiki.net/index.php/Main_ContribExamples/AxiSymmetric
"""
from math import radians
from blockmeshbuilder import BlockMeshDict, TubeBlockStruct, Boundary

# Wedge dimensions
wd = radians(10.0)
r = 0.19

# Prepare blockmeshbuilder to gather projection geometries, boundary faces, and block structures.
bmd = BlockMeshDict()

# Set metrics
bmd.set_metric('mm')

rs = [0, r]
ts = [-wd / 2, wd / 2]
zs = [0, 1]

nr = [10]
nt = [1]
nz = [10]

wedge = TubeBlockStruct(rs, ts, zs, nr, nt, nz, 'wedge')
wFaces = wedge['faces']

# Front and back boundaries
front_bnd = Boundary('patch', 'front', faces=[wFaces[0, 0, 0, 2]])
bmd.add_boundary(front_bnd)

back_bnd = Boundary('patch', 'back', faces=[wFaces[0, 0, 1, 2]])
bmd.add_boundary(back_bnd)

wedge.write(bmd)

# Output
with open(r'OF_case/system/blockMeshDict', 'w') as infile:
	infile.write(bmd.format())
