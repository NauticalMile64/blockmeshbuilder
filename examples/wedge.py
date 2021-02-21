"""
Reproduce the wedge shown at
https://openfoamwiki.net/index.php/Main_ContribExamples/AxiSymmetric
"""
from math import radians
from blockmeshbuilder import BlockMeshDict, TubeBlockStruct, BoundaryTag, ZoneTag

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

wedge = TubeBlockStruct(rs, ts, zs, nr, nt, nz, zone_tag=ZoneTag('wedge'))
wFaces = wedge['faces']

# Tag front and back boundaries
wedge['boundary_tags'][..., 0, 2] = BoundaryTag('front')
wedge['boundary_tags'][..., -1, 2] = BoundaryTag('back')

wedge.write(bmd)

# Output
with open(r'OF_case/system/blockMeshDict', 'w') as infile:
	infile.write(bmd.format())
