"""
Reproduce the wedge shown at
https://openfoamwiki.net/index.php/Main_ContribExamples/AxiSymmetric
"""
from math import radians
from blockmeshbuilder import BlockMeshDict, TubeBlockStruct, BoundaryTag

# Wedge dimensions
wedge_angle = radians(10.0)
wedge_radius = 0.19

# Co-ordinates of block vertices (r, theta, z)
rs = [0, wedge_radius]
ts = [-wedge_angle / 2, wedge_angle / 2]
zs = [0, 1]

# Number of cell divisions of the block in each direction
nr = 10
nt = 1
nz = 10

# Create the block structure
wedge = TubeBlockStruct(rs, ts, zs, nr, nt, nz, zone_tag='wedge')

# Tag front and back boundaries
wedge.boundary_tags[..., 0, 2] = BoundaryTag('front')
wedge.boundary_tags[..., -1, 2] = BoundaryTag('back')

# Initialize blockmeshbuilder to gather block structures.
block_mesh_dict = BlockMeshDict(metric='mm', of_dist='.org')

# Write the wedge to the block_mesh_dict object
wedge.write(block_mesh_dict)

# Output
block_mesh_dict.write_file('OF_case')
