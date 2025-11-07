"""
This example generates a blockMeshDict to build a square cavity mesh which can be used for running a lid-driven-cavity
CFD problem.

Two variations on the traditional lid-driven cavity have been implemented here:
1. The inner blocks have been twisted, so the mesh is no longer rectilinear
2. The center block has been assigned to a different zone (representing a solid square).
These zones are not automatically incorporated into the final mesh using the blockMesh command.
Therefore an additional step is required: after meshing, type the command `splitMeshRegions -cellZones -overwrite`.

If you wish to simply omit the center block from the mesh instead of writing it, change the 'block_mask' of the 2,2,0
block to True instead of editing the 'zone_tags' array.
"""
import numpy as np
from blockmeshbuilder import BlockMeshDict, CartBlockStruct, \
	SimpleGradingElement, BoundaryTag, ZoneTag
from blockmeshbuilder.transform import Transform, rotation_z

xs = ys = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0]) - 0.5
zs = np.array([0.0, 0.01])

ndx = ndy = 14
ndz = 1

cavity = CartBlockStruct(xs, ys, zs, ndx, ndy, ndz, zone_tag=ZoneTag('fluid_zone'))

GD = cavity.grading
edge_grading = SimpleGradingElement(4)
GD[:, 0, :, 1] = edge_grading  # bottom
GD[:, -2, :, 1] = -edge_grading  # top
GD[0, :, :, 0] = edge_grading  # left
GD[-2, :, :, 0] = -edge_grading  # right

# Successively twist smaller groups of blocks about the z-axis
# Note the innermost block has the transformation applied twice, for a 45 deg twist
transform = Transform(rotation=rotation_z(np.radians(22.5)))
vts_twist = cavity.vertices
for i in range(2):
	vts_twist = vts_twist[1:-1, 1:-1, :]
	vts_twist[:] = transform.apply(vts_twist)

# Set middle block to solid
cavity.zone_tags[2, 2, 0] = ZoneTag('solid_zone')

# Label the lid
cavity.boundary_tags[:, -1, :, 1] = BoundaryTag('lid')

block_mesh_dict = BlockMeshDict(metric='mm')
cavity.write(block_mesh_dict)
block_mesh_dict.write_file('OF_case', run_blockMesh=True)
