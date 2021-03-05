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
from blockmeshbuilder import BlockMeshDict, cart_to_cyl, cyl_to_cart, CartBlockStruct, \
	SimpleGradingElement, BoundaryTag, ZoneTag

xs = ys = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0]) - 0.5
zs = np.array([0.0, 0.01])

ndx = ndy = 14
ndz = 1

cavity = CartBlockStruct(xs, ys, zs, ndx, ndy, ndz, zone_tag=ZoneTag('fluid_zone'))

GD = cavity.grading
edge_grd = 4
GD[:, 0, :, 1] = SimpleGradingElement(edge_grd)  # bottom
GD[:, -2, :, 1] = SimpleGradingElement(1 / edge_grd)  # top
GD[0, :, :, 0] = SimpleGradingElement(edge_grd)  # left
GD[-2, :, :, 0] = SimpleGradingElement(1 / edge_grd)  # right

# Rotate vertices
vts_cyl = cart_to_cyl(cavity.vertices)  # Create new array of vertices in a cylindrical co-ordinate system

theta_coordinates = vts_cyl[..., 1]  # Select the angle co-ordinate theta
theta_coordinates[1:-1, 1:-1, :] += np.pi / 8  # Rotate interior blocks by 22.5 degrees
theta_coordinates[2:-2, 2:-2, :] += np.pi / 8  # Rotate innermost block by a further 22.5 degrees

cavity.vertices[:] = cyl_to_cart(vts_cyl)  # Transform back into the Cartesian system and re-assign

# Set middle block to solid
cavity.zone_tags[2, 2, 0] = ZoneTag('solid_zone')

# Label the lid
cavity.boundary_tags[:, -1, :, 1] = BoundaryTag('lid')

block_mesh_dict = BlockMeshDict(metric='mm')
cavity.write(block_mesh_dict)
block_mesh_dict.write_file('OF_case/system/')
