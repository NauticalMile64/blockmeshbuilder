# Builds a mesh with thin baffle structures

import numpy as np
from blockmeshbuilder import BlockMeshDict, cart_to_cyl, cyl_to_cart, CartBlockStruct, \
	SimpleGradingElement, BoundaryTag, ZoneTag

xs = ys = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0]) - 0.5
zs = np.array([0.0, 0.01])

ndx = ndy = 14
ndz = 1

cavity = CartBlockStruct(xs, ys, zs, ndx, ndy, ndz)

GD = cavity.grading
edge_grading = SimpleGradingElement(4)
GD[:, 0, :, 1] = edge_grading  # bottom
GD[:, -2, :, 1] = -edge_grading  # top
GD[0, :, :, 0] = edge_grading  # left
GD[-2, :, :, 0] = -edge_grading  # right

# # Rotate vertices
# vts_cyl = cart_to_cyl(cavity.vertices)  # Create new array of vertices in a cylindrical co-ordinate system
#
# theta_coordinates = vts_cyl[..., 1]  # Select the angle co-ordinate theta
# theta_coordinates[1:-1, 1:-1, :] += np.pi / 8  # Rotate interior blocks by 22.5 degrees
# theta_coordinates[2:-2, 2:-2, :] += np.pi / 8  # Rotate innermost block by a further 22.5 degrees
#
# cavity.vertices[:] = cyl_to_cart(vts_cyl)  # Transform back into the Cartesian system and re-assign
#
# # Set middle block to solid
# cavity.zone_tags[2, 2, 0] = ZoneTag('solid_zone')

# Label the lid
cavity.boundary_tags[:, -1, :, 1] = BoundaryTag('lid')

block_mesh_dict = BlockMeshDict(metric='mm')
cavity.write(block_mesh_dict)
block_mesh_dict.write_file('OF_case', run_blockMesh=True)
