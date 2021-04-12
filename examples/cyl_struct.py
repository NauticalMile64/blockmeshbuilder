# Builds a structured O-grid mesh

import numpy as np
from blockmeshbuilder import BlockMeshDict, CylBlockStructContainer, BoundaryTag, Point

rs = np.array([0.3, 0.6, 1.0])
num_side_blocks = 3  # Must be a positive integer
ts = np.linspace(0, 2 * np.pi, num_side_blocks * 4 + 1, endpoint=True)
zs = np.array([0.0, 0.5, 1.5])

ndr = 6

# Specify node counts of circumferential blocks; as with other structures len(ndt) == len(ts) or len(ndt) = len(ts) - 1
nd1 = [2, 3, 4]  # Node counts for successive blocks in x direction
nd2 = [5, 6, 7]  # Node counts for successive blocks in y direction
ndt = nd1 + nd2 + nd1[::-1] + nd2[::-1]  # Specifying the node counts in this way preserves the o-grid anti-symmetry

ndz = 8

iac = 0.45
cyl = CylBlockStructContainer(rs, ts, zs, ndr, ndt, ndz, zone_tag='ts', is_core_aligned=True, inner_arc_curve=iac,
							  offset=Point((1., 1., 0.)))

# Increase size of back half
scale = 1.15
cyl.tube_struct.vertices[:, :, 1:, 0] *= scale		# Scale tube radii
cyl.core_struct.vertices[:, :, 1:, [0, 1]] *= scale	# Scale core x and y co-ordinates

# Twist the outermost shell of the tube block structure
cyl.tube_struct.vertices[-1, :-1, -1, 1] += np.pi / 16

cyl.tube_struct.boundary_tags[-1, ..., 0] = BoundaryTag('wall')

block_mesh_dict = BlockMeshDict(metric='mm')
cyl.write(block_mesh_dict)
block_mesh_dict.write_file('OF_case', run_blockMesh=True)
