# Builds a structured O-grid mesh

import numpy as np
from blockmeshbuilder import BlockMeshDict, CylBlockStructContainer, BoundaryTag, ZoneTag

rs = np.array([0.3, 0.6, 1.0])
ts = np.linspace(0, 2 * np.pi, 13, endpoint=True)  # Try any 4*n+1, where n is a positive int.
zs = np.array([0.0, 0.5, 1.5])

ndr = np.full_like(rs, 6)
ndt = np.full_like(ts, 6)
ndz = np.full_like(zs, 8)

iac = 0.45
cyl = CylBlockStructContainer(rs, ts, zs, ndr, ndt, ndz, zone_tag=ZoneTag('ts'), is_core_aligned=True, inner_arc_curve=iac)

# Increase size of back half
scale = 1.15
cyl.tube_struct['vertices'][:, :, 1:, 0] *= scale		# Scale tube radii
cyl.core_struct['vertices'][:, :, 1:, [0, 1]] *= scale	# Scale core x and y co-ordinates

# Twist the outermost shell of the tube block structure
cyl.tube_struct['vertices'][-1, :-1, -1, 1] += np.pi / 16

cyl.tube_struct['boundary_tags'][-1, ..., 0] = BoundaryTag('wall')

block_mesh_dict = BlockMeshDict(metric='mm')
cyl.write(block_mesh_dict)
block_mesh_dict.write_file('OF_case/system/')
