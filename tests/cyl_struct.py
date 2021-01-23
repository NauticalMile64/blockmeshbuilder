# Builds a structured O-grid mesh

import numpy as np
from ofblockmeshdicthelper import BlockMeshDict, CylBlockStructContainer, Boundary

bmd = BlockMeshDict()
bmd.set_metric('mm')

rs = np.array([0.3, 0.6, 1.0])
ts = np.linspace(0, 2 * np.pi, 13, endpoint=True)  # Try any 4*n+1, where n is a positive int.
zs = np.array([0.0, 0.5, 1.5])

ndr = np.full_like(rs, 6)
ndt = np.full_like(ts, 6)
ndz = np.full_like(zs, 8)

iac = 0.001
cyl = CylBlockStructContainer(rs, ts, zs, ndr, ndt, ndz, zone='ts', is_core_aligned=True, inner_arc_curve=iac)

# Increase size of back half
scale = 1.15
cyl.tube_struct['vertices'][:, :, 1:, 0] *= scale		# Scale tube radii
cyl.core_struct['vertices'][:, :, 1:, [0, 1]] *= scale	# Scale core x and y co-ordinates

# Twist the outermost shell of the tube block structure
cyl.tube_struct['vertices'][-1, :-1, -1, 1] += np.pi / 16

wall_faces = cyl.tube_struct['faces'][-1, :-1, :-1, 0].flatten()
wall_bnd = Boundary('patch', 'wall', faces=wall_faces)
bmd.add_boundary(wall_bnd)

cyl.write(bmd)

with open(r'OF_case/system/blockMeshDict', 'w') as infile:
	infile.write(bmd.format())
