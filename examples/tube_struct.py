# Builds a structured tube mesh

import numpy as np
from blockmeshbuilder import BlockMeshDict, TubeBlockStruct, ZoneTag

bmd = BlockMeshDict(of_dist='.org')
bmd.set_metric('mm')

is_complete = True
rs = np.array([0., 0.5, 0.8, 1.0])
ts = np.linspace(0., 2 * np.pi, 9, endpoint=is_complete) - 3 * np.pi / 8
zs = np.linspace(0., 3., 6, endpoint=True)

ndr = np.full_like(rs, 6)
ndt = np.full_like(ts, 5)
ndz = np.full_like(zs, 8)

tube = TubeBlockStruct(rs, ts, zs, ndr, ndt, ndz, ZoneTag('ts'), is_complete=is_complete)

alternate_blocks = tube[1:, 1::2, 3:]
alternate_blocks['edge_mask'][..., 1:, 1] = True
alternate_blocks['face_mask'][..., 0] = True

tube['vertices'][-1, :-1, 0, 1] += np.pi / 16
tube['vertices'][-1, :-1, [2, 3], 1] -= np.pi / 16
tube['vertices'][-1, :-1, [2, 3], 0] += 0.2		# Produces correct conical surface for of_dist='.com'

tube.write(bmd)

with open(r'OF_case/system/blockMeshDict', 'w') as infile:
	infile.write(bmd.format())
