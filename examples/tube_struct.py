# Builds a structured tube mesh

import numpy as np
from blockmeshbuilder import BlockMeshDict, TubeBlockStruct

is_complete = True
rs = np.array([0., 0.5, 0.8, 1.0])
ts = np.linspace(0., 2 * np.pi, 9, endpoint=is_complete) - 3 * np.pi / 8
zs = np.linspace(0., 3., 6, endpoint=True)

ndr = 6
ndt = 5
ndz = 8

tube = TubeBlockStruct(rs, ts, zs, ndr, ndt, ndz, zone_tag='tube', is_complete=is_complete)

alternate_blocks = tube[1:, 1::2, 3:]
alternate_blocks.edge_mask[..., 1:, 1] = True
alternate_blocks.face_mask[..., 0] = True

tube.vertices[-1, :-1, 0, 1] += np.pi / 16
tube.vertices[-1, :-1, [2, 3], 1] -= np.pi / 16
tube.vertices[-1, :-1, [2, 3], 0] += 0.2		# Produces correct conical surface for of_dist='.com'

block_mesh_dict = BlockMeshDict(metric='mm', of_dist='.org')
tube.write(block_mesh_dict)
block_mesh_dict.write_file('OF_case/system/')
