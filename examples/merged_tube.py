"""
This example demonstrates an instance where a `mergePatchPairs` operation can be done cleanly.
The outer tube circumferential block divisions do not align with the inner pipe, but because the total number of
divisions are made to match and the spacing is uniform, the merge faces are identical.
"""

import numpy as np
from blockmeshbuilder import BlockMeshDict, TubeBlockStruct, CylBlockStructContainer, BoundaryTag

pipe_radius = 0.5

rs_pipe = np.array([0.3, 0.6, 1.0]) * pipe_radius
num_side_blocks = 1  # Must be a positive integer
ts_pipe = np.linspace(0, 2 * np.pi, num_side_blocks * 4 + 1, endpoint=True)
zs_pipe = np.array([0.0, 0.1, 1])

pipe_divs = 10

rs_tube = np.array([1.0, 1.5, 2.0]) * pipe_radius
ts_tube = np.linspace(0, 2 * np.pi, pipe_divs + 1, endpoint=True)
zs_tube = zs_pipe[:-1]

ndr_pipe = 5
ndt_pipe = pipe_divs
ndz_pipe = np.array((3, 10))

ndr_tube = 5
ndt_tube = int(ndt_pipe * num_side_blocks * 4 / pipe_divs)
ndz_tube = ndz_pipe[0]

print(f'pipe divisions = {4 * num_side_blocks * ndt_pipe}')
print(f'tube divisions = {pipe_divs * ndt_tube}')

pipe = CylBlockStructContainer(rs_pipe, ts_pipe, zs_pipe, ndr_pipe, ndt_pipe, ndz_pipe, is_core_aligned=False)
tube = TubeBlockStruct(rs_tube, ts_tube, zs_tube, ndr_tube, ndt_tube, ndz_tube, is_complete=True)

tube.block_mask[1:, ::2] = True

tube_tag = BoundaryTag('tube_inner')
tube.boundary_tags[0, :, :, 0] = tube_tag

pipe_tag = BoundaryTag('pipe_outer')
pipe.tube_struct.boundary_tags[-1, :, 0, 0] = pipe_tag

end_tag = BoundaryTag('merged_end')
tube.boundary_tags[:, :, 0, 2] = end_tag
pipe.tube_struct.boundary_tags[:, :, 0, 2] = end_tag
pipe.core_struct.boundary_tags[:, :, 0, 2] = end_tag

block_mesh_dict = BlockMeshDict(metric='m')
tube.write(block_mesh_dict)
pipe.write(block_mesh_dict)
block_mesh_dict.add_mergePatchPair(tube_tag, pipe_tag)
block_mesh_dict.write_file('OF_case', run_blockMesh=True)
