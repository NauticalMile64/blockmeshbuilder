
import numpy as np
from blockmeshbuilder import BlockMeshDict, CartBlockStruct, CylBlockStructContainer, BoundaryTag, ZoneTag

xs = ys = np.array([0.0, 1.0]) - 0.5
ndx = ndy = 30

cavity = CartBlockStruct(xs, ys, [0.0, 0.1], ndx, ndy, 4, zone_tag=ZoneTag('fluid_zone'))
pipe = CylBlockStructContainer((0.125, 0.25), np.linspace(0., 2 * np.pi, 9, endpoint=True), (0.1, 0.5), 5, 5, 5)

cavity_front_tag = BoundaryTag('cavity_front')
cavity.boundary_tags[:, :, -1, 2] = cavity_front_tag

pipe_entrance_tag = BoundaryTag('pipe_entrance')
pipe.core_struct.boundary_tags[:, :, 0, 2] = pipe_entrance_tag
pipe.tube_struct.boundary_tags[:, :, 0, 2] = pipe_entrance_tag

block_mesh_dict = BlockMeshDict(metric='mm')
cavity.write(block_mesh_dict)
pipe.write(block_mesh_dict)
block_mesh_dict.add_mergePatchPair(cavity_front_tag, pipe_entrance_tag)
block_mesh_dict.write_file('OF_case', run_blockMesh=True)
