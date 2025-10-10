# Builds a structured O-grid mesh

import numpy as np
from blockmeshbuilder import BlockMeshDict, CartBlockStruct, TubeBlockStruct, BoundaryTag, Point, utilities

X = R = 0
Y = THETA = 1
Z = 2

channel_width = 6.
channel_height = 4.5
channel_length = 12.
partition_thickness = 0.65
partition_start = 0.75

b = partition_thickness
a = 1.2 * b
# a = b * 3 / 2

na = 2
nb = 3

dxs = np.full(na, a)
dys = np.full(nb, b)
xs = np.r_[0, np.cumsum(dxs)] - na * a / 2
ys = np.r_[0, np.cumsum(dys)] - nb * b / 2

rs = np.array([0.2, 0.35, 0.6, 1.0]) * channel_width

dls = np.tile(np.r_[dxs, dys], 2)
ls = np.cumsum(dls)
# ts = (np.r_[0, ls] * 2 / ls[-1] - 3 / 4) * np.pi

rab = np.sqrt((nb * b / 2)**2 + a**2)
ta = np.arcsin(a/rab)
tb = (np.pi - 2 * ta) / 3
dts = np.tile(np.r_[np.full(na, ta), np.full(nb, tb)], 2)
ts = np.cumsum(np.r_[0, dts]) - np.pi / 2 - ta
zs = np.array([0.0, channel_height])

nda = 5
ndx = np.full(na, nda)
ndy = np.full(nb, np.round(nda * a / b))
ndr = np.full_like(rs, 6)
ndt = np.tile(np.r_[ndx, ndy], 2)
ndz = 8

core = CartBlockStruct(xs, ys, zs, ndx, ndy, ndz)

iac = 0.45

tube = TubeBlockStruct(rs, ts, zs, ndr, ndt, ndz, zone_tag='ts', is_complete=True)

core_b_vts = core.baked_vertices
tube_b_vts = tube.baked_vertices

# Attach core and tube structures
tube_b_vts[0, [0, 1], :] = np.rot90(core_b_vts, k=0)[:-1, 0, :]
tube_b_vts[0, [2, 3, 4], :] = np.rot90(core_b_vts, k=-1)[:-1, 0, :]
tube_b_vts[0, [5, 6], :] = np.rot90(core_b_vts, k=-2)[:-1, 0, :]
tube_b_vts[0, [7, 8, 9], :] = np.rot90(core_b_vts, k=-3)[:-1, 0, :]
tube_b_vts[0, -1, :] = tube_b_vts[0, 0, :]

# Mask tube inner edges
tube.edge_mask[0] = True
tube.vertex_mask[0] = True

# Edge masking
tube.edge_mask[1, 3, :, THETA] = True
tube.edge_mask[2, 2:5, :, THETA] = True
tube.vertices[2, 3:5, :, R] -= 0.5

# Shrink fit
dummy_struct_xs = np.r_[partition_start, tube.vertices[1:3, 3, 0, R], channel_length]
dummy_struct_ys = np.array((-1, 1)) * partition_thickness / 2

dummy_struct = CartBlockStruct(dummy_struct_xs, dummy_struct_ys, zs, 1, 1, 1)

tube.baked_vertices[:, 3:5] = dummy_struct.baked_vertices[:]
core.baked_vertices[-1, 1:3] = dummy_struct.baked_vertices[0]

# Partition
pxs = [0, 2.5, channel_length]
pys_half = np.r_[partition_thickness / 2, rs[2:]]
pys = np.r_[-pys_half[::-1], pys_half]

npx = [ndt[5], 35]
npy_half = (ndt[4], ndr[-2])
npy = np.r_[npy_half[::-1], ndt[3], npy_half]
# npy = np.r_[ndr[2:][::-1], ndt[3], ndr[2:]]

baffle_loc = 0.6 * channel_width
bl_ind, bh_ind, py_ind = utilities.get_indices(pys, (-baffle_loc, baffle_loc, -partition_thickness / 2))

part = CartBlockStruct(pxs, pys, zs, npx, npy, ndz)

# Mate tube and part
tube.block_mask[:, 3, :] = True
tube.block_mask[2:, 1:6, :] = True
part.block_mask[:-2, bl_ind:bh_ind, :] = True
part.block_mask[:, py_ind] = True

# Mate outer bend
part.baked_vertices[0, -2:] = tube.baked_vertices[-2:, 6]
part.baked_vertices[0, :2] = tube.baked_vertices[-2:, 1][::-1]

# Mate inner bend
part.baked_vertices[1, 1:-1] = tube.baked_vertices[-2, 2:6]

# ##############
# Boundaries
# ##############

stationary_wall_tag = BoundaryTag('wallBoundaryStationary', type_='wall')
water_surface_bc_tag = BoundaryTag('waterFreeSurface', type_='wall')

part.boundary_tags[-1, :2, :, X] = BoundaryTag('inlet', type_='patch')
part.boundary_tags[-1, 3:, :, X] = BoundaryTag('outlet', type_='patch')

structs = [core, tube, part]
for struct in structs:
	struct.boundary_tags[:, :, -1, Z] = water_surface_bc_tag

block_mesh_dict = BlockMeshDict(metric='m')
for struct in structs:
	struct.write(block_mesh_dict)

block_mesh_dict.write_file('OF_case', run_blockMesh=True, density_scale=1.2, default_boundary_tag=BoundaryTag('wall'), block_structure_only=False)
