"""
Usage examples showing the Transform interface for CartBlockStruct and TubeBlockStruct

All transformations are visualizable:
- Translation: move by vector
- Rotation: axis-angle
- Scale: stretch/compress along axes
"""

from blockmeshbuilder import CartBlockStruct, BlockMeshDict, BoundaryTag
from blockmeshbuilder.transform import Transform, Quaternion
from blockmeshbuilder.transform import rotation_x, rotation_y, rotation_z, rotation_about_axis
import numpy as np
from math import radians

# ============================================================================
# Example 1: Simple translation
# ============================================================================

xs = [0, 1, 2]
ys = [0, 1]
zs = [0, 1]
nx, ny, nz = 10, 5, 5

# Translate along x-axis
transform = Transform(translation=[5, 0, 0])
block1 = CartBlockStruct(xs, ys, zs, nx, ny, nz, transform=transform, zone_tag="B1")
block1.boundary_tags[:] = BoundaryTag("B1")

# Can also use the old offset method
block = CartBlockStruct(xs, ys, zs, nx, ny, nz, offset=Point((5, 0, 0)))
block0.boundary_tags[:] = BoundaryTag("B0")

# ============================================================================
# Example 2: Pure rotation about Z-axis by 45 degrees
# ============================================================================

# Rotate about Z-axis (vertical)
angle = radians(45)
rotation = rotation_z(angle)
transform = Transform(rotation=rotation)
block2 = CartBlockStruct(xs, ys, zs, nx, ny, nz, transform=transform, zone_tag="B2")
block2.boundary_tags[:] = BoundaryTag("B2")

# ============================================================================
# Example 3: Rotate about arbitrary axis
# ============================================================================

# Rotate 30 degrees about axis pointing (1, 1, 0) direction
axis = [1, 1, 0]  # Will be normalized automatically
angle = radians(30)
rotation = rotation_about_axis(axis, angle)
transform = Transform(rotation=rotation)
block3 = CartBlockStruct(xs, ys, zs, nx, ny, nz, transform=transform)
block3.boundary_tags[:] = BoundaryTag("B3")

# ============================================================================
# Example 4: Combined rotation + translation
# ============================================================================

# First rotate 90° about X, then move up by 3 units
rotation = rotation_x(radians(90))
translation = [0, 0, 3]
transform = Transform(translation=translation, rotation=rotation)
block4 = CartBlockStruct(xs, ys, zs, nx, ny, nz, transform=transform)
block4.boundary_tags[:] = BoundaryTag("B4")

# ============================================================================
# Example 5: Align block with a direction vector
# ============================================================================

# Rotate block so its initial X-direction ([1,0,0]) aligns with [1,1,1]
from_vec = [1, 0, 0]  # Initial direction
to_vec = [1, 1, 1]  # Target direction (will be normalized)
rotation = Quaternion.align_vectors(from_vec, to_vec)
transform = Transform(rotation=rotation)
block5 = CartBlockStruct(xs, ys, zs, nx, ny, nz, transform=transform)
block5.boundary_tags[:] = BoundaryTag("B5")

# ============================================================================
# Example 6: Rotate + translate + scale (full transformation)
# ============================================================================

# Useful for creating repeated structures with variations
rotation = rotation_z(radians(60))
translation = [10, 5, 0]
scale = 1.5  # Uniform scaling by 1.5x
transform = Transform(translation=translation, rotation=rotation, scale=scale)
block6 = CartBlockStruct(xs, ys, zs, nx, ny, nz, transform=transform)
block6.boundary_tags[:] = BoundaryTag("B6")

# ============================================================================
# Example 7: Non-uniform scaling (anisotropic)
# ============================================================================

# Stretch 2x in X, compress 0.5x in Y, keep Z unchanged
scale = [2.0, 0.5, 1.0]
transform = Transform(scale=scale)
block7 = CartBlockStruct(xs, ys, zs, nx, ny, nz, transform=transform)
block7.boundary_tags[:] = BoundaryTag("B7")

# ============================================================================
# Example 8: Multiple rotations (quaternion composition)
# ============================================================================

# First rotate 45° about Z, then 30° about the new X-axis
rot1 = rotation_z(radians(45))
rot2 = rotation_x(radians(30))
combined_rotation = rot2 * rot1  # Apply rot1 first, then rot2
transform = Transform(rotation=combined_rotation)
block8 = CartBlockStruct(xs, ys, zs, nx, ny, nz, transform=transform)
block8.boundary_tags[:] = BoundaryTag("B8")

# ============================================================================
# Example 9: Create ring of blocks (circular pattern)
# ============================================================================

blocks = []
n_blocks = 8
radius = 5.0

for i in range(n_blocks):
	angle = 2 * np.pi * i / n_blocks

	# Rotate block to point outward
	rotation = rotation_z(angle)

	# Translate to circle
	translation = [radius * np.cos(angle), radius * np.sin(angle), 0]

	transform = Transform(translation=translation, rotation=rotation)
	block = CartBlockStruct(xs, ys, zs, nx, ny, nz, transform=transform)
	blocks.append(block)

# ============================================================================
# Example 10: Tube-like structure with wedge blocks
# ============================================================================

from blockmeshbuilder import TubeBlockStruct

rs = [0.5, 1.0]
thetas = [0, radians(30)]  # 30-degree wedge
zs = [0, 1]
nr, nt, nz = 5, 10, 20

# Create base wedge
wedge_base = TubeBlockStruct(rs, thetas, zs, nr, nt, nz)

# Create 12 wedges arranged in a circle
wedges = []
for i in range(12):
	angle = radians(30) * i  # Each wedge is 30°, so space by 30°
	rotation = rotation_z(angle)
	translation = [0, 0, 0.5 * i]
	transform = Transform(rotation=rotation, translation=translation)

	wedge = TubeBlockStruct(rs, thetas, zs, nr, nt, nz, transform=transform)
	wedges.append(wedge)

block_mesh_dict = BlockMeshDict(metric='m')

# block0.write(block_mesh_dict)
block1.write(block_mesh_dict)
block2.write(block_mesh_dict)
block3.write(block_mesh_dict)
block4.write(block_mesh_dict)
block5.write(block_mesh_dict)
block6.write(block_mesh_dict)
block7.write(block_mesh_dict)
block8.write(block_mesh_dict)

for block in blocks:
	block.write(block_mesh_dict)

for wedge in wedges:
	wedge.write(block_mesh_dict)

block_mesh_dict.write_file(r'OF_case', run_blockMesh=True)
