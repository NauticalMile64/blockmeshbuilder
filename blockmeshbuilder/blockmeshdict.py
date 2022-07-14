from io import StringIO
import warnings
import subprocess
from pathlib import Path
from .blockelements import Face
from .geometry import _of_distribution_geometries
from .boundary_tags import BoundaryTag, _Boundary, _of_distribution_constraints
from .merge_patch_pairs import _MergePatchPair
from .check_names import BoundaryNameClashError
from .version import __version__
import numpy as np

init_pos = np.arange(3)
init_pos.setflags(write=False)

_of_distributions = ('.org', '.com')
features = ('constraints', 'geometries')
_of_distribution_features = dict()
for dist in _of_distributions:
	_of_distribution_features[dist] = dict()
	for dict_name, feature_dict in zip(features, (_of_distribution_constraints, _of_distribution_geometries)):
		_of_distribution_features[dist][dict_name] = feature_dict[dist]


def _format_section(section_name, section_items):
	buf = StringIO()
	if section_name == 'geometry':
		brackets = ['{', '}']
	else:
		brackets = ['(', ')']

	buf.write(f'{section_name}\n')
	buf.write(f'{brackets[0]}\n')
	for item in section_items:
		buf.write(f'    {item.format()}\n')
	buf.write(f'{brackets[1]};')

	return buf.getvalue()


class BlockMeshDict:
	metric_conversion_dict = {
		'km': 1000,
		'm': 1,
		'cm': 0.01,
		'mm': 0.001,
		'um': 1e-6,
		'nm': 1e-9,
		'A': 1e-10,
		'Angstrom': 1e-10
	}

	def __init__(self, metric='m', of_dist='.org'):
		if of_dist not in _of_distributions:
			warnings.warn(
				f'Unknown OpenFOAM distribution {of_dist}. The available options are {_of_distributions}. '
				f'\nSwitching to .org distribution.')
			of_dist = '.org'
		self.of_distribution = of_dist
		self.of_distribution_features = _of_distribution_features[of_dist]
		self.convert_to_meters = self.metric_conversion_dict[metric]
		self.blocks = set()
		self.edges = set()
		self.boundaries = {}
		self.geometries = set()
		self.faces = set()
		self.mergePatchPairs = set()

	def add_hexblock(self, block):
		self.blocks.add(block)

	def add_edge(self, edge):
		self.edges.add(edge)

	def add_boundary_face(self, boundary_tag, face):

		if boundary_tag not in self.boundaries:
			for other_boundary_tag in self.boundaries:
				if boundary_tag.name == other_boundary_tag.name:
					raise BoundaryNameClashError(f'The boundary named {boundary_tag.name} of type '
												 f'{boundary_tag._type} shares a name with another boundary of type '
												 f'{other_boundary_tag._type}. Please rename one of the boundaries.')

			self.boundaries[boundary_tag] = _Boundary(boundary_tag)

		self.boundaries[boundary_tag].add_face(face)

	def add_geometries(self, other_geometries):
		for geometry in other_geometries:
			if type(geometry) not in self.of_distribution_features['geometries']:
				raise TypeError(
					f'Geometry of type {type(geometry)} is not implemented in the {self.of_distribution} distribution.')

		self.geometries.update(other_geometries)

	def add_face(self, face):
		self.faces.add(face)

	def add_mergePatchPair(self, large_boundary, small_boundary):
		self.mergePatchPairs.add(_MergePatchPair(large_boundary, small_boundary))

	def _assign_vertexid(self):
		valid_vertices = []

		i = 0
		for b in self.blocks:
			for v in b.vertices.flatten():
				if v not in valid_vertices:
					valid_vertices.append(v)
					v.index = i
					i += 1

		self.valid_vertices = valid_vertices

	def format(self, block_structure_only=False, density_scale=1.0, default_boundary_tag=None):

		if block_structure_only:
			for block in self.blocks:
				block.cells = (1, 1, 1)
		elif density_scale != 1.0:
			for block in self.blocks:
				block.cells = tuple(int(max(round(block.cells[i] * density_scale), 1)) for i in range(3))

		if default_boundary_tag is not None:
			if not isinstance(default_boundary_tag, BoundaryTag):
				raise TypeError('default_boundary must be a BoundaryTag.')

			# Get list of boundary faces already in use
			face_vertices_dict = dict()
			for boundary in self.boundaries.values():
				for face in boundary.faces:
					face_vertices_dict[face] = sorted(face.vertices.flatten(), key=id)

			# Find external faces not assigned to boundaries
			for b in self.blocks:
				block_remaining_faces = dict()
				for s in range(3):
					block_vertices = np.moveaxis(b.vertices, init_pos, np.roll(init_pos, s))
					for face_vertices in block_vertices:
						face_vertices_sorted = sorted(face_vertices.flatten(), key=id)
						match_found = False
						for face, vertices in face_vertices_dict.items():
							if face_vertices_sorted == vertices:
								del face_vertices_dict[face]
								match_found = True
								break

						if not match_found:
							block_remaining_faces[Face(face_vertices)] = face_vertices_sorted

				face_vertices_dict.update(block_remaining_faces)

			if not (default_boundary_tag in self.boundaries.keys()):
				self.boundaries[default_boundary_tag] = _Boundary(default_boundary_tag)

			self.boundaries[default_boundary_tag].faces |= face_vertices_dict.keys()

		self._assign_vertexid()
		return f'''
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  8.0.0                                 |
|   \\\\  /    A nd           | Web:      www.OpenFOAM{self.of_distribution}                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
// This blockMeshDict was generated using blockmeshbuilder version {__version__},
// a Python module available at https://github.com/NauticalMile64/blockmeshbuilder
FoamFile
{{
	version     2.0;
	format      ascii;
	class       dictionary;
	object      blockMeshDict;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

convertToMeters {self.convert_to_meters};

{_format_section('geometry', self.geometries)}

{_format_section('vertices', self.valid_vertices)}

{_format_section('edges', self.edges)}

{_format_section('blocks', self.blocks)}

{_format_section('faces', self.faces)}

{_format_section('boundary', list(self.boundaries.values()))}

{_format_section('mergePatchPairs', self.mergePatchPairs)}

// ************************************************************************* //
'''

	def write_file(self, of_case_path=Path(), file_name='blockMeshDict', run_blockMesh=False, run_renumberMesh=True, **kwargs):
		of_case_path = Path(of_case_path)
		local_bmd_path = Path('system') / file_name

		with open(of_case_path / local_bmd_path, 'w') as infile:
			infile.write(self.format(**kwargs))

		if run_blockMesh:
			try:
				subprocess.run(["blockMesh", "-case", of_case_path, "-dict", local_bmd_path])
			# Note that in the cases when the blockMeshDict file cannot be found, blockMesh should throw it's own error
			except FileNotFoundError:
				warnings.warn("The system couldn't find the blockMesh application.")

			if run_renumberMesh:
				try:
					subprocess.run(["renumberMesh", "-case", of_case_path, "-overwrite"])
				except FileNotFoundError:
					warnings.warn("The system couldn't find the renumberMesh application.")
