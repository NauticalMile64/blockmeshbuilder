from io import StringIO
import warnings
import subprocess
from pathlib import Path
from .geometry import _of_geometries
from .tags import _Boundary


def _format_section(name, section_items):
	buf = StringIO()
	if name == 'geometry':
		brackets = ['{', '}']
	else:
		brackets = ['(', ')']

	buf.write(f'{name}\n')
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

	def __init__(self, metric='m', of_dist='.org', block_structure_only=False):
		if of_dist not in _of_geometries:
			warnings.warn(
				f'Unknown OpenFOAM distribution {of_dist}. The available options are {_of_geometries.keys()}. '
				f'\nSwitching to .org distribution.')
			of_dist = '.org'
		self.of_dist = of_dist
		self.of_available_geometries = _of_geometries[of_dist]
		self.convert_to_meters = self.metric_conversion_dict[metric]
		self.blocks = set()
		self.edges = set()
		self.boundaries = {}
		self.geometries = set()
		self.faces = set()
		self.block_structure_only = block_structure_only

	def add_hexblock(self, block):
		self.blocks.add(block)

	def add_edge(self, edge):
		self.edges.add(edge)

	def add_boundary_face(self, boundary_tag, face):
		if boundary_tag not in self.boundaries:
			self.boundaries[boundary_tag] = _Boundary(boundary_tag)

		self.boundaries[boundary_tag].add_face(face)

	def add_geometries(self, other_geometries):
		for geometry in other_geometries:
			if type(geometry) not in self.of_available_geometries:
				raise TypeError(
					f'Geometry of type {type(geometry)} is not implemented in the {self.of_dist} distribution.')

		self.geometries.update(other_geometries)

	def add_face(self, face):
		self.faces.add(face)

	def _assign_vertexid(self):
		valid_vertices = []

		i = 0
		for b in self.blocks:
			for v in b.vertices:
				if v not in valid_vertices:
					valid_vertices.append(v)
					v.index = i
					i += 1

		self.valid_vertices = valid_vertices

	def format(self, block_structure_only=False):
		if block_structure_only or self.block_structure_only:
			for block in self.blocks:
				block.cells = (1, 1, 1)

		self._assign_vertexid()
		return f'''
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

mergePatchPairs
(
);

// ************************************************************************* //
'''

	def write_file(self, of_case_path=Path(), file_name='blockMeshDict', block_structure_only=False,
				   run_blockMesh=False, density_scale=1.0):
		of_case_path = Path(of_case_path)
		local_bmd_path = Path('system') / file_name

		if density_scale != 1.0 and not block_structure_only:
			for block in self.blocks:
				block.cells = tuple(round(block.cells[i] * density_scale) for i in range(3))

		with open(of_case_path / local_bmd_path, 'w') as infile:
			infile.write(self.format(block_structure_only))

		if run_blockMesh:
			try:
				subprocess.run(["blockMesh", "-case", of_case_path, "-dict", local_bmd_path])
			except FileNotFoundError as error:
				warnings.warn("The system couldn't find the blockMesh application.")


'''
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  5.0.0                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
'''
