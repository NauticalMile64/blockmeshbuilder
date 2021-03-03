# -*- coding: future_fstrings -*-
# for compatibility with Py2.7
from __future__ import print_function

from sys import version_info

if version_info[0] == 2:
	from StringIO import StringIO
else:
	from io import StringIO

from six import string_types
import numpy as np
from numpy import sin, cos
import warnings


def cart_to_cart(crds):
	return crds.copy()


def cyl_to_cart(crds):
	ncrds = crds.copy()
	ncrds[0], ncrds[1] = crds[0] * cos(crds[1]), crds[0] * sin(crds[1])
	return ncrds


def cart_to_cyl(crds):
	ncrds = crds.copy()
	ncrds[0] = np.sqrt(crds[0] ** 2 + crds[1] ** 2)
	ncrds[1] = np.arctan2(crds[1], crds[0])
	return ncrds


class Point(object):
	def __init__(self, crds, conv_func=cart_to_cart):
		self.crds = np.asarray(crds, dtype=np.float64)
		self.conv_func = conv_func

	def format(self):
		ccrds = self.conv_func(self.crds)
		return f'( {ccrds[0]:18.15g} {ccrds[1]:18.15g} {ccrds[2]:18.15g} )'

	@staticmethod
	def _get_oth_coordinates(rhs):
		return rhs.crds if issubclass(type(rhs), Point) else rhs

	def __add__(self, rhs):
		return type(self)(self.crds + Point._get_oth_coordinates(rhs), self.conv_func)

	def __radd__(self, lhs):
		return self + lhs

	def __sub__(self, rhs):
		return type(self)(self.crds - Point._get_oth_coordinates(rhs), self.conv_func)

	def __rsub__(self, lhs):
		return type(self)(Point._get_oth_coordinates(lhs) - self.crds, self.conv_func)

	def __mul__(self, rhs):
		return type(self)(self.crds * Point._get_oth_coordinates(rhs), self.conv_func)

	def __rmul__(self, lhs):
		return self * lhs

	def __pow__(self, rhs):
		return type(self)(self.crds ** rhs, self.conv_func)

	def __truediv__(self, rhs):
		return type(self)(self.crds / Point._get_oth_coordinates(rhs), self.conv_func)

	def __isub__(self, rhs):
		self.crds -= Point._get_oth_coordinates(rhs)
		return self

	def __iadd__(self, rhs):
		self.crds += Point._get_oth_coordinates(rhs)
		return self

	def __imul__(self, rhs):
		self.crds *= Point._get_oth_coordinates(rhs)
		return self

	def __itruediv__(self, rhs):
		self.crds /= Point._get_oth_coordinates(rhs)
		return self

	def __lt__(self, rhs):
		return self.crds < self._get_oth_coordinates(rhs)

	def __gt__(self, lhs):
		return lhs < self

	def __neg__(self):
		return type(self)(-self.crds, self.conv_func)

	def __abs__(self):
		return abs(self.crds)

	def __getitem__(self, key):
		return self.crds[key]

	def __setitem__(self, key, value):
		self.crds[key] = value


class Projectable(object):
	def __init__(self, geometries=None):
		self.proj_g = set(geometries) if geometries else set()

	def proj_geom(self, geometry):
		self.proj_g.add(geometry)

	def format_geom(self):
		if self.proj_g:
			geom_names = ' '.join(geom.name for geom in self.proj_g)
			return 'project ', f'({geom_names}) '
		else:
			return '', ''


class Vertex(Point, Projectable):
	def __init__(self, crds, conv_func=cart_to_cart, name='', geometries=None):
		Point.__init__(self, crds, conv_func)
		Projectable.__init__(self, geometries)
		self.name = name
		self.index = None

	def format(self):
		com = f'{self.index} {self.name}'
		vertex_str = Point.format(self)
		proj_str, geom_name = self.format_geom()
		return f'{proj_str}{vertex_str} {geom_name}// {com:s}'


def _dict_format(dict_name, data_dict, indent_level=1):
	buf = StringIO()
	buf.write(dict_name)
	indent = '\t' * indent_level
	buf.write(f'\n{indent}{{')

	for label, data in data_dict.items():
		buf.write(f'\n{indent}\t{label} {data};')

	buf.write(f'\n{indent}}}')
	return buf.getvalue()


class Geometry(object):
	_unique_id = 0

	def __init__(self, name=''):
		self.name = f'{name}_id-{self._unique_id}'
		Geometry._unique_id += 1

	def do_format(self, data_dict):
		return _dict_format(self.name, data_dict)


class Plane(Geometry):
	plane_type = 'Plane Base Class'
	plane_data = {}

	def __init__(self, name):
		Geometry.__init__(self, name)

	def do_format(self, plane_data):
		return Geometry.do_format(self,
								  {
									  'type': 'searchablePlane',
									  'planeType': self.plane_type,
									  _dict_format(f'{self.plane_type}Dict', plane_data, indent_level=2): ''
								  })


class PlanePointAndNormal(Plane):
	plane_type = 'pointAndNormal'

	def __init__(self, point, normal, name):
		Plane.__init__(self, name)
		self.point = point
		self.normal = normal

	def format(self):
		return Plane.do_format(self,
							   {
								   'point': self.point.format(),
								   'normal': self.normal.format()
							   })


class PlaneEmbeddedPoints(Plane):
	plane_type = 'embeddedPoints'

	def __init__(self, point1, point2, point3, name):
		Plane.__init__(self, name)
		self.points = (point1, point2, point3)

	def format(self):
		return Plane.do_format(self,
							   {
								   'point1': self.points[0].format(),
								   'point2': self.points[1].format(),
								   'point3': self.points[2].format()
							   })


class PlaneEquation(Plane):
	plane_type = 'planeEquation'

	def __init__(self, a, b, c, d, name):
		Plane.__init__(self, name)
		self.equation_coeffs = (a, b, c, d)

	def format(self):
		return Plane.do_format(self,
							   {
								   'a': f'{self.equation_coeffs[0]:18.15g}',
								   'b': f'{self.equation_coeffs[1]:18.15g}',
								   'c': f'{self.equation_coeffs[2]:18.15g}',
								   'd': f'{self.equation_coeffs[3]:18.15g}',
							   })


class Sphere(Geometry):
	def __init__(self, center, radius, name):
		Geometry.__init__(self, name)
		self.center = center
		self.radius = radius

	def format(self):
		return Geometry.do_format(self,
								  {
									  'type': 'searchableSphere',
									  'centre': self.center.format(),
									  'radius': f'{self.radius:18.15g}'
								  })


class Cylinder(Geometry):
	def __init__(self, point1, point2, radius, name):
		Geometry.__init__(self, name)
		self.point1 = point1
		self.point2 = point2
		self.radius = radius

	def format(self):
		return Geometry.do_format(self,
								  {
									  'type': 'searchableCylinder',
									  'point1': self.point1.format(),
									  'point2': self.point2.format(),
									  'radius': f'{self.radius:18.15g}'
								  })


class Cone(Geometry):
	small_value = np.finfo(np.float64).eps

	def __init__(self, point1, point2, radius1, radius2, name, inner_radius1=small_value, inner_radius2=small_value):
		Geometry.__init__(self, name)
		self.point1 = point1
		self.point2 = point2
		self.radius1 = radius1
		self.radius2 = radius2
		self.inner_radius1 = inner_radius1
		self.inner_radius2 = inner_radius2

	def format(self):
		return Geometry.do_format(self,
								  {
									  'type': 'searchableCone',
									  'point1': self.point1.format(),
									  'point2': self.point2.format(),
									  'radius1': f'{self.radius1:18.15g}',
									  'radius2': f'{self.radius2:18.15g}',
									  'inner_radius1': f'{self.inner_radius1:18.15g}',
									  'inner_radius2': f'{self.inner_radius2:18.15g}'
								  })


class Face(object):
	def __init__(self, vertices, name=''):

		# vertices is a 2x2 array
		self.vertices = vertices
		self.name = name
		self.proj_g = None

	def proj_geom(self, geometry):
		if self.proj_g:
			warnings.warn(
				f'Face-{self.name} has already been projected to {self.proj_g.name}; over-writing with {geometry.name}.')

		self.proj_g = geometry

	def format(self, write_proj=True):

		v_arr = self.vertices
		vts = [v_arr[0][0], v_arr[0][1], v_arr[1][0], v_arr[1][1]]

		index = ' '.join(str(v.index) for v in vts)
		com = ' '.join(v.name for v in vts)

		proj_str, geom_name = '', ''
		if write_proj:
			proj_str = 'project '
			geom_name = self.proj_g.name

		res_str = f'{proj_str}({index:s}) {geom_name}// {self.name:s} ({com:s})'

		# If either vertex has not been included in any blocks the edge is meaningless
		if 'None' in index:
			res_str = '// ' + res_str

		return res_str


class GradingElement(object):
	pass


class SimpleGradingElement(GradingElement):
	def __init__(self, d):
		self.d = d

	def format(self):
		return str(self.d)


class MultiGradingElement(GradingElement):
	def __init__(self, len_pcts, cell_pcts, exp_ratios):
		self.len_pcts = len_pcts
		self.cell_pcts = cell_pcts
		self.exp_ratios = exp_ratios

	def format(self):
		return '({})'.format(' '.join(f'({lp:.4f} {nc:.4f} {ex:.4f})'
									  for lp, nc, ex in zip(self.len_pcts, self.cell_pcts, self.exp_ratios)))


class Grading(object):
	def __init__(self, grading_elements):
		self.grading_elements = grading_elements

	def format(self):
		return '{{0}}Grading ({0})'.format(' '.join(ge.format() for ge in self.grading_elements))


class SimpleGrading(Grading):
	def __init__(self, grading_elements):
		assert (len(grading_elements) == 3)
		Grading.__init__(self, grading_elements)

	def format(self):
		return Grading.format(self).format('simple')


class EdgeGrading(Grading):
	def __init__(self, grading_elements):
		assert (len(grading_elements) == 12)
		Grading.__init__(self, grading_elements)

	def format(self):
		return Grading.format(self).format('edge')


# Helper function
def get_grading_info(len_pcts, dens):
	cum_lens = np.insert(np.cumsum(len_pcts), 0, 0.0)
	adens = np.array([np.trapz(dens[s:s + 2], cum_lens[s:s + 2]) for s in range(len_pcts.size)])
	exp_ratios = np.divide(dens[:-1], dens[1:])
	return len_pcts, adens, exp_ratios


uniformGradingElement = SimpleGradingElement(1)
uniformGrading = SimpleGrading([uniformGradingElement] * 3)


class ZoneTag(object):
	def __init__(self, name):
		self.name = name


DEFAULT_ZONE_TAG = ZoneTag('DEFAULT')


class HexBlock(object):
	def __init__(self, vertices, cells, zone_tag=DEFAULT_ZONE_TAG, grading=uniformGrading):
		self.vertices = vertices
		self.cells = cells
		self.zone_tag = zone_tag
		self.grading = grading

	def format(self):
		index = ' '.join(str(v.index) for v in self.vertices)
		vcom = ' '.join(v.name for v in self.vertices)
		cls = self.cells
		zone_name = self.zone_tag.name
		return f'hex ({index:s}) {zone_name:s} ({cls[0]:d} {cls[1]:d} {cls[2]:d}) {self.grading.format():s}  // {zone_name:s} ({vcom:s})'

	def get_face(self, index, name=None):
		"""Generate Face object
		index is number or keyword to identify the face of Hex
			0 = 'w' = 'xm' = '-100' = (0 4 7 3)
			1 = 'e' = 'xp' = '100' = (1 2 5 6)
			2 = 's' = 'ym' = '0-10' = (0 1 5 4)
			3 = 'n' = 'yp' = '010' = (2 3 7 6)
			4 = 'b' = 'zm' = '00-1' = (0 3 2 1)
			5 = 't' = zp' = '001' = (4 5 6 7)
		name is given to Face instance. If omitted, name is automatically
			genaratied like ('f-' + self.name + '-w')
		"""
		kw_to_index = {
			'w': 0, 'xm': 0, '-100': 0,
			'e': 1, 'xp': 1, '100': 1,
			's': 2, 'ym': 2, '0-10': 2,
			'n': 3, 'yp': 3, '010': 3,
			'b': 4, 'zm': 4, '00-1': 4,
			't': 5, 'zp': 5, '001': 5}
		index_to_vertex = [
			((0, 4), (7, 3)),
			((1, 2), (6, 5)),
			((0, 1), (5, 4)),
			((2, 3), (7, 6)),
			((0, 3), (2, 1)),
			((4, 5), (6, 7))]
		index_to_defaultsuffix = [
			'f-{}-w',
			'f-{}-n',
			'f-{}-s',
			'f-{}-n',
			'f-{}-b',
			'f-{}-t']

		if isinstance(index, string_types):
			index = kw_to_index[index]

		face_vertices = [[self.vertices[i] for i in row] for row in index_to_vertex[index]]

		if name is None:
			name = index_to_defaultsuffix[index].format(self.name)

		return Face(face_vertices, name)


class Edge(object):
	def __init__(self, vertices, name=''):
		# http://www.openfoam.org/docs/user/mesh-description.php

		self.vertices = vertices
		self.name = name

	@staticmethod
	def is_relevant():
		return True

	def format(self):
		indices = [v.index for v in self.vertices]
		index = ' '.join(str(ind) for ind in indices)
		vcom = ' '.join(str(v.name) for v in self.vertices)
		res_str = f'{{0}} {index:s} {{1}} // {self.name:s} ({vcom:s})'

		# If either vertex has not been included in any blocks the edge is meaningless
		if None in indices:
			res_str = '// ' + res_str

		return res_str


class ArcEdge(Edge):
	def __init__(self, vertices, arc_mid_point, name=''):
		Edge.__init__(self, vertices, name)
		self.arc_mid_point = arc_mid_point

	def format(self):
		return Edge.format(self).format('arc', self.arc_mid_point.format())


class CurvedEdge(Edge):
	edge_type = 'Curved Edge Base Class'

	def __init__(self, vertices, points, name=''):
		Edge.__init__(self, vertices, name)
		self.points = points

	def format(self):
		buf = StringIO()

		buf.write(Edge.format(self).format(self.edge_type, ''))
		buf.write('\n     (\n')
		for p in self.points:
			buf.write(f'         {p.format()}\n')
		buf.write('\n     )\n')
		buf.write('')
		return buf.getvalue()


class SplineCurvedEdge(CurvedEdge):
	edge_type = 'spline'


class PolyLineCurvedEdge(CurvedEdge):
	edge_type = 'polyLine'


class BSplineCurvedEdge(CurvedEdge):
	edge_type = 'BSpline'


class ProjectionEdge(Edge, Projectable):
	def __init__(self, vertices, name='', geometries=[]):
		Edge.__init__(self, vertices, name)
		Projectable.__init__(self, geometries)

	def is_relevant(self):
		return bool(self.proj_g)

	def format(self):
		return Edge.format(self).format(*self.format_geom())


class BoundaryTag(object):
	def __init__(self, name, type_='patch'):
		self.type_ = type_
		self.name = name


class _Boundary(object):
	def __init__(self, boundary_tag):
		assert (isinstance(boundary_tag, BoundaryTag))
		self.boundary_tag = boundary_tag
		self.faces = set()

	def add_face(self, face):
		self.faces.add(face)

	def format(self):
		buf = StringIO()

		buf.write(self.boundary_tag.name + '\n')
		buf.write('{\n')
		buf.write(f'    type {self.boundary_tag.type_};\n')
		buf.write('    faces\n')
		buf.write('    (\n')
		for f in self.faces:
			buf.write(f'        {f.format(False)}\n')
		buf.write('    );\n')
		buf.write('}')
		return buf.getvalue()


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


class BlockMeshDict(object):
	metricsym_to_conversion = {
		'km': 1000,
		'm': 1,
		'cm': 0.01,
		'mm': 0.001,
		'um': 1e-6,
		'nm': 1e-9,
		'A': 1e-10,
		'Angstrom': 1e-10}

	_of_geometries = {
		'.org': {
			PlanePointAndNormal, PlaneEmbeddedPoints, PlaneEquation, Cylinder, Sphere,
			# closedTriSurfaceMesh, Box, Disc, ExtrudedCircle, Plate, SurfaceCollection, SurfaceWithGaps, triSurfaceMesh
		},
		'.com': {
			PlanePointAndNormal, PlaneEmbeddedPoints, PlaneEquation, Cylinder, Sphere, Cone
			# distributedTriSurfaceMesh, Box, Disc, ExtrudedCircle, Plate, SurfaceCollection, SurfaceWithGaps, triSurfaceMesh
		}
	}

	def __init__(self, metric='m', of_dist='.org', block_structure_only=False):
		if of_dist not in self._of_geometries:
			warnings.warn(
				f'Unknown OpenFOAM distribution {of_dist}. The available options are {self._of_geometries.keys()}. Switching to .org distribution.')
			of_dist = '.org'
		self.of_dist = of_dist
		self.of_available_geometries = self._of_geometries[of_dist]
		self.convert_to_meters = self.metricsym_to_conversion[metric]
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

	def write_file(self, of_case_path, block_structure_only=False):
		if of_case_path[-1] in r"/\\":  # Test whether a directory or the full file path is specified
			of_case_path += 'blockMeshDict'
		with open(of_case_path, 'w') as infile:
			infile.write(self.format(block_structure_only))


'''
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  5.0.0                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
'''
