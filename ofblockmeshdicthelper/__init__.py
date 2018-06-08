# for compatibility to Py2.7
from __future__ import unicode_literals, print_function
from six import string_types

import io
from collections import Iterable
from string import Template


class Point(object):
	def __init__(self, x, y, z):
		self.x = x
		self.y = y
		self.z = z

	def format(self):
		return '( {0:18.15g} {1:18.15g} {2:18.15g} )'.format(
			self.x, self.y, self.z)

	def __lt__(self, rhs):
		return (self.z, self.y, self.x) < (rhs.z, rhs.y, rhs.z)
		
	def __add__(self, rhs):
		return Point(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
	
	def __sub__(self, rhs):
		return Point(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)

class Vertex(Point):
	def __init__(self, x, y, z, name, index=None):
		
		Point.__init__(self, x, y, z)
		
		self.name = name  # identical name
		self.alias = set([name])  # aliasname, self.name should be included

		# seqential index which is assigned at final output
		# for blocks, edges, boundaries
		self.index = None
		self.proj_g = None
	
	def format(self):
		com = str(self.index) + ' ' + self.name
		if len(self.alias) > 1:
			com += ' : '
			com += ' '.join(self.alias)
		
		vertex_str = Point.format(self)
		
		proj_str, geom_name = '',''
		if self.proj_g:
			proj_str = 'project'
			geom_name = '({})'.format(self.proj_g.name)
		
		return '{0} {1} {2} // {3:s}'.format(
			proj_str, vertex_str, geom_name, com)
	
	def proj_geom(self, geometry):
		self.proj_g = geometry

class Geometry(object):
	def __init__(self, name):
		self.name = name

class Sphere(Geometry):
	def __init__(self, name, center, radius):
		Geometry.__init__(self, name)
		self.center = center
		self.radius = radius
	
	def format(self):
		return '''{0}
	{{
		type searchableSphere;
		centre {1};
		radius {2:18.15g};
	}}
'''.format(self.name, self.center.format(), self.radius)

class Face(object):
	def __init__(self, vnames, name):
		"""
		vname is list or tuple of vertex names
		"""
		self.vnames = vnames
		self.name = name

	def format(self, prj_geom=''):
		"""Format instance to dump
		vertices is dict of name to Vertex
		"""
		index = ' '.join(str(vertices[vn].index) for vn in self.vnames)
		com = ' '.join(self.vnames)  # for comment
		return '({0:s}) {3} // {1:s} ({2:s})'.format(index, self.name, com, prj_geom)


class Grading(object):
	"""base class for Simple- and Edge- Grading"""
	pass


class SimpleGradingElement(object):
	"""x, y or z Element of simpleGrading. adopted to multi-grading
	"""
	def __init__(self, d):
		"""initialization
		d is single number for expansion ratio
		  or iterative object consits (dirction ratio, cell ratio, expansion ratio)
		"""
		self.d = d

	def format(self):
		if isinstance(self.d, Iterable):
			s = io.StringIO()
			s.write('( ')
			for e in self.d:
				s.write('( {0:g} {1:g} {2:g} ) '.format(e[0], e[1], e[2]))
			s.write(')')
			return s.getvalue()
		else:
			return str(self.d)


class SimpleGrading(Grading):
	"""configutation for 'simpleGrading'
	multi-grading is not implemented yet
	"""
	def __init__(self, x, y, z):
		if not isinstance(x, SimpleGradingElement):
			self.x = SimpleGradingElement(x)
		else:
			self.x = x
		if not isinstance(y, SimpleGradingElement):
			self.y = SimpleGradingElement(y)
		else:
			self.y = y
		if not isinstance(z, SimpleGradingElement):
			self.z = SimpleGradingElement(z)
		else:
			self.z = z

	def format(self):
		return 'simpleGrading ({0:s} {1:s} {2:s})'.format(self.x.format(), self.y.format(), self.z.format())


class HexBlock(object):
	def __init__(self, vertices, cells, name, grading=SimpleGrading(1, 1, 1)):
		"""Initialize HexBlock instance
		vnames is the vertex names in order descrived in
			http://www.openfoam.org/docs/user/mesh-description.php
		cells is number of cells devied into in each direction
		name is the uniq name of the block
		grading is grading method.
		"""
		self.vertices = vertices
		self.cells = cells
		self.name = name
		self.grading = grading

	def format(self):
		"""Format instance to dump
		vertices is dict of name to Vertex
		"""
		index = ' '.join(str(v.index) for v in self.vertices)
		vcom = ' '.join(v.name for v in self.vertices)  # for comment
		return 'hex ({0:s}) {2:s} ({1[0]:d} {1[1]:d} {1[2]:d}) '\
			   '{4:s}  // {2:s} ({3:s})'.format(
					index, self.cells, self.name, vcom, self.grading.format())

	def face(self, index, name=None):
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
			(0, 4, 7, 3),
			(1, 2, 6, 5),
			(0, 1, 5, 4),
			(2, 3, 7, 6),
			(0, 3, 2, 1),
			(4, 5, 6, 7)]
		index_to_defaultsuffix = [
			'f-{}-w',
			'f-{}-n',
			'f-{}-s',
			'f-{}-n',
			'f-{}-b',
			'f-{}-t']

		if isinstance(index, string_types):
			index = kw_to_index[index]

		face_vertices = tuple([self.vertices[i] for i in index_to_vertex[index]])
		if name is None:
			name = index_to_defaultsuffix[index].format(self.name)
		return Face(face_vertices, name)


class Edge(object):
	def __init__(self, vertices, name):
	
		#http://www.openfoam.org/docs/user/mesh-description.php
		
		self.vertices = vertices
		self.name = name

	def format(self):
		index = ' '.join(str(v.index) for v in self.vertices)
		vcom = ' '.join(str(v.name) for v in self.vertices)  # for comment
		return '{{0}} {0:s} {{1}}'\
				'// {1:s} ({2:s})'.format(
						index, self.name, vcom)
		
class ArcEdge(Edge):
	def __init__(self, vertices, name, interVertex):
		
		Edge.__init__(self, vertices, name)
		self.interVertex = interVertex

	def format(self):
		
		return Edge.format(self).format('arc',
				'({0.x:f} {0.y:f} {0.z:f}) '.format(self.interVertex))

class SplineEdge(Edge):
	def __init__(self, vertices, name, points):
		"""
		  http://www.openfoam.org/docs/user/mesh-description.php
		"""
		
		Edge.__init__(self,vertices,name)
		self.points = points

	def format(self):
		"""Format instance to dump
		vertices is dict of name to Vertex
		"""
		
		buf = io.StringIO()

		buf.write(Edge.format(self).format('spline',''))
		buf.write('\n     (\n')
		for p in self.points:
			buf.write('         '+p.format()+'\n')
		buf.write('\n     )\n')
		buf.write('')
		return buf.getvalue()

class ProjectionEdge(Edge):
	def __init__(self, vertices, name, geom):
		Edge.__init__(self, vertices, name)
		self.proj_geom = geom

	def format(self):
		return Edge.format(self).format('project',
				'({})'.format(self.proj_geom.name))

class Boundary(object):
	def __init__(self, type_, name, faces=[]):
		""" initialize boundary
		type_ is type keyword (wall, patch, empty, ..)
		name is name of boundary emelment
		faces is faces which are applied with this boundary conditions
		"""
		self.type_ = type_
		self.name = name
		self.faces = faces

	def add_face(self, face):
		"""add face instance
		face is a Face instance (not name) to be added
		"""
		self.faces.append(face)

	def format(self):
		"""Format instance to dump
		vertices is dict of name to Vertex
		"""
		buf = io.StringIO()

		buf.write(self.name + '\n')
		buf.write('{\n')
		buf.write('    type {};\n'.format(self.type_))
		buf.write('    faces\n')
		buf.write('    (\n')
		for f in self.faces:
			s = f.format()
			buf.write('        {}\n'.format(s))
		buf.write('    );\n')
		buf.write('}')
		return buf.getvalue()


class BlockMeshDict(object):
	def __init__(self):
		self.convert_to_meters = 1.0
		self.vertices = {}  # mapping of uniq name to Vertex object
		self.blocks = {}
		self.edges = {}
		self.boundaries = {}
		self.geometries = {}
		self.proj_faces = {}

	def set_metric(self, metric):
		"""set self.comvert_to_meters by word"""
		metricsym_to_conversion = {
			'km': 1000,
			'm': 1,
			'cm': 0.01,
			'mm': 0.001,
			'um': 1e-6,
			'nm': 1e-9,
			'A': 1e-10,
			'Angstrom': 1e-10}
		self.convert_to_meters = metricsym_to_conversion[metric]

	def reduce_vertex(self, name1, *names):
		"""treat name1, name2, ... as same point.

		name2.alias, name3.alias, ... are merged with name1.alias
		the key name2, name3, ... in self.vertices are kept and mapped to
		same Vertex instance as name1
		"""
		v = self.vertices[name1]
		for n in names:
			w = self.vertices[n]
			v.alias.update(w.alias)
			# replace mapping from n w by to v
			self.vertices[n] = v

	def add_hexblock(self, block, name):
		self.blocks[name] = block
	
	def add_edge(self, edge, name):
		self.edges[name] = edge

	def add_boundary(self, boundary, name):
		self.boundaries[name] = boundary
	
	def add_geometry(self, name, geometry):
		self.geometries[name] = geometry
		return geometry
	
	def add_proj_face(self, name, face, proj_geometry_name):
		self.proj_faces[name] = {'face' : face, 'proj_geom' : proj_geometry_name}
		return face
	
	def assign_vertexid(self):
		"""1. create list of Vertex which are referred by blocks only.
		2. sort vertex according to (x, y, z)
		3. assign sequence number for each Vertex
		4. sorted list is saved as self.valid_vertices
		"""

		# gather 'uniq' names which are refferred by blocks
		valid_vertices = []
		i = 0
		
		for b in self.blocks.values():
			for v in b.vertices:
				if v.name == 'cv-3':
					print('why')
					print(i)
					for vr in valid_vertices:
						if v == vr:
							print(v)
							print(vr)
							print(vr.name)
					print(v in valid_vertices)
					
				if v not in valid_vertices:
					valid_vertices.append(v)
					v.index = i
					i += 1

		self.valid_vertices = valid_vertices
		#for i, v in enumerate(self.valid_vertices):
		#    v.index = i

	def format_vertices_section(self):
		"""format vertices section.
		assign_vertexid() should be called before this method, because
		self.valid_vertices should be available and member self.valid_vertices
		should have valid index.
		"""
		buf = io.StringIO()
		buf.write('vertices\n')
		buf.write('(\n')
		for v in self.valid_vertices:
			buf.write('    ' + v.format() + '\n')
		buf.write(');')
		return buf.getvalue()

	def format_geometry_section(self):
		"""format geometry section.
		"""
		buf = io.StringIO()
		buf.write('geometry\n')
		buf.write('{\n')
		for g in self.geometries.values():
			buf.write('    ' + g.format() + '\n')
		buf.write('};')
		return buf.getvalue()
	
	def format_blocks_section(self):
		"""format blocks section.
		assign_vertexid() should be called before this method, because
		vertices referred by blocks should have valid index.
		"""
		buf = io.StringIO()
		buf.write('blocks\n')
		buf.write('(\n')
		for b in self.blocks.values():
			buf.write('    ' + b.format() + '\n')
		buf.write(');')
		return buf.getvalue()

	def format_edges_section(self):
		"""format edges section.
		assign_vertexid() should be called before this method, because
		vertices referred by blocks should have valid index.
		"""
		buf = io.StringIO()
		buf.write('edges\n')
		buf.write('(\n')
		for e in self.edges.values():
			buf.write('  ' + e.format() + '\n')
		buf.write(');')
		return buf.getvalue()
	
	def format_faces_section(self):
		"""format faces section.
		assign_vertexid() should be called before this method, because
		vertices refered by blocks should have valid index.
		"""
		buf = io.StringIO()
		buf.write('faces\n')
		buf.write('(\n')
		for pFace in self.proj_faces.values():
			buf.write(' project {0} \n'.format(
				pFace['face'].format(pFace['proj_geom'])))
		
		buf.write(');')
		return buf.getvalue()
	
	def format_boundary_section(self):
		"""format boundary section.
		assign_vertexid() should be called before this method, because
		vertices referred by faces should have valid index.
		"""
		buf = io.StringIO()
		buf.write('boundary\n')
		buf.write('(\n')
		for b in self.boundaries.values():
			# format Boundary instance and add indent
			indent = ' ' * 4
			s = b.format().replace('\n', '\n'+indent)
			buf.write(indent + s + '\n')
		buf.write(');')
		return buf.getvalue()

	def format_mergepatchpairs_section(self):
		return '''\
mergePatchPairs
(
);'''

	def format(self):
		template = Template(r'''/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  5.0.0                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
	version     2.0;
	format      ascii;
	class       dictionary;
	object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

convertToMeters $metricconvert;

$geometry

$vertices

$edges

$blocks

$faces

$boundary

$mergepatchpairs

// ************************************************************************* //
''')

		return template.substitute(
			metricconvert=str(self.convert_to_meters),
			geometry=self.format_geometry_section(),
			vertices=self.format_vertices_section(),
			edges=self.format_edges_section(),
			blocks=self.format_blocks_section(),
			faces=self.format_faces_section(),
			boundary=self.format_boundary_section(),
			mergepatchpairs=self.format_mergepatchpairs_section())