# -*- coding: future_fstrings -*-
# for compatibility with Py2.7
from __future__ import unicode_literals, print_function
from six import string_types

import io

class Point(object):
	def __init__(self, x, y, z):
		self.x = x
		self.y = y
		self.z = z
	
	def format(self):
		return f'( {self.x:18.15g} {self.y:18.15g} {self.z:18.15g} )'
	
	def __lt__(self, rhs):
		return (self.z, self.y, self.x) < (rhs.z, rhs.y, rhs.z)
	
	def __add__(self, rhs):
		return Point(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
	
	def __sub__(self, rhs):
		return Point(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)


class Vertex(Point):
	def __init__(self, x, y, z, name, index=None):
		
		Point.__init__(self, x, y, z)
		
		self.name = name
		self.alias = set([name])  # aliasname, self.name should be included
		
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
			geom_name = f'({self.proj_g.name})'
		
		return f'{proj_str} {vertex_str} {geom_name} // {com:s}'
	
	def proj_geom(self, geometry):
		self.proj_g = geometry


class Geometry(object):
	def __init__(self, name):
		assert(name)
		self.name = name
	
	def format(self,type,data_dict):
		buf = io.StringIO()
		
		buf.write(self.name)
		buf.write('\n\t{')
		buf.write(f'\n\t\ttype {type};')
		for label, data in data_dict.items():
			buf.write(f'\n\t\t{label} {data};')
		
		buf.write('\n\t}')
		return buf.getvalue()


class Sphere(Geometry):
	def __init__(self, name, center, radius):
		Geometry.__init__(self, name)
		self.center = center
		self.radius = radius
	
	def format(self):
		return Geometry.format(self,'searchableSphere',
			{
				'centre' : self.center.format(),
				'radius' : f'{self.radius:18.15g}'
			})


class Cylinder(Geometry):
	def __init__(self, name, point1, point2, radius):
		Geometry.__init__(self, name)
		self.point1 = point1
		self.point2 = point2
		self.radius = radius
	
	def format(self):
		return Geometry.format(self,'searchableCylinder',
			{
				'point1' : self.point1.format(),
				'point2' : self.point2.format(),
				'radius' : f'{self.radius:18.15g}'
			})


class Face(object):
	def __init__(self, vertices, name):
		self.vertices = vertices
		self.name = name
	
	def format(self, prj_geom=''):
		index = ' '.join(str(v.index) for v in self.vertices)
		com = ' '.join(v.name for v in self.vertices)
		return f'({index:s}) {prj_geom} // {self.name:s} ({com:s})'


class GradingElement(object):
	pass


class SimpleGradingElement(GradingElement):
	def __init__(self, d):
		self.d = d
	
	def format(self):
		return str(self.d)


class MultiGradingElement(GradingElement):
	def __init__(self,len_pcts,cell_pcts,exp_ratios):
		self.len_pcts = len_pcts
		self.cell_pcts = cell_pcts
		self.exp_ratios = exp_ratios
	
	def format(self):
		return '({})'.format(' '.join(f'({lp} {nc} {ex})' 
				for lp,nc,ex in zip(self.len_pcts, self.cell_pcts, self.exp_ratios)))


class Grading(object):
	def __init__(self,gradingElements):
		self.gradingElements = gradingElements
	
	def format(self):
		return '{{0}}Grading ({0})'.format(' '.join(ge.format() for ge in self.gradingElements))


class SimpleGrading(Grading):
	def __init__(self, gradingElements):
		assert(len(gradingElements) == 3)
		Grading.__init__(self, gradingElements)
	
	def format(self):
		return Grading.format(self).format('simple')


class EdgeGrading(Grading):
	def __init__(self, gradingElements):
		assert(len(gradingElements) == 12)
		Grading.__init__(self, gradingElements)
	
	def format(self):
		return Grading.format(self).format('edge')


uniformGradingElement = SimpleGradingElement(1)
uniformGrading = SimpleGrading([uniformGradingElement]*3)
class HexBlock(object):
	def __init__(self, vertices, cells, name, grading=uniformGrading):
		"""Initialize HexBlock instance
		vnames is the vertex names in order described in
			http://www.openfoam.org/docs/user/mesh-description.php
		cells is number of cells devied into in each direction
		name is the unique name of the block
		grading is grading method.
		"""
		self.vertices = vertices
		self.cells = cells
		self.name = name
		self.grading = grading
	
	def format(self):
		index = ' '.join(str(v.index) for v in self.vertices)
		vcom = ' '.join(v.name for v in self.vertices)
		cls = self.cells
		return f'hex ({index:s}) {self.name:s} ({cls[0]:d} {cls[1]:d} {cls[2]:d}) {self.grading.format():s}  // {self.name:s} ({vcom:s})'
	
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
		
		indices = [v.index for v in self.vertices]
		index = ' '.join(str(ind) for ind in indices)
		vcom = ' '.join(str(v.name) for v in self.vertices)
		res_str = f'{{0}} {index:s} {{1}} // {self.name:s} ({vcom:s})'
		
		# If either vertex has not been included in any blocks the edge is meaningless
		if None in indices:
			res_str = '// ' + res_str
		
		return res_str

class ArcEdge(Edge):
	def __init__(self, vertices, name, interVertex):
		
		Edge.__init__(self, vertices, name)
		self.interVertex = interVertex
	
	def format(self):
		iv = self.interVertex
		return Edge.format(self).format('arc',
				f'({iv.x:f} {iv.y:f} {iv.z:f}) ')


class SplineEdge(Edge):
	def __init__(self, vertices, name, points):
		"""
		  http://www.openfoam.org/docs/user/mesh-description.php
		"""
		
		Edge.__init__(self,vertices,name)
		self.points = points
	
	def format(self):
		buf = io.StringIO()
		
		buf.write(Edge.format(self).format('spline',''))
		buf.write('\n     (\n')
		for p in self.points:
			buf.write(f'         {p.format()}\n')
		buf.write('\n     )\n')
		buf.write('')
		return buf.getvalue()


class ProjectionEdge(Edge):
	def __init__(self, vertices, name, geom):
		Edge.__init__(self, vertices, name)
		self.proj_geom = geom
	
	def format(self):
		return Edge.format(self).format('project',
				f'({self.proj_geom.name})')


class Boundary(object):
	def __init__(self, type_, name, faces=[]):
		self.type_ = type_
		self.name = name
		self.faces = faces
	
	def add_face(self, face):
		self.faces.append(face)
	
	def format(self):
		buf = io.StringIO()

		buf.write(self.name + '\n')
		buf.write('{\n')
		buf.write(f'    type {self.type_};\n')
		buf.write('    faces\n')
		buf.write('    (\n')
		for f in self.faces:
			s = f.format()
			buf.write('        {}\n'.format(s))
		buf.write('    );\n')
		buf.write('}')
		return buf.getvalue()

def _format_section(name, secList):
	buf = io.StringIO()
	buf.write(f'{name}\n')
	buf.write('(\n')
	for item in secList:
		buf.write(f'    {item.format()}\n')
	buf.write(');')
	return buf.getvalue()

class BlockMeshDict(object):
	def __init__(self):
		self.convert_to_meters = 1.0
		self.vertices = {}
		self.blocks = {}
		self.edges = {}
		self.boundaries = {}
		self.geometries = {}
		self.proj_faces = {}
	
	def set_metric(self, metric):
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
	
	def add_boundary(self, boundary):
		self.boundaries[boundary.name] = boundary
	
	def add_geometry(self, name, geometry):
		self.geometries[name] = geometry
	
	def add_proj_face(self, name, face, proj_geometry_name):
		self.proj_faces[name] = {'face' : face, 'proj_geom' : proj_geometry_name}
	
	def assign_vertexid(self):
		valid_vertices = []
		i = 0
		
		for b in self.blocks.values():
			for v in b.vertices:
				if v not in valid_vertices:
					valid_vertices.append(v)
					v.index = i
					i += 1

		self.valid_vertices = valid_vertices
	
	def format_faces_section(self):
		buf = io.StringIO()
		buf.write('faces\n')
		buf.write('(\n')
		for pFace in self.proj_faces.values():
			buf.write(' project {0} \n'.format(
				pFace['face'].format(pFace['proj_geom'].name)))
		
		buf.write(');')
		return buf.getvalue()
	
	def format(self):
		return f'''/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  5.0.0                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{{
	version     2.0;
	format      ascii;
	class       dictionary;
	object      blockMeshDict;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

convertToMeters {self.convert_to_meters};

{_format_section('geometry',list(self.geometries.values()))}

{_format_section('vertices',self.valid_vertices)}

{_format_section('edges',list(self.edges.values()))}

{_format_section('blocks',list(self.blocks.values()))}

{self.format_faces_section()}

{_format_section('blocks',list(self.blocks.values()))}

mergePatchPairs
(
);

// ************************************************************************* //
'''