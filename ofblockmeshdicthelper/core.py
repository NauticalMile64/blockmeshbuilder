# -*- coding: future_fstrings -*-
# for compatibility with Py2.7
from __future__ import unicode_literals, print_function
from six import string_types
from math import sin,cos

import io

def cart_to_cart(crds):
	return crds.copy()

def cyl_to_cart(crds):
	ncrds = crds.copy()
	ncrds[0],ncrds[1] = crds[0]*cos(crds[1]),crds[0]*sin(crds[1])
	return ncrds

class Point(object):
	def __init__(self, crds, conv_func=cart_to_cart):
		self.crds = crds
		self.conv_func = conv_func
	
	def format(self):
		ccrds = self.conv_func(self.crds)
		return f'( {ccrds[0]:18.15g} {ccrds[1]:18.15g} {ccrds[2]:18.15g} )'
	
	def __lt__(self, rhs):
		return (self.z, self.y, self.x) < (rhs.z, rhs.y, rhs.z)
	
	def __neg__(self):
		return Point([-self.x, -self.y, -self.z], self.conv_func)
	
	def __add__(self, rhs):
		return Point([self.x + rhs.x, self.y + rhs.y, self.z + rhs.z], self.conv_func)
	
	def __sub__(self, rhs):
		return Point([self.x - rhs.x, self.y - rhs.y, self.z - rhs.z], self.conv_func)
	
	def __mul__(self, rhs):
		return Point([self.x*rhs, self.y*rhs, self.z*rhs], self.conv_func)
	
	def __rmul__(self, lhs):
		return self*lhs
	
	def __truediv__(self, rhs):
		return Point([self.x/rhs, self.y/rhs, self.z/rhs], self.conv_func)
	
	def __getitem__(self, key):
		return self.crds[key]
		
	def __setitem__(self, key, value):
		self.crds[key] = value
	
	def __getattr__(self, name):
		if name == 'x':
			return self.crds[0]
		elif name == 'y':
			return self.crds[1]
		elif name == 'z':
			return self.crds[2]
	
	def __setattr__(self, name, value):
		if name == 'x':
			self.crds[0] = value
		elif name == 'y':
			self.crds[1] = value
		elif name == 'z':
			self.crds[2] = value


class Projectable(object):
	def __init__(self, geometries=[]):
		self.proj_g = geometries.copy()
	
	def proj_geom(self, geometry):
		self.proj_g.append(geometry)
	
	def format_geom(self):
		if self.proj_g:
			geom_names = ' '.join(geom.name for geom in self.proj_g)
			return 'project ', f'({geom_names}) '
		else:
			return '',''


class Vertex(Point,Projectable):
	def __init__(self, crds, conv_func=cart_to_cart, name='', geometries=[]):
		Point.__init__(self, crds, conv_func)
		Projectable.__init__(self,geometries)
		self.name = name
		self.index = None
	
	def format(self):
		com = f'{self.index} {self.name}'
		vertex_str = Point.format(self)
		proj_str, geom_name = self.format_geom()
		return f'{proj_str}{vertex_str} {geom_name}// {com:s}'


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
	def __init__(self, center, radius, name):
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
	def __init__(self, point1, point2, radius, name):
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
	def __init__(self, vertices, name=''):

		#vertices is a 2x2 array
		self.vertices = vertices
		self.name = name
		self.proj_g = None
	
	def proj_geom(self, geometry):
		if self.proj_g:
			print(f'WARNING: face-{self.name} has already been projected to {self.proj_geom.name}, over-writing with {geometry.name}.')
		
		self.proj_g = geometry
	
	def is_projected(self):
		return (self.proj_g is not None)
	
	def format(self, write_proj=True):
		
		v_arr = self.vertices
		vts = [v_arr[0][0], v_arr[0][1], v_arr[1][0], v_arr[1][1]]
		
		index = ' '.join(str(v.index) for v in vts)
		com = ' '.join(v.name for v in vts)
		
		proj_str, geom_name = '',''
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
	def __init__(self, vertices, cells, name='', grading=uniformGrading):
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
		#http://www.openfoam.org/docs/user/mesh-description.php
		
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
	def __init__(self, vertices, arcMidPoint, name=''):
		
		Edge.__init__(self, vertices, name)
		self.arcMidPoint = arcMidPoint
	
	def format(self):
		return Edge.format(self).format('arc',self.arcMidPoint.format())


class SplineEdge(Edge):
	def __init__(self, vertices, points, name=''):
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


class ProjectionEdge(Edge,Projectable):
	def __init__(self, vertices, name='', geometries=[]):
		Edge.__init__(self, vertices, name)
		Projectable.__init__(self, geometries)
	
	def is_relevant(self):
		return bool(self.proj_g)
	
	def format(self):
		return Edge.format(self).format(*self.format_geom())


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
			buf.write(f'        {f.format(False)}\n')
		buf.write('    );\n')
		buf.write('}')
		return buf.getvalue()

def _format_section(name, secList):
	buf = io.StringIO()
	if name == 'geometry':
		brackets = ['{','}']
	else:
		brackets = ['(',')']
	
	buf.write(f'{name}\n')
	buf.write(f'{brackets[0]}\n')
	for item in secList:
		buf.write(f'    {item.format()}\n')
	buf.write(f'{brackets[1]};')
	
	return buf.getvalue()

class BlockMeshDict(object):
	def __init__(self):
		self.convert_to_meters = 1.0
		self.blocks = []
		self.edges = []
		self.boundaries = []
		self.geometries = []
		self.faces = []
	
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
	
	def add_hexblock(self, block):
		self.blocks.append(block)
	
	def add_edge(self, edge):
		self.edges.append(edge)
	
	def add_boundary(self, boundary):
		self.boundaries.append(boundary)
	
	def add_geometry(self, geometry):
		self.geometries.append(geometry)
	
	def add_face(self, face):
		self.faces.append(face)
	
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
	
	def format(self):
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

{_format_section('geometry',self.geometries)}

{_format_section('vertices',self.valid_vertices)}

{_format_section('edges',self.edges)}

{_format_section('blocks',self.blocks)}

{_format_section('faces',self.faces)}

{_format_section('boundary',self.boundaries)}

mergePatchPairs
(
);

// ************************************************************************* //
'''

'''
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  5.0.0                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
'''