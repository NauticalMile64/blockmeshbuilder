from io import StringIO
import numpy as np
import warnings
from .tags import DEFAULT_ZONE_TAG
from .grading import uniformGrading


def cart_to_cart(crds):
	return np.asarray(crds).copy()


def cyl_to_cart(crds):
	crds = np.asarray(crds, dtype=np.floating)
	ncrds = crds.copy()
	ncrds[..., 0] = np.multiply(crds[..., 0], np.cos(crds[..., 1]))
	ncrds[..., 1] = np.multiply(crds[..., 0], np.sin(crds[..., 1]))
	return ncrds


def cart_to_cyl(crds):
	crds = np.asarray(crds, dtype=np.floating)
	ncrds = crds.copy()
	ncrds[..., 0] = np.hypot(crds[..., 0], crds[..., 1])
	ncrds[..., 1] = np.arctan2(crds[..., 1], crds[..., 0])
	return ncrds


cart_conv_pair = (cart_to_cart, cart_to_cart)
cyl_conv_pair = (cyl_to_cart, cart_to_cyl)


class Point:
	def __init__(self, crds, conv_funcs=cart_conv_pair):
		self.crds = np.asarray(crds, dtype=np.float64)
		self.conv_funcs = conv_funcs

	def get_cart_crds(self):
		return self.conv_funcs[0](self.crds)

	def format(self):
		ccrds = self.get_cart_crds()
		return f'( {ccrds[0]:18.15g} {ccrds[1]:18.15g} {ccrds[2]:18.15g} )'

	@staticmethod
	def _get_oth_coordinates(rhs):
		return rhs.get_cart_crds() if issubclass(type(rhs), Point) else rhs

	def __add__(self, rhs):
		return type(self)(self.conv_funcs[1](self.get_cart_crds() + Point._get_oth_coordinates(rhs)), self.conv_funcs)

	def __radd__(self, lhs):
		return self + lhs

	def __sub__(self, rhs):
		return type(self)(self.conv_funcs[1](self.get_cart_crds() - Point._get_oth_coordinates(rhs)), self.conv_funcs)

	def __rsub__(self, lhs):
		return type(self)(self.conv_funcs[1](Point._get_oth_coordinates(lhs) - self.get_cart_crds()), self.conv_funcs)

	def __mul__(self, rhs):
		return type(self)(self.conv_funcs[1](self.get_cart_crds() * Point._get_oth_coordinates(rhs)), self.conv_funcs)

	def __rmul__(self, lhs):
		return self * lhs

	def __pow__(self, rhs):
		return type(self)(self.conv_funcs[1](self.get_cart_crds() ** rhs), self.conv_funcs)

	def __truediv__(self, rhs):
		return type(self)(self.conv_funcs[1](self.get_cart_crds() / Point._get_oth_coordinates(rhs)), self.conv_funcs)

	@staticmethod
	def _convert_oth_coordinates(rhs, conv_func):
		return conv_func(rhs.get_cart_crds()) if issubclass(type(rhs), Point) else rhs

	def __iadd__(self, rhs):
		self.crds += Point._convert_oth_coordinates(rhs, self.conv_funcs[1])
		return self

	def __isub__(self, rhs):
		self.crds -= Point._convert_oth_coordinates(rhs, self.conv_funcs[1])
		return self

	def __imul__(self, rhs):
		self.crds *= Point._convert_oth_coordinates(rhs, self.conv_funcs[1])
		return self

	def __itruediv__(self, rhs):
		self.crds /= Point._convert_oth_coordinates(rhs, self.conv_funcs[1])
		return self

	def __lt__(self, rhs):
		return self.crds < Point._convert_oth_coordinates(rhs, self.conv_funcs[1])

	def __gt__(self, lhs):
		return lhs < self

	def __neg__(self):
		return type(self)(self.conv_funcs[1](-self.get_cart_crds()), self.conv_funcs)

	def __abs__(self):
		return abs(self.get_cart_crds())

	def __getitem__(self, key):
		return self.crds[key]

	def __setitem__(self, key, value):
		self.crds[key] = value


class Projectable:
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


class Edge:
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


class Face:
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


class HexBlock:
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

		if isinstance(index, str):
			index = kw_to_index[index]

		face_vertices = [[self.vertices[i] for i in row] for row in index_to_vertex[index]]

		if name is None:
			name = index_to_defaultsuffix[index].format(self.name)

		return Face(face_vertices, name)
