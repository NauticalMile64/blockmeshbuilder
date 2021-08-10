from io import StringIO
import numpy as np
from .blockelements import Point
from .check_names import check_name


def _dict_format(dict_name, data_dict, indent_level=1):
	buf = StringIO()
	buf.write(dict_name)
	indent = '\t' * indent_level
	buf.write(f'\n{indent}{{')

	for label, data in data_dict.items():
		buf.write(f'\n{indent}\t{label} {data};')

	buf.write(f'\n{indent}}}')
	return buf.getvalue()


class Geometry:
	_unique_id = 0

	def __init__(self, name):
		check_name(name)
		self.name = f'{name}_id-{self._unique_id}'
		Geometry._unique_id += 1

	def do_format(self, data_dict):
		return _dict_format(self.name, data_dict)


class Plane(Geometry):
	plane_type = 'Plane Base Class'
	plane_data = {}

	def __init__(self, name='Plane'):
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
								   'point': Point.format(self.point),
								   'normal': Point.format(self.normal)
							   })


class PlaneEmbeddedPoints(Plane):
	plane_type = 'embeddedPoints'

	def __init__(self, point1, point2, point3, name):
		Plane.__init__(self, name)
		self.points = (point1, point2, point3)

	def format(self):
		return Plane.do_format(self,
							   {
								   'point1': Point.format(self.points[0]),
								   'point2': Point.format(self.points[1]),
								   'point3': Point.format(self.points[2])
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
	def __init__(self, center, radius, name='Sphere'):
		Geometry.__init__(self, name)
		self.center = center
		self.radius = radius

	def format(self):
		return Geometry.do_format(self,
								  {
									  'type': 'searchableSphere',
									  'centre': Point.format(self.center),
									  'radius': f'{self.radius:18.15g}'
								  })


class Cylinder(Geometry):
	def __init__(self, point1, point2, radius, name='Cylinder'):
		Geometry.__init__(self, name)
		self.point1 = point1
		self.point2 = point2
		self.radius = radius

	def format(self):
		return Geometry.do_format(self,
								  {
									  'type': 'searchableCylinder',
									  'point1': Point.format(self.point1),
									  'point2': Point.format(self.point2),
									  'radius': f'{self.radius:18.15g}'
								  })


class Cone(Geometry):
	small_value = np.finfo(np.float64).eps

	def __init__(self, point1, point2, radius1, radius2, name='Cone',
				 inner_radius1=small_value, inner_radius2=small_value):
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
									  'point1': Point.format(self.point1),
									  'point2': Point.format(self.point2),
									  'radius1': f'{self.radius1:18.15g}',
									  'radius2': f'{self.radius2:18.15g}',
									  'inner_radius1': f'{self.inner_radius1:18.15g}',
									  'inner_radius2': f'{self.inner_radius2:18.15g}'
								  })


_of_distribution_geometries = {
	'.org': {
		PlanePointAndNormal, PlaneEmbeddedPoints, PlaneEquation, Cylinder, Sphere,
		# closedTriSurfaceMesh, Box, Disc, ExtrudedCircle, Plate, SurfaceCollection, SurfaceWithGaps, triSurfaceMesh
	},
	'.com': {
		PlanePointAndNormal, PlaneEmbeddedPoints, PlaneEquation, Cylinder, Sphere, Cone
		# distributedTriSurfaceMesh, Box, Disc, ExtrudedCircle, Plate, SurfaceCollection, SurfaceWithGaps, triSurfaceMesh
	}
}
