from io import StringIO
from .blockelements import Point


def _rotational_specs(rotational, rotation_axis, rotation_centre):
	if rotational:
		if not (isinstance(rotation_axis, Point) and isinstance(rotation_centre, Point)):
			raise TypeError('Since the boundary is rotational, you must provide rotation_axis and rotation_centre '
							'objects of type blockmeshbuilder.Point')

		return {
			'transformType': 'rotational',
			'rotationAxis': rotation_axis,
			'rotationCentre': rotation_centre
		}
	else:
		return {}


def _match_tolerance_spec(match_tolerance):
	if match_tolerance is not None:
		return {'matchTolerance': match_tolerance}
	else:
		return {}


def _transform_patch_spec(transform_patch):
	if not isinstance(transform_patch, BoundaryTag):
		raise TypeError('repeat=True and something other than a BoundaryTag was passed as transform_patch. '
						 '\n A BoundaryTag must be passed to transform_patch.')
	return {'transformPatch': transform_patch.name}


class BoundaryTag:
	def __init__(self, name, type_='patch', info_dict=None):
		self.type_ = type_
		self.name = name
		if info_dict is None:
			info_dict = dict()
		info_dict['type'] = type_
		self.info_dict = info_dict

	@classmethod
	def empty_tag(cls, name):
		return cls(name, type_='empty')

	@classmethod
	def _cyclic_tags(cls, name_a, name_b, type_, common_dict):
		return cls(name_a, type_, {**{'neighbourPatch': name_b}, **common_dict}), \
			   cls(name_b, type_, {**{'neighbourPatch': name_a}, **common_dict})

	@classmethod
	def cyclic_tags(cls, name_a, name_b, rotational=False, rotation_axis=None, rotation_centre=None):
		common_dict = {}
		common_dict.update(_rotational_specs(rotational, rotation_axis, rotation_centre))
		return cls._cyclic_tags(name_a, name_b, 'cyclic', common_dict)

	@classmethod
	def cyclicAMI_tags(cls, name_a, name_b, rotational=False, rotation_axis=None, rotation_centre=None,
					   match_tolerance=None):
		common_dict = {}
		common_dict.update(_rotational_specs(rotational, rotation_axis, rotation_centre))
		common_dict.update(_match_tolerance_spec(match_tolerance))
		return cls._cyclic_tags(name_a, name_b, 'cyclicAMI', common_dict)

	@classmethod
	def cyclicRepeatAMI_tags(cls, name_a, name_b, rotational=False, rotation_axis=None, rotation_centre=None,
							 match_tolerance=None, transform_patch=None):
		common_dict = {}
		common_dict.update(_rotational_specs(rotational, rotation_axis, rotation_centre))
		common_dict.update(_match_tolerance_spec(match_tolerance))
		common_dict.update(_transform_patch_spec(transform_patch))
		return cls._cyclic_tags(name_a, name_b, 'cyclicRepeatAMI', common_dict)

	@classmethod
	def symmetry_tag(cls, name):
		return cls(name, 'symmetry')

	@classmethod
	def symmetryPlane_tag(cls, name):
		return cls(name, 'symmetryPlane')

	@classmethod
	def wedge_tag(cls, name):
		return cls(name, 'wedge')

	def format(self):
		buf = StringIO()
		for keyword, data in self.info_dict.items():
			if isinstance(data, Point):
				data = data.format()
			buf.write(f'\t\t{keyword}\t{data};\n')
		return buf.getvalue()


class _Boundary:
	def __init__(self, boundary_tag):
		assert (isinstance(boundary_tag, BoundaryTag))
		self.boundary_tag = boundary_tag
		self.faces = set()

	def add_face(self, face):
		self.faces.add(face)

	def format(self):
		buf = StringIO()
		buf.write(self.boundary_tag.name + '\n')
		buf.write('\t{\n')
		buf.write(self.boundary_tag.format())
		buf.write('\n\t\tfaces\n')
		buf.write('\t\t(\n')
		for f in self.faces:
			buf.write(f'\t\t\t{f.format(False)}\n')
		buf.write('\t\t);\n')
		buf.write('\t}')
		return buf.getvalue()


_of_distribution_constraints = {
	'.org': {
		'cyclic',  'cyclicAMI', 'cyclicRepeatAMI',  'empty', 'symmetry', 'symmetryPlane', 'wedge'
		# 'cyclicACMI', 'cyclicSlip', 'jumpCyclic', 'jumpCyclicAMI', 'processor', 'processorCyclic',
	},
	'.com': {
		'cyclic', 'cyclicAMI', 'empty', 'symmetry', 'symmetryPlane', 'wedge'
		# 'cyclicACMI', 'cyclicSlip', 'jumpCyclic', 'jumpCyclicAMI', 'processor', 'processorCyclic',
		# 'nonUniformTransformCyclic'
	}
}
