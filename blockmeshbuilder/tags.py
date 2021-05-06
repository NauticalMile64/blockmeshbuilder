from io import StringIO


class ZoneTag:
	def __init__(self, name):
		self.name = name


DEFAULT_ZONE_TAG = ZoneTag('DEFAULT')


class BoundaryTag:
	def __init__(self, name, type_='patch', info_dict=None):
		self.type_ = type_
		self.name = name
		if info_dict is None:
			info_dict = dict()
		info_dict['type'] = type_
		self.info_dict = info_dict

	def format(self):
		buf = StringIO()
		for keyword, data in self.info_dict.items():
			if not isinstance(data, str):
				data = data.format()
			buf.write(f'\t\t{keyword}\t{data};')
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
		# buf.write(f'    type {self.boundary_tag.type_};\n')
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
		'cyclic', 'cyclicACMI', 'cyclicAMI', 'cyclicRepeatAMI', 'cyclicSlip', 'empty', 'jumpCyclic',
		'jumpCyclicAMI', 'processor', 'processorCyclic', 'symmetry', 'symmetryPlane', 'wedge'
	},
	'.com': {
		'cyclic', 'cyclicACMI', 'cyclicAMI', 'cyclicSlip', 'empty', 'jumpCyclic', 'jumpCyclicAMI',
		'nonUniformTransformCyclic', 'processor', 'processorCyclic', 'symmetry', 'symmetryPlane', 'wedge'
	}
}
