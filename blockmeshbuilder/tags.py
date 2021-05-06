from io import StringIO


class ZoneTag:
	def __init__(self, name):
		self.name = name


DEFAULT_ZONE_TAG = ZoneTag('DEFAULT')


class BoundaryTag:
	def __init__(self, name, type_='patch'):
		self.type_ = type_
		self.name = name


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
		buf.write('{\n')
		buf.write(f'    type {self.boundary_tag.type_};\n')
		buf.write('    faces\n')
		buf.write('    (\n')
		for f in self.faces:
			buf.write(f'        {f.format(False)}\n')
		buf.write('    );\n')
		buf.write('}')
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
