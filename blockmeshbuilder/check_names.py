import re


class IllegalDictionaryNameError(Exception):
	pass


class BoundaryNameClashError(Exception):
	pass


_illegal_char_check = re.compile(r'[$*{};\'"/\\]').search


def check_name(name):
	if not isinstance(name, str):
		raise TypeError(f'The OpenFOAM dictionary name or zone label "{name}" is not a string. '
						f'Please provide a string.')

	if name == '':
		raise IllegalDictionaryNameError(f'An empty string was provided as a name for an OpenFOAM dictionary name or '
										 f'zone label. Please provide a name.')

	if name[0].isdigit():
		raise IllegalDictionaryNameError(f'The OpenFOAM dictionary name or zone label "{name}" contains a leading '
										 f'digit, which will yield an error during meshing with blockMesh. Please '
										 f'provide a name without any leading digits.')

	illegal_char_result = _illegal_char_check(name)
	if illegal_char_result:
		raise IllegalDictionaryNameError(f'The OpenFOAM dictionary name or zone label "{name}" contains an illegal '
										 f'character "{illegal_char_result[0]}", which will produce unexpected results '
										 f'or errors during meshing with blockMesh or running with OpenFOAM. Please '
										 f'provide a label without any of the following characters: $*{{}};\'"/\\')
