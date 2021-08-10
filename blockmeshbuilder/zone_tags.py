from .check_names import check_name


class ZoneTag:
	def __init__(self, name):
		check_name(name)
		self.name = name


DEFAULT_ZONE_TAG = ZoneTag('DEFAULT')
