import numpy as np


class GradingElement:
	pass


class SimpleGradingElement(GradingElement):
	def __init__(self, d):
		self.d = d

	def __neg__(self):
		return SimpleGradingElement(1 / self.d)

	def format(self):
		return str(self.d)


class MultiGradingElement(GradingElement):
	def __init__(self, len_pcts, cell_pcts, exp_ratios):
		self.len_pcts = np.asarray(len_pcts)
		self.cell_pcts = np.asarray(cell_pcts)
		self.exp_ratios = np.asarray(exp_ratios)

	def __neg__(self):
		return MultiGradingElement(self.len_pcts[::-1], self.cell_pcts[::-1], 1. / self.exp_ratios[::-1])

	def format(self):
		return '({})'.format(' '.join(f'({lp:.4f} {nc:.4f} {ex:.4f})'
									  for lp, nc, ex in zip(self.len_pcts, self.cell_pcts, self.exp_ratios)))


class Grading:
	def __init__(self, grading_elements):
		self.grading_elements = grading_elements

	def format(self):
		return '{{0}}Grading ({0})'.format(' '.join(ge.format() for ge in self.grading_elements))


class SimpleGrading(Grading):
	def __init__(self, grading_elements):
		assert (len(grading_elements) == 3)
		Grading.__init__(self, grading_elements)

	def format(self):
		return Grading.format(self).format('simple')


class EdgeGrading(Grading):
	def __init__(self, grading_elements):
		assert (len(grading_elements) == 12)
		Grading.__init__(self, grading_elements)

	def format(self):
		return Grading.format(self).format('edge')


# Helper function
def get_grading_info(len_pcts, dens):
	len_pcts = np.asarray(len_pcts)
	dens = np.asarray(dens)
	cum_lens = np.insert(np.cumsum(len_pcts), 0, 0.0)
	adens = np.array([np.trapz(dens[s:s + 2], cum_lens[s:s + 2]) for s in range(len_pcts.size)])
	exp_ratios = np.divide(dens[:-1], dens[1:])
	return len_pcts, adens, exp_ratios


uniformGradingElement = SimpleGradingElement(1)
uniformGrading = SimpleGrading([uniformGradingElement] * 3)
