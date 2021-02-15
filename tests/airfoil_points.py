import numpy as np


class NACA4:

	def __init__(self, foil, close_trailing_edge=False):
		m = float(foil[0]) / 100.  # max camber
		p = float(foil[1]) / 10.   # chord-wise position of max camber
		t = float(foil[2]) / 100.  # thickness
		self.foil = (m, p, t)
		if close_trailing_edge:
			self.c4 = -0.1036
		else:
			self.c4 = -0.1015

	def calculate_camber_line(self, x):
		x = np.asarray(x)
		m, p = self.foil[0], self.foil[1]
		if m:
			yc = np.asarray(2 * p * x - x ** 2)
			sel_a = x <= p
			yc[sel_a] *= m / p ** 2
			sel_b = np.logical_not(sel_a)
			yc[sel_b] = m / (1 - p) ** 2 * (1 - 2 * p + yc[sel_b])
			return yc
		else:
			return np.zeros_like(x)

	def calculate_surface_points(self, x):
		x = np.asarray(x)
		yt = 5 * self.foil[2] * (0.2969 * np.sqrt(x) - 0.126 * x - 0.3516 * x ** 2 + 0.2843 * x ** 3 + self.c4 * x ** 4)

		# If there is any camber
		m = self.foil[0]
		if m:
			yc = self.calculate_camber_line(x)

			p = self.foil[1]
			dydx = np.asarray(2 * m * (p - x) / p ** 2)
			dydx[x > p] *= p ** 2 / (1 - p) ** 2
			theta = np.arctan(dydx)

			yt_sin = yt * np.sin(theta)
			yt_cos = yt * np.cos(theta)
			xu = x - yt_sin
			yu = yc + yt_cos
			xl = x + yt_sin
			yl = yc - yt_cos
			return xu, yu, xl, yl

		else:
			return x, yt, x, -yt

	def calculate_leading_edge_radius(self):
		return 1.1019*(float(self.foil[2])/100.)**2


if __name__ == "__main__":

	NACA_code = (2, 4, 12)
	x_points = (1 - np.cos(np.linspace(0, np.pi, 100))) / 2
	a2412 = NACA4(NACA_code, close_trailing_edge=False)
	xu1, yu1, xl1, yl1 = a2412.calculate_surface_points(x_points)
	xu2, yu2, xl2, yl2 = NACA4(NACA_code, close_trailing_edge=True).calculate_surface_points(x_points)

	leading_edge_radius = 1.1019*(12./100)**2
	leading_edge_y = a2412.calculate_camber_line(leading_edge_radius)
	print(a2412.calculate_surface_points(0.05))

	import matplotlib.pyplot as plt

	plt.plot(xu1, yu1)
	plt.scatter(xl1, yl1)
	circle = plt.Circle((leading_edge_radius, leading_edge_y), leading_edge_radius, color='blue')
	plt.plot(xu2, yu2, color='orange')
	plt.scatter(xl2, yl2, color='orange')
	ax = plt.gca()
	ax.add_patch(circle)
	ax.set_aspect('equal')
	plt.show()
