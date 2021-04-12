import os
from setuptools import setup, find_packages

ver_file = os.path.join('blockmeshbuilder', 'version.py')
vars = {}
exec(open(ver_file).read(), vars)

setup(
	name='blockmeshbuilder',
	version=vars['__version__'],
	description='Helper utilities for OpenFOAM blockMeshDict generation.',
	long_description=open('README.md', 'rt').read(),
	author='Nolan Dyck',
	url='https://github.com/NauticalMile64/blockmeshbuilder',
	packages=find_packages(),
	package_dir={'blockmeshbuilder': 'blockmeshbuilder'},
	install_requires=['numpy'],
	classifiers=[
		"Programming Language :: Python :: 3.6",
		"Programming Language :: Python :: 3.7",
		"Programming Language :: Python :: 3.8",
		"Programming Language :: Python :: 3.9",
		"Programming Language :: Python :: 3.10",
		"License :: OSI Approved :: MIT License",
		"Development Status :: 4 - Beta",
		"Environment :: Other Environment",
		"Intended Audience :: Science/Research",
		"Operating System :: OS Independent",
		"Topic :: Scientific/Engineering :: Physics"])
