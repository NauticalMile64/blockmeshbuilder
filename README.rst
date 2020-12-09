
ofblockmeshdicthelper
=============================
This module helps to automate the generation of hexahedral block-structured meshes using OpenFOAMs blockMesh tool.

How does it work?
---
The user writes a python script using ofblockmeshdicthelper to create a clear, descriptive representation of the desired block structure(s). At the end of the python script, the user can write this information into a blockMeshDict file with a single line of code. Then the user can run OpenFOAM's blockMesh tool to generate the mesh.

How is this better than writing the blockMeshDict file myself?
---
blockMeshDict files contain lists of geometries, blocks, faces, etc... which all reference a master list of vertices by index or by name. This is reasonable for simple geometries, but quickly becomes untennable for complex ones. I have even found it difficult to debug blockMeshDict files with two blocks.

Installation
---
Clone this repository and execute `python -m pip install .` for the regular installation, or `python -m pip install -e .` for developer mode.

Example
---
Here is an example which generate wedged model shown at 
https://openfoamwiki.net/index.php/Main_ContribExamples/AxiSymmetric

.. code-block:: python

	""" example of ofblockmeshdicthelper
	Reproduce the wedge shown at
	https://openfoamwiki.net/index.php/Main_ContribExamples/AxiSymmetric
	"""
	import numpy as np
	from math import radians
	from ofblockmeshdicthelper import BlockMeshDict, TubeBlockStruct, Boundary

	# wedge dimensions
	wd = radians(10.0)
	r = 0.19
	l = 1.1

	# prepare ofblockmeshdicthelper.BlockMeshDict instance to gather vertices, blocks, faces and boundaries.
	bmd = BlockMeshDict()

	# set metrics
	bmd.set_metric('mm')

	rs = np.array([0,r])
	ts = np.array([-wd/2,wd/2])
	zs = np.array([0,1])

	nr = np.array([10])
	nt = np.array([1])
	nz = np.array([10])

	wedge = TubeBlockStruct(rs,ts,zs,nr,nt,nz,'wedge')
	wFaces = wedge['faces']

	#Front and back boundaries
	front_bnd = Boundary('patch', 'front', faces = [wFaces[0,0,0,2]])
	bmd.add_boundary(front_bnd)

	back_bnd = Boundary('patch', 'back', faces = [wFaces[0,0,1,2]])
	bmd.add_boundary(back_bnd)

	wedge.write(bmd)

	# output
	with open(r'OF_case/system/blockMeshDict','w') as infile:
	infile.write(bmd.format())
