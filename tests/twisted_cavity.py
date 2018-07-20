#This example generates a blockMeshDict to build a square cavity mesh which can be used for running a lid-driven-cavity CFD problem

import numpy as np
from ofblockmeshdicthelper import BlockMeshDict, CartBlockStruct, SimpleGradingElement, Boundary

bmd = BlockMeshDict()
bmd.set_metric('mm')

xs = np.array([0.0,0.2,0.4,0.6,0.8,1.0]) - 0.5
ys = xs.copy()
zs = np.array([0.0,0.01])

ndx = np.full_like(xs,14)
ndy = ndx.copy()
ndz = np.array([1,0])

cavity = CartBlockStruct(xs,ys,zs,ndx,ndy,ndz,name='ts')

GD = cavity['grading']
edge_grd = 4
GD[:,0,:,1] = SimpleGradingElement(edge_grd)	#bottom
GD[:,-2,:,1] = SimpleGradingElement(1/edge_grd)	#top
GD[0,:,:,0] = SimpleGradingElement(edge_grd)	#left
GD[-2,:,:,0] = SimpleGradingElement(1/edge_grd)	#right

#Rotate vertices
vts = cavity['vertices']
XS,YS = vts[:,:,:,0],vts[:,:,:,1]
RS = np.sqrt(XS**2 + YS**2)
TS = np.arctan2(YS,XS)

TS[1:5,1:5,:] += np.pi/8 #Rotate interior blocks by 22.5 degrees
TS[2:4,2:4,:] += np.pi/8 #Rotate innermost block by a further 22.5 degrees

XS[:] = RS*np.cos(TS)
YS[:] = RS*np.sin(TS)

cavity.create_blocks()
cavity.write_blocks(bmd)

lid_bnd = Boundary('patch', 'lid')

for block in cavity['blocks'][:-1,-2,0]:
	lid_bnd.add_face(block.face('n'))

bmd.add_boundary(lid_bnd)

with open(r'OF_case/system/blockMeshDict','w') as infile:
	infile.write(bmd.format())