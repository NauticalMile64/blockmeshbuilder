#Builds a structured tube mesh

import numpy as np
from ofblockmeshdicthelper import BlockMeshDict, TubeBlockStruct

bmd = BlockMeshDict()
bmd.set_metric('mm')

is_complete = True
rs = np.array([0.5,0.8,1.0])
ts = np.linspace(0,2*np.pi,9,endpoint=is_complete)
zs = np.array([0.0,0.5,1.5])

ndr = np.full_like(rs,6)
ndt = np.full_like(ts,5)
ndz = np.full_like(zs,8)

tube = TubeBlockStruct(rs,ts,zs,ndr,ndt,ndz,'ts',is_complete=is_complete)

tube.write(bmd)

with open(r'OF_case/system/blockMeshDict','w') as infile:
	infile.write(bmd.format())