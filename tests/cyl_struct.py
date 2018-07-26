#Builds a structured O-grid mesh

import numpy as np
from ofblockmeshdicthelper import BlockMeshDict, CylBlockStructContainer

bmd = BlockMeshDict()
bmd.set_metric('mm')

rs = np.array([0.3,0.6,1.0])
ts = np.linspace(0,2*np.pi,9,endpoint=True)
zs = np.array([0.0,0.5,1.5])

ndr = np.full_like(rs,6)
ndt = np.full_like(ts,6)
ndz = np.full_like(zs,8)

tube = CylBlockStructContainer(rs,ts,zs,ndr,ndt,ndz,'ts')

tube.write(bmd)

with open(r'OF_case/system/blockMeshDict','w') as infile:
	infile.write(bmd.format())