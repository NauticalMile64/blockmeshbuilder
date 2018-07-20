from .core import *

import numpy as np

headers = ['vertices', 'num_divisions', 'grading', 'baked_vertices', 'block_mask', 'blocks']
formats = ['3f4','3u4','3O','O','?','O']
dtype_dict = {'names' : headers, 'formats' : formats}
struct_type = np.dtype(dtype_dict)

class BaseBlockStruct:
	
	def __init__(self, x0, x1, x2, nd0, nd1, nd2, conv_func=cart_to_cart, name=''):
		#Assume x0,x1,x2 are ascending 1D numpy arrays with dtype=np.float32, minimum 2 elements each
		#n0,n1,n2 are 1D numpy arrays of the number of divisions in each direction
		
		shape = (x0.size,x1.size,x2.size)
		self.str_arr = np.empty(shape,dtype=struct_type)
		self.rshape = (shape[0]-1,shape[1]-1,shape[2]-1)
		
		#Initialize vertices
		X0,X1,X2 = np.meshgrid(x0,x1,x2,indexing='ij')
		
		vts = self['vertices']
		vts[...,0] = X0
		vts[...,1] = X1
		vts[...,2] = X2
		
		b_vts = self['baked_vertices']
		for ind in np.ndindex(self.shape):
			b_vts[ind] = Vertex(vts[ind],conv_func)
		
		#Initialize number of divisions
		ND0,ND1,ND2 = np.meshgrid(nd0,nd1,nd2,indexing='ij')
		
		nds = self['num_divisions']
		nds[...,0] = ND0
		nds[...,1] = ND1
		nds[...,2] = ND2
		
		#Initialize grading
		self['grading'][:] = uniformGradingElement
		
		self.name = name
	
	@staticmethod
	def _get_grading(gt):
		
		#Get relevant edges
		grd_arr = np.array([
			np.moveaxis(gt[...,s],s,0)[0].T for s in range(3)
		])
		
		#Get simplest grading type
		if np.all(grd_arr == uniformGradingElement):
			return uniformGrading
			
		elif np.all(np.array([grd_arr[s] == grd_arr[s,0,0] for s in range(3)])):
			return SimpleGrading(grd_arr[:,0,0])
			
		else:
			for gc in grd_arr:
				gc[1,0],gc[1,1] = gc[1,1],gc[1,0]
			
			return EdgeGrading(grd_arr.flatten())
	
	@staticmethod
	def _get_block_vertices(c_vs):
		return tuple((c_vs[0,0,0],c_vs[1,0,0],
			c_vs[1,1,0],c_vs[0,1,0],
			
			c_vs[0,0,1],c_vs[1,0,1],
			c_vs[1,1,1],c_vs[0,1,1]))
	
	def create_blocks(self):
		
		#Not very elegant right now, maybe I'll replace these loops with something more pythonic later
		for i in range(self.rshape[0]):
			for j in range(self.rshape[1]):
				for k in range(self.rshape[2]):
					
					if not self['block_mask'][i,j,k]:
						#Get subarray
						blockData = self[i:i+2,j:j+2,k:k+2]
						
						gt = blockData['grading'].copy()
						grading = self._get_grading(gt)
						
						nd = blockData['num_divisions'][0,0,0]
						
						vts = self._get_block_vertices(blockData['baked_vertices'])
						
						block_name = f'{self.name}-{i}-{j}-{k}'
						
						block = HexBlock(vts, nd, block_name, grading)
						blockData['blocks'][0,0,0] = block
	
	def write_blocks(self,block_mesh_dict):
		
		block_arr = self['blocks']
		for ind in np.ndindex(self.rshape):
			block = block_arr[ind]
			if block:
				block_mesh_dict.add_hexblock(block,f'{self.name} {ind}')
	
	#Default to underlying structured array
	def __getattr__(self, name):
		return getattr(self.str_arr, name)
	
	def __getitem__(self, key):
		return self.str_arr[key]


class CartBlockStruct(BaseBlockStruct):
	pass


def wrapRadians(values):
	return values % (2*np.pi)


class TubeBlockStruct(BaseBlockStruct):
	
	def __init__(self, rs, ts, zs, nr, nt, nz, name='', is_complete=False):
		
		if is_complete and ~np.isclose(wrapRadians(ts[0]),wrapRadians(ts[-1])):
			print(f'WARNING -- TubeBlockStruct {name} is marked as complete, while the first and last angles are unequal; make sure these are separated by 2*pi')
		
		BaseBlockStruct.__init__(self, rs, ts, zs, nr, nt, nz, cyl_to_cart, name)
		
		b_vts = self['baked_vertices']
		if is_complete:
			b_vts[:,-1] = b_vts[:,0]
		
		self.is_complete = is_complete
