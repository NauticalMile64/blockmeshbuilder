from .core import *

import numpy as np

headers = ['vertices', 'num_divisions', 'grading', 'baked_vertices', 'faces', 'proj_vts', 'proj_edges', 'proj_faces', 'block_mask']
formats = ['3f4','3u4','3O','O','3O','O','3O','3O','?']
struct_type = np.dtype({'names' : headers, 'formats' : formats})

class BaseBlockStruct(object):
	
	def __init__(self, x0, x1, x2, nd0, nd1, nd2, conv_func=cart_to_cart, name=''):
		#Assume x0,x1,x2 are ascending 1D numpy arrays with dtype=np.float32, minimum 2 elements each
		#n0,n1,n2 are 1D numpy arrays of the number of divisions in each direction
		
		shape = (x0.size,x1.size,x2.size)
		self.str_arr = np.empty(shape,dtype=struct_type)
		rshape = (shape[0]-1,shape[1]-1,shape[2]-1)
		self.rshape = rshape
		
		#Initialize vertices
		X0,X1,X2 = np.meshgrid(x0,x1,x2,indexing='ij')
		
		vts = self['vertices']
		vts[...,0] = X0
		vts[...,1] = X1
		vts[...,2] = X2
		
		for ind in np.ndindex(self.shape):
			self['baked_vertices'][ind] = Vertex(vts[ind],conv_func)
		
		#Initialize number of divisions
		ND0,ND1,ND2 = np.meshgrid(nd0,nd1,nd2,indexing='ij')
		
		nds = self['num_divisions']
		nds[...,0] = ND0
		nds[...,1] = ND1
		nds[...,2] = ND2
		
		#Initialize grading
		self['grading'][:] = uniformGradingElement
		
		#Initialize faces
		init_pos = np.arange(3)
		for s in range(3):
			roll_pos = np.roll(init_pos,s)
			d_faces = np.moveaxis(self['faces'][...,s],init_pos,roll_pos)
			d_vts = np.moveaxis(self['baked_vertices'],init_pos,roll_pos)
			for i in range(shape[s]):
				for j in range(rshape[(s+1)%3]):
					for k in range(rshape[(s+2)%3]):
						d_faces[i,j,k] = Face(d_vts[i,j:j+2,k:k+2])
		
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
	
	def write(self,block_mesh_dict):
		
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
						block_mesh_dict.add_hexblock(block)
		
		#Project vertices, edges, and faces
		shape = self.shape
		
		#Project vertices
		proj_vts = self['proj_vts']
		if np.any(proj_vts):
			for ind in np.ndindex(shape):
				if proj_vts[ind]:
					self['baked_vertices'][ind].proj_geom(proj_vts[ind])
		
		#Edges and Faces
		if not (np.any(self['proj_edges']) or np.any(self['proj_faces'])):
			return
		
		rshape = self.rshape
		init_pos = np.arange(3)
		for s in range(3):
			roll_pos = np.roll(init_pos,s)
			d_pe = np.moveaxis(self['proj_edges'][...,s],init_pos,roll_pos)
			d_pf = np.moveaxis(self['proj_faces'][...,s],init_pos,roll_pos)
			
			d_faces = np.moveaxis(self['faces'][...,s],init_pos,roll_pos)
			d_vts = np.moveaxis(self['baked_vertices'],init_pos,roll_pos)
			d_blkmsk = np.moveaxis(self['block_mask'],init_pos,roll_pos)
			
			#Project edges
			for i in range(rshape[s]):
				for j in range(shape[(s+1)%3]):
					for k in range(shape[(s+2)%3]):
						if d_pe[i,j,k]:
							block_mesh_dict.add_edge(ProjectionEdge(d_vts[i:i+2,j,k],d_pe[i,j,k]))
			
			#Project faces
			for i in range(shape[s]):
				for j in range(rshape[(s+1)%3]):
					for k in range(rshape[(s+2)%3]):
						if d_pf[i,j,k]:
							d_faces[i,j,k].proj_geom(d_pf[i,j,k])
							block_mesh_dict.add_face(d_faces[i,j,k])
	
	
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
	
	def __init__(self, rs, ts, zs, nr, nt, nz, name='', is_complete=False, inner_arc_comp=0.0):
		
		if is_complete and ~np.isclose(wrapRadians(ts[0]),wrapRadians(ts[-1])):
			print(f'WARNING -- TubeBlockStruct {name} is marked as complete, while the first and last angles are unequal; make sure these are separated by 2*pi')
		
		BaseBlockStruct.__init__(self, rs, ts, zs, nr, nt, nz, cyl_to_cart, name)
		
		b_vts = self['baked_vertices']
		if is_complete:
			b_vts[:,-1] = b_vts[:,0]
		
		self.is_complete = is_complete
		self.inner_arc_comp = inner_arc_comp
	
	def write(self, block_mesh_dict):
		
		shape = self.shape
		shp = tuple((shape[0],shape[1]-1,shape[2]))
		iac = self.inner_arc_comp
		skip_inner_arc = np.isclose(iac,1.0)
		
		vts = self['vertices']
		b_vts = self['baked_vertices']
		ang_vts = vts[...,1]
		for ind in np.ndindex(shp):
			end_pts = vts[ind[0],ind[1]:ind[1]+2,ind[2]]
			end_vts = b_vts[ind[0],ind[1]:ind[1]+2,ind[2]]
			mid_pt = Point((end_pts[0] + end_pts[1])/2,cyl_to_cart)
			if ind[0] == 0:
				if skip_inner_arc:
					continue
				
				sweep_angle = (end_pts[1,1] - end_pts[0,1])/2
				print(sweep_angle)
				mid_pt.crds[0] *= iac*np.cos(sweep_angle) + (1-iac)
			
			block_mesh_dict.add_edge(ArcEdge(end_vts,mid_pt))
		
		BaseBlockStruct.write(self, block_mesh_dict)