from .core import *

import numpy as np

headers = ['vertices', 'num_divisions', 'grading', 'baked_vertices', 'edges', 'faces', 'block_mask', 'vertex_mask', 'edge_mask', 'face_mask']
formats = ['3f4','3u4','3O','O','3O','3O','?','?','3?','3?']
struct_type = np.dtype({'names' : headers, 'formats' : formats})

init_pos = np.arange(3)

def wrapRadians(values):
	return values % (2*np.pi)

def np_cyl_to_cart(crds):
	ncrds = crds.copy()
	ncrds[...,0] = np.multiply(crds[...,0],np.cos(crds[...,1]))
	ncrds[...,1] = np.multiply(crds[...,0],np.sin(crds[...,1]))
	return ncrds

def np_cart_to_cyl(crds):
	ncrds = crds.copy()
	ncrds[...,0] = np.linalg.norm(crds[...,:-1],axis=-1)
	ncrds[...,1] = np.arctan2(crds[...,1],crds[...,0])
	return ncrds

class BaseBlockStruct(object):
	
	def __init__(self, x0, x1, x2, nd0, nd1, nd2, conv_func=cart_to_cart, name=''):
		#Assume x0,x1,x2 are ascending 1D numpy arrays with dtype=np.float32, minimum 2 elements each
		#n0,n1,n2 are 1D numpy arrays of the number of divisions in each direction
		
		for x in [x0,x1,x2]:
			less_zero = np.diff(x) < 0.0
			if np.any(less_zero):
				print('ERROR -- A dimension array is in non-ascending order. Blocks will be inside out.')
				print(less_zero)
		
		shape = (x0.size,x1.size,x2.size)
		self.str_arr = np.empty(shape,dtype=struct_type)
		self.rshape = rshape = (shape[0]-1,shape[1]-1,shape[2]-1)
		
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
		
		#Initialize edges and faces
		for s in range(3):
			roll_pos = np.roll(init_pos,s)
			
			d_edges = np.moveaxis(self['edges'][...,s],init_pos,roll_pos)
			d_faces = np.moveaxis(self['faces'][...,s],init_pos,roll_pos)
			d_vts = np.moveaxis(self['baked_vertices'],init_pos,roll_pos)
			
			for i in range(rshape[s]):
				for j in range(shape[(s+1)%3]):
					for k in range(shape[(s+2)%3]):
						d_edges[i,j,k] = ProjectionEdge(d_vts[i:i+2,j,k])
			
			for i in range(shape[s]):
				for j in range(rshape[(s+1)%3]):
					for k in range(rshape[(s+2)%3]):
						d_faces[i,j,k] = Face(d_vts[i,j:j+2,k:k+2])
		
		self.name = name
	
	def project_structure(self, dir, face_ind, geometry):
		
		#Get the subarray relevant to the face being projected
		roll_pos = np.roll(init_pos,dir)
		struct = np.moveaxis(self.str_arr,init_pos,roll_pos)[face_ind]
		rshape = np.roll(np.array(self.rshape),-dir)[1:]
		shape = struct.shape
		
		#Project vertices
		b_vts = struct['baked_vertices']
		vt_mask = struct['vertex_mask']
		for ind in np.ndindex(shape):
			if not vt_mask[ind]:
				b_vts[ind].proj_geom(geometry)
		
		#Project edges
		edges = np.roll(struct['edges'],-dir,axis=-1)[...,1:]
		edge_mask = np.roll(struct['edge_mask'],-dir,axis=-1)[...,1:]
		for j in range(shape[0]):
			for k in range(shape[1]):
				for s in range(2):
					edge = edges[j,k,s]
					if (not edge_mask[j,k,s]) and isinstance(edge, ProjectionEdge):
						edge.proj_geom(geometry)
		
		#Project faces
		faces = struct['faces'][...,dir]
		face_mask = struct['face_mask'][...,dir]
		for j in range(rshape[0]):
			for k in range(rshape[1]):
				if not face_mask[j,k]:
					faces[j,k].proj_geom(geometry)
	
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
		
		#Reset block_mask edge
		bmask = self['block_mask']
		bmask[-1] = True
		bmask[:,-1] = True
		bmask[:,:,-1] = True
		
		for i in range(self.rshape[0]):
			for j in range(self.rshape[1]):
				for k in range(self.rshape[2]):
					
					if not bmask[i,j,k]:
						#Get subarray
						blockData = self[i:i+2,j:j+2,k:k+2]
						
						gt = blockData['grading'].copy()
						grading = self._get_grading(gt)
						
						nd = blockData['num_divisions'][0,0,0]
						
						vts = self._get_block_vertices(blockData['baked_vertices'])
						
						block_name = f'{self.name}-{i}-{j}-{k}'
						block = HexBlock(vts, nd, block_name, grading)
						block_mesh_dict.add_hexblock(block)
		
		#write relevant edges and faces
		shape = self.shape
		rshape = self.rshape
		
		for s in range(3):
			roll_pos = np.roll(init_pos,s)
			
			d_bmask = np.moveaxis(bmask,init_pos,roll_pos)
			d_edges = np.moveaxis(self['edges'][...,s],init_pos,roll_pos)
			d_faces = np.moveaxis(self['faces'][...,s],init_pos,roll_pos)
			
			for i in range(rshape[s]):
				for j in range(shape[(s+1)%3]):
					for k in range(shape[(s+2)%3]):
						edge = d_edges[i,j,k]
						
						if not edge.is_relevant():
							continue
						
						lmsk = d_bmask[i,max(j-1,0):j+1,max(k-1,0):k+1]
						if not np.all(lmsk):
							block_mesh_dict.add_edge(edge)
			
			for i in range(shape[s]):
				for j in range(rshape[(s+1)%3]):
					for k in range(rshape[(s+2)%3]):
						face = d_faces[i,j,k]
						
						if not face.is_projected():
							continue
						
						lmsk = d_bmask[max(i-1,0):i+1,j,k]
						if not np.all(lmsk):
							block_mesh_dict.add_face(face)
	
	#Default to underlying structured array
	def __getattr__(self, name):
		return getattr(self.str_arr, name)
	
	def __getitem__(self, key):
		return self.str_arr[key]


class CartBlockStruct(BaseBlockStruct):
	pass


class TubeBlockStruct(BaseBlockStruct):
	
	def __init__(self, rs, ts, zs, nr, nt, nz, name='', is_complete=False):
		
		if is_complete and ~np.isclose(wrapRadians(ts[0]),wrapRadians(ts[-1])):
			print(f'WARNING -- TubeBlockStruct {name} is marked as complete, while the first and last angles are unequal; make sure these are separated by 2*pi')
		
		BaseBlockStruct.__init__(self, rs, ts, zs, nr, nt, nz, cyl_to_cart, name)
		
		b_vts = self['baked_vertices']
		if is_complete:
			b_vts[:,-1] = b_vts[:,0]
		
		self.is_complete = is_complete
	
	def write(self, block_mesh_dict):
		
		shape = self.shape
		shp = tuple((shape[0],shape[1]-1,shape[2]))
		
		vts = self['vertices']
		b_vts = self['baked_vertices']
		
		cyls  = {}
		s_pt = Point([0,0,-1e5])
		e_pt = Point([0,0,1e5])
		for i,r in np.ndenumerate(np.unique(vts[...,0])):
			cyl = Cylinder(s_pt,e_pt,r,f'{self.name}-blockcyl-{i}')
			cyls[r] = cyl
			block_mesh_dict.add_geometry(cyl)
		
		b_vts = self['baked_vertices']
		vertex_mask = self['vertex_mask']
		proj_rcrds = self['vertices'][...,0]
		edges = self['edges'][...,1]
		edge_mask = self['edge_mask'][...,1]
		
		for ind in np.ndindex(shp):
			if not vertex_mask[ind]:
				b_vts[ind].proj_geom(cyls[proj_rcrds[ind]])
		
		for ind in np.ndindex(edges.shape):
			edge = edges[ind]
			if (not edge_mask[ind]) and isinstance(edge, ProjectionEdge):
				edge.proj_geom(cyls[proj_rcrds[ind]])
		
		BaseBlockStruct.write(self, block_mesh_dict)

dummy_vertex = Vertex(0,0,0)
dummy_edge = Edge([dummy_vertex]*2,name = 'dummy')

class CylBlockStructContainer(object):
	
	def __init__(self, rs, ts, zs, nr, nt, nz, name='', inner_arc_comp=0.25):
		
		self.tube_struct = TubeBlockStruct(rs, ts, zs, nr, nt, nz, name=name+'-tube', is_complete=True)
		
		self.inner_arc_comp = inner_arc_comp
		Ng = ((ts.size-1) // 4) + 1 #Assume integer number of divisions
		
		xs = np.linspace(-rs[0],rs[0],Ng)
		ys = xs.copy()
		
		nx = nt[:Ng].copy()
		ny = nt[Ng:2*Ng].copy()
		
		self.core_struct = CartBlockStruct(xs, ys, zs, nx, ny, nz, name=name+'-core')
		
		cyl_vts = np_cart_to_cyl(self.core_struct['vertices'])
		cyl_vts[...,1] -= 5/4*np.pi
		self.core_struct['vertices'][:] = np_cyl_to_cart(cyl_vts)
		
		core_b_vts = self.core_struct['baked_vertices']
		tube_b_vts = self.tube_struct['baked_vertices']
		
		#Connect the outer tube structure to the core
		tInds = np.arange(ts.size-1).reshape(4,Ng-1)
		
		for s in range(4):
			np.rot90(core_b_vts,k=-s)[:-1,0,:] = tube_b_vts[0,tInds[s],:]
	
	def write(self,block_mesh_dict):
		
		iac = self.inner_arc_comp
		if not np.isclose(iac,1.0):
			
			tube = self.tube_struct
			shape = tube.shape
			shp = tuple((shape[0],shape[1]-1,shape[2]))
			
			vts = tube['vertices'][0]
			b_vts = tube['baked_vertices'][0]
			edge_mask = tube['edge_mask'][0,...,1]
			tube['edges'][0,...,1:] = dummy_edge
			
			for ind in np.ndindex(shp[1:]):
				
				if edge_mask[ind]:
					continue
				
				end_pts = vts[ind[0]:ind[0]+2,ind[1]]
				end_vts = b_vts[ind[0]:ind[0]+2,ind[1]]
				mid_pt = Point((end_pts[0] + end_pts[1])/2,cyl_to_cart)
				sweep_angle = (end_pts[1,1] - end_pts[0,1])/2
				mid_pt.crds[0] *= iac*np.cos(sweep_angle) + (1-iac)
				block_mesh_dict.add_edge(ArcEdge(end_vts,mid_pt))
		
		self.tube_struct.write(block_mesh_dict)
		self.core_struct.write(block_mesh_dict)