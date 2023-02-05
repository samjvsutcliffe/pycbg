import gmsh, os, json, pickle, csv, runpy, sys, shutil
import multiprocessing
import numpy as np
import itertools as it  
from pycbg import __version__ as pycbg_version

_ED2Q36BS_nodes = {
            0 :(-2,-2), 1 :(-1,-2), 2 :(0,-2), 3 :(1,-2), 4 :(2,-2), 5 :(3,-2),
            11 :(-2,-1), 10 :(-1,-1), 9 :(0,-1), 8 :(1,-1), 7:(2,-1), 6:(3,-1),
            12:(-2, 0), 13:(-1, 0), 14:(0, 0), 15:(1, 0), 16:(2, 0), 17:(3, 0),
            23:(-2, 1), 22:(-1, 1), 21:(0, 1), 20:(1, 1), 19:(2, 1), 18:(3, 1),
            24:(-2, 2), 25:(-1, 2), 26:(0, 2), 27:(1, 2), 28:(2, 2), 29:(3, 2),
            35:(-2, 3), 34:(-1, 3), 33:(0, 3), 32:(1, 3), 31:(2, 3), 30:(3, 3)
        }

class Mesh():
    """Create and write to a file a mesh using gmsh.

    Parameters
    ----------
    dimensions : tuple of floats
        Dimensions of the mesh. Its length should be 3, with `dimensions[n]` the dimension of the mesh on the axis `n`.
    ncells : tuple of ints
        Number of cells in each direction. Its length should be 3, with `ncells[n]` the number of cells on the axis `n`.
    origin : tuple of floats
        Origin of the mesh. Default is `(0.,0.,0.)`.
    directory : str, optional
        Directory in which the mesh file will be saved. If the directory doesn't already exist, it will be created. It is set by default to the current working directory.
    check_duplicates : bool, optional  
        See CB-Geo documentation for informations on this parameter. Default is `True`.
    cell_type : {'ED3H8', 'ED3H20', 'ED3H64G'}, optional
        Type of cell. Only 3D Hexahedrons are supported. The number of nodes can be 8, 20 or 64. Default is 'ED3H8'.

    Attributes
    ----------
    nodes : numpy array
        Positions of all nodes in the mesh. The id of a node is the index of its line in this array. Noting `nnodes` the number of nodes, the shape of `nodes` is ``(nnodes,3)``.
    cells : numpy array
        Connections between cells and nodes. Each line corresponds to a cell, its index is the cell's id. The columns correspond to the ids of the nodes composing a cell. Noting `nnode_pcell` the number of nodes per cell (8, 20 or 64), the shape of `cells` is ``(ncells,nnode_pcell)``.
    dimensions : tuple of floats
        Dimensions of the mesh.
    l0, l1, l2 : floats
        Dimensions of the mesh (``self.l0, self.l1, self.l2 = self.dimensions``).
    ncells : tuple of ints
        Number of cells in each direction.
    nc0, nc1, nc2 : ints
        Number of cells in each direction (``self.nc0, self.nc1, self.nc2 = self.ncells``).
    origin : tuple of floats
        Origin of the mesh.
    directory : str
        Directory in which the mesh file will be saved.
    check_duplicates : bool
        See CB-Geo documentation.
    cell_type : {'ED3H8', 'ED3H20', 'ED3H64G', 'ED2Q4', 'ED2Q8', 'ED2Q9', 'ED2Q16G', 'ED2Q36BS'}
        Type of cell. 
    n_dims : int
        Number of dimensions (2 for 2D and 3 for 3D), automatically determined from the cell type.
    round_decimal : int or None
        Rounds nodes coordinates to the specified decimal (`round_decimal` is directly passed to the built-in function `round`). This is useful when using `ED3H64G` or 'ED2Q16G'. Default to None. 
    params : dict
        Dictionary containing all necessary parameters (except the `directory` parameter) to create a copy of this `Mesh` object. For instance, one can do: `mesh_copy = Mesh(**existing_mesh.params)`.

    Notes
    -----
     - The mesh file is written upon creating the object.
     
    Examples
    --------
    Creating a cubic mesh of 1000 cells : 

    >>> mesh = Mesh((1.,1.,1.), (10,10,10))
    >>> mesh.nc0 * mesh.nc1 * mesh.nc2
    1000
    """
    ## TODO: - Test 'ED3H20' and 'ED3H64G'
    ##       - Avoid having to write the mesh file from gmsh for rewritting it again
    ##       - Make crete_mesh usable by the user 

    def __init__(self, dimensions, ncells, origin=(0.,0.,0.), directory="", check_duplicates=True, cell_type="ED3H8", round_decimal=None):
        self.cell_type = cell_type
        if cell_type=='ED3H8': self.nn_percell, self.n_dims = 8, 3
        elif cell_type=='ED3H20': self.nn_percell, self.n_dims = 20, 3
        elif cell_type=='ED3H64G': self.nn_percell, self.n_dims = 64, 3
        elif cell_type=='ED2Q4': self.nn_percell, self.n_dims = 4, 2
        elif cell_type=='ED2Q8': self.nn_percell, self.n_dims = 8, 2
        elif cell_type=='ED2Q9': self.nn_percell, self.n_dims = 9, 2
        elif cell_type=='ED2Q16G': self.nn_percell, self.n_dims = 16,2 #16, 2
        elif cell_type=='ED2Q36BS': self.nn_percell, self.n_dims = 36, 2
        elif cell_type=='ED2GIMP': self.nn_percell, self.n_dims = 4,2 #16, 2
        else : raise ValueError("cell_type is set to '{:}' while it should be one of the following: 'ED3H8', 'ED3H20', 'ED3H64G', 'ED2Q4', 'ED2Q8', 'ED2Q9', 'ED2Q16G' or 'ED2Q36BS'".format(cell_type))
        
        if self.n_dims == 2 and origin == (0.,0.,0.): origin = (0.,0.)
        self.set_parameters(dimensions, ncells, origin)
        if directory == '' : directory = './'
        if directory[-1] != '/' : directory += '/'
        if not os.path.isdir(directory): os.mkdir(directory)
        self.directory = directory

        self.check_duplicates = check_duplicates

        self._isoparametric = False # Shouldn't have to be set to another value
        self._io_type = "Ascii{:d}D".format(self.n_dims) # Shouldn't have to be set to another value
        self._node_type = "N{:d}D".format(self.n_dims) # Shouldn't have to be set to another value
        self.round_decimal = round_decimal

        self.dimensions, self.ncells = dimensions, ncells
        self.__reset_params()
        if cell_type == 'ED2GIMP':
            gen_gimp_mesh(self.dimensions,self.ncells,None,self)
        else:
            self.write_file()

    def set_parameters(self, dimensions, ncells, origin):
        """Set the dimensions and number of cells of the mesh.

        Parameters
        ----------
        dimensions : tuple of floats
            Dimensions of the mesh. Its length should be 3, with `dimensions[n]` the dimension of the mesh on the axis `n`.
        ncells : tuple of ints
            Number of cells in each direction. Its length should be 3, with `ncells[n]` the number of cells on the axis `n`.
        origin : tuple of floats
            Origin of the mesh. Default is `(0.,0.,0.)`.
        """
        if self.n_dims==2: 
            self.l0, self.l1 = dimensions
            self.nc0, self.nc1 = ncells
        elif self.n_dims==3: 
            self.l0, self.l1, self.l2 = dimensions
            self.nc0, self.nc1, self.nc2 = ncells
        else: raise RuntimeError("Number of dimensions coulnd't be detected, please check if `cell_type` is correctly set.")

        self.origin = origin

    def create_mesh(self):
        """Create the mesh in gmsh.

        Notes
        -----
         - This method calls `gmsh.initialize` but doesn't call `gmsh.finalize``
         - `cells` and `nodes` attributes are not created by this method
         - User shouldn't have to use this method as it is called by `write_file`
        """
        gmsh.initialize()
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.1)
        
        if self.n_dims==2: origin = list(self.origin) + [0]
        elif self.n_dims==3: origin = list(self.origin)
        else: raise RuntimeError("Number of dimensions coulnd't be detected, please check if `cell_type` is correctly set.")


        p = gmsh.model.geo.addPoint(*origin)
        l = gmsh.model.geo.extrude([(0, p)], self.l0, 0, 0, [self.nc0], [1])
        s = gmsh.model.geo.extrude([l[1]], 0, self.l1, 0, [self.nc1], [1], recombine=True)
        if self.n_dims==3:
            v = gmsh.model.geo.extrude([s[1]], 0, 0, self.l2, [self.nc2], [1], recombine=True)
            group = v[1][1]
        elif self.n_dims==2: group = s[1][1]
        gmsh.model.geo.synchronize()
        gmsh.model.addPhysicalGroup(self.n_dims, [group])
        gmsh.model.mesh.generate(self.n_dims)
        if self.cell_type in ['ED3H20', 'ED2Q8', 'ED2Q9']: 
            if self.cell_type!='ED2Q9': gmsh.option.setNumber("Mesh.SecondOrderIncomplete", 1)
            gmsh.model.mesh.setOrder(2)
        elif self.cell_type in ['ED3H64G', 'ED2Q16G']: 
        #elif self.cell_type in ['ED3H64G']: 
            gmsh.model.mesh.setOrder(3)
        elif self.cell_type in ['ED2GIMP']: 
            gmsh.model.mesh.setOrder(2)

    def write_file(self, filename="mesh"):
        """
        Write the mesh file formated for CB-Geo MPM.
        
        Parameters
        ----------
        filename : str, optional
            Name of the mesh file, the extension '.txt' is automatically added. Default is 'mesh'.
        """
        self._gmsh_filename = self.directory + filename + ".msh"
        self._filename = self.directory + filename + ".txt"
        self.create_mesh()
        
        
        gmsh.write(self._gmsh_filename)
        gmsh.finalize()
        self.__reformat_from_gmsh()
        self.__reset_params()
        self.cells, self.nodes = np.array(self.cells), np.array(self.nodes)
        if self.cell_type == 'ED2GIMP':
            print("Creating GIMP mesh")
            #print("outputdir:{}".format(self.directory))
            gen_gimp_mesh(self.dimensions,self.ncells,self.directory,self)
    
    def __reformat_from_gmsh(self): # Not meant for the user
        """Reads the mesh file generated by gmsh and reformats it for CB-
        Geo."""
        with open(self._gmsh_filename, 'r') as fil: lines = fil.readlines()
        os.remove(self._gmsh_filename)
        for i, line in enumerate(lines):
            if "$Nodes" in line : nn, start_nodes = int(float(lines[i+1])), i+2
            elif "$EndNodes" in line : end_nodes = i
            elif "$Elements" in line : ne, start_ele = int(float(lines[i+1])), i+2
            elif "$EndElements" in line : end_ele = i
        
        self.cells, self.nodes = [], []
        with open(self._filename, 'w') as fil:
            fil.write("{:d}\t{:d}\n".format(nn, ne))
            for line in lines[start_nodes:end_nodes]: 
                sl = line.split(' ')
                def wrapped_round(x): return round(x, self.round_decimal) if self.round_decimal is not None else x
                node = [wrapped_round(float(c)) for c in sl[-3:][:self.n_dims]]
                self.nodes.append(node)

                out_line = ""
                for coord in node: out_line += str(coord) + " "
                out_line += "\n" #sl[-1]
                fil.write(out_line)
            for line in lines[start_ele:end_ele]: 
                sl = line.split(' ')
                if self.cell_type=="ED2Q36BS": 
                    cell_nodes = [int(float(node)-1) for node in sl[-4:]]
                    local_nodes = [-1]*36 # With CB-Geo's numbering
                    
                    steps = [self.l0/self.nc0, self.l1/self.nc1]
                    ref_node = self.nodes[cell_nodes[0]]

                    for iloc, nsteps in _ED2Q36BS_nodes.items():
                        current_node = [round(x+n*step, 10) for x,n,step in zip(ref_node, nsteps, steps)]
                        if current_node in self.nodes: local_nodes[iloc] = self.nodes.index(current_node)
                    
                    self.cells.append(local_nodes)

                else: self.cells.append([int(float(node)-1) for node in sl[-self.nn_percell:]])
                
                out_line = ""
                for node in self.cells[-1]: out_line += str(node) + " "
                out_line = out_line[:-1] + "\n"
                fil.write(out_line)

    def __reset_params(self):
        self.params = {"dimensions": self.dimensions, "ncells": self.ncells, "origin": self.origin, "check_duplicates": self.check_duplicates, "cell_type": self.cell_type, "round_decimal": self.round_decimal}

def in_bounds(x,y,cells):
    return (x>=0) and (y>=0) and (x<=cells[0]) and (y<=cells[1])
def id_from_pos(x,y,cells):
    return x + (y*(cells[0]+1))
def id_from_pos_cell(x,y,cells):
    return x + (y*(cells[0]))
gimp_order = [
        [0, 0],
        [1 ,0] ,
        [1 ,1]  ,
        [0, 1] ,
        [-1, 0],
        [2 ,0 ],
        [1 ,-1] ,
        [1 ,2 ] ,
        [2 ,1 ] ,
        [-1, 1] ,
        [0, 2 ],
        [0, -1],
        [-1, -1],
        [2 ,-1 ],
        [2 ,2  ],
        [-1, 2 ]
        ]
def gen_gimp_mesh(length,cells,outdir="consol",mesh=None):
    #print(length)
    #print(cells)
    nd = 2
    length = np.array(length)
    cells = np.array(cells,dtype=int)
    cellvs = cells+1
    resolution = length/(cells)
    verts = np.zeros((cellvs[0]*cellvs[1],2))
    for x in range(cellvs[0]):
        for y in range(cellvs[1]):
            verts[id_from_pos(x,y,cells),0] = x*resolution[0]
            verts[id_from_pos(x,y,cells),1] = y*resolution[1]
    elements = np.zeros(((cells[0])*(cells[1]),16))
    for e_x in range(cells[0]):
        for e_y in range(cells[1]):
            ni = 0
            ei = id_from_pos_cell(e_x,e_y,cells)
            #for dx in range(-1,3):
            #    for dy in range(-1,3):
            for ni in range(16):
                dx,dy = gimp_order[ni]
                #print("{} {}".format(e_x+dx,e_y+dy))
                if in_bounds(e_x+dx,e_y+dy,cells):
                    elements[ei,ni] = id_from_pos(e_x+dx,e_y+dy,cells)
                else:
                    elements[ei,ni] = -1
    #print(verts)
    #print(elements)
    if outdir != None:
        with open("{}/mesh.txt".format(outdir),"w") as f:
            f.write("{}   {}\n".format(len(verts),len(elements)))
            for v in verts:
                f.write("{} {}\n".format(v[0],v[1]))
            for e in elements:
                for ni in e:
                    f.write("{} ".format(int(ni)))
                f.write("\n")
    if mesh:
        mesh.cells = []
        for e in elements:
            mesh.cells.append(e.astype(int).tolist()[0:4])
        mesh.nodes = verts
    return elements
