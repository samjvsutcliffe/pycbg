import gmsh, os, json, pickle, csv, runpy, sys, shutil
import multiprocessing
import numpy as np
import itertools as it  
import __main__ as main

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
    nc1, nc2, nc3 : ints
        Number of cells in each direction (``self.nc0, self.nc1, self.nc2 = self.ncells``).
    origin : tuple of floats
        Origin of the mesh.
    directory : str
        Directory in which the mesh file will be saved.
    check_duplicates : bool
        See CB-Geo documentation.
    cell_type : {'ED3H8', 'ED3H20', 'ED3H64G', 'ED2Q4', 'ED2Q8', 'ED2Q9', 'ED2Q16G'}
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
        elif cell_type=='ED2Q16G': self.nn_percell, self.n_dims = 16, 2
        else : raise ValueError("cell_type is set to '{:}' while it should be one of the following: 'ED3H8', 'ED3H20', 'ED3H64G', 'ED2Q4', 'ED2Q8', 'ED2Q9' or 'ED2Q16G'".format(cell_type))
        
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

        self._dimensions, self._ncells = dimensions, ncells
        self.__reset_params()

        self.write_file()
        self.cells, self.nodes = np.array(self.cells), np.array(self.nodes)

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
        else: raise RuntimeError("Number of dimensions coulnd't be detected, please check `cell_type` is correctly set.")

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
        else: raise RuntimeError("Number of dimensions coulnd't be detected, please check `cell_type` is correctly set.")


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
            gmsh.model.mesh.setOrder(3)

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
                self.cells.append([int(float(node)-1) for node in sl[-self.nn_percell:]])
                out_line = ""
                for node in sl[-self.nn_percell:-1]: out_line += str(int(float(node)-1)) + " "
                out_line += str(int(float(sl[-1])-1)) + "\n"
                fil.write(out_line)

    def __reset_params(self):
        self.params = {"dimensions": self._dimensions, "ncells": self._ncells, "origin": self.origin, "check_duplicates": self.check_duplicates, "cell_type": self.cell_type, "round_decimal": self.round_decimal}

        

class Particles():
    """For defining particles on a :class:`~pycbg.preprocessing.Mesh` and ultimately generate the particles file expected by CB-Geo MPM .

    Parameters
    ----------
    mesh : :class:`~pycbg.preprocessing.Mesh` object
        Mesh in which the particles will be generated.
    npart_perdim_percell : int, optional
        Number of particles along each dimension in a cell. All cells will thus contain ``npart_perdim_percell**3`` equally spaced particles (note that particles are equally spaced within a cell, but not between cells). Default is 1. It is possible to subsequently override the default behavior, manually tuning particles number and locations, see below.
    positions : numpy array, optional
        Particles initial positions, if specified. If `positions=None` (default value), particles are automatically generated using `automatic_generation` and `npart_perdim_percell`. 
    directory : str, optional
        Directory in which the particles file will be saved. If the directory doesn't already exist, it will be created. It is set by default to the current working directory.
    check_duplicates : bool, optional  
        See CB-Geo documentation for informations on this parameter. Default is `True`.
    automatic_generation : {'pycbg', 'cbgeo'}, optional
        Use PyCBG or CB-Geo for automatic material points' generation. CB-Geo will generate the materials points at the Gauss' points of the cells. Default is `'pycbg'`.

    Attributes
    ----------
    positions : numpy array
        Positions of all particles written into the particles file. The id of a particle is the index of its line in this array. Noting `npart` the number of particles, the shape of `positions` is ``(npart,3)``.
    npart_perdim_percell : int
        Number of particles for each dimensions in one cell.
    mesh : :class:`~pycbg.preprocessing.Mesh` object
        Mesh used to generate particles (i.e. the one passed to the constructor).
    directory : str
        Directory in which the particles file will be saved.
    check_duplicates : bool 
        See CB-Geo documentation.
    automatic_generation : str
        Use PyCBG or CB-Geo for automatic material points' generation. CB-Geo will generate the materials points at the Gauss' points of the cells.
    params : dict
        Dictionary containing all necessary parameters (except the `directory` parameter) to create a copy of this `Particles` object. For instance, one can do: `particles_copy = Particles(**existing_particles.params)`.


    Notes
    -----
     - One can manually define the particles (number and positions) by direct modification of the `positions` attribute, after object instantiation.

    Examples
    --------
    Generating automatically 8 particles in a one cell mesh :

    >>> mesh = Mesh((1.,1.,1.), (1,1,1))
    >>> particles = Particles(mesh, 2)
    >>> particles.positions
    array([[0.33333333, 0.33333333, 0.33333333],
           [0.33333333, 0.33333333, 0.66666667],
           [0.33333333, 0.66666667, 0.33333333],
           [0.33333333, 0.66666667, 0.66666667],
           [0.66666667, 0.33333333, 0.33333333],
           [0.66666667, 0.33333333, 0.66666667],
           [0.66666667, 0.66666667, 0.33333333],
           [0.66666667, 0.66666667, 0.66666667]])

    Manually generating four particles :

    >>> mesh = Mesh((1.,1.,1.), (1,1,1))
    >>> particles = Particles(mesh)
    >>> particles.positions = np.array([[.02, .02, .02],
    ...                                 [.05, .02, .02],
    ...                                 [.02, .05, .02],
    ...                                 [.02, .02, .05]])

    Note that a mesh has to be specified even if it isn't used.
    """
    ## TODO: Make the empty initialisation of Particles object possible (without specifying a mesh)
    def __init__(self, mesh, npart_perdim_percell=1, positions=None, directory="", check_duplicates=True, automatic_generation="pycbg"):
        if directory == '' : directory = '/'
        if directory[-1] != '/' : directory += '/'
        if not os.path.isdir(directory): os.mkdir(directory)
        self.directory = directory
        
        self.positions = positions if positions is not None else []
        self.automatic_generation = automatic_generation
        if positions is None: self.create_particles(mesh, npart_perdim_percell, automatic_generation)
        else: self.type = "file"
        
        self.check_duplicates = check_duplicates
        self.n_dims = mesh.n_dims
        self._io_type = "Ascii{:d}D".format(mesh.n_dims) # Shouldn't have to be set to another value
        self._particle_type = "P{:d}D".format(mesh.n_dims) # Shouldn't have to be set to another value

        self.mesh, self.npart_perdim_percell = mesh, npart_perdim_percell

        self.__reset_params()

    def create_particles(self, mesh, npart_perdim_percell=1, automatic_generation="pycbg"):
        """Create the particles using the given mesh.

        Parameters
        ----------
        mesh : :class:`~pycbg.preprocessing.Mesh` object
            Mesh in which the particles will be generated.
        npart_perdim_percell : int, optional
            Number of particles for each dimensions in one cell. All cells will contain ``npart_perdim_percell**3`` equally spaced particles. Note that particles are equally spaced within a cell, not between cells. Default is 1 .
        automatic_generation : {'pycbg', 'cbgeo'}, optional
            Use PyCBG or CB-Geo for automatic material points' generation. CB-Geo will generate the materials points at the Gauss' points of the cells. Default is `'pycbg'`.
        """
        if automatic_generation == "pycbg":
            self.type = "file"
            for ie, e in enumerate(mesh.cells):
                coors = np.array([mesh.nodes[i] for i in e])
                mins, maxs = coors.min(axis=0), coors.max(axis=0)
                step = (maxs-mins)/(npart_perdim_percell)
                poss = [mins + step*(.5+i) for i in range(npart_perdim_percell)]
                coor_ss = []
                for i in range(mesh.n_dims): 
                    coor_ss.append([p[i] for p in poss])
                self.positions += list(it.product(*coor_ss))
            self.positions = np.array(self.positions)
        elif automatic_generation == "cbgeo":
            self.type = "gauss"
            self.cset_id = -1
            self.nparticles_per_dir = npart_perdim_percell
        else: raise ValueError("automatic_generation is set to '{:}' while it should be 'pycbg' or 'cbgeo'".format(automatic_generation))


    def write_file(self, filename="particles"):
        """
        Write the particles file formatted for CB-Geo MPM.
        
        Parameters
        ----------
        filename : str, optional
            Name of the particles file, the extension '.txt' is automatically added. Default is 'particles'.
        """
        self._filename = self.directory + filename + '.txt'
        pfile = open(self._filename, "w") 
        pfile.write("{:d}\n".format(len(self.positions)))   
        for p in self.positions: pfile.write("\t".join(["{:e}"]*self.n_dims).format(*p)+"\n")
        self.__reset_params()

    def __reset_params(self):
        self.params = {"mesh": self.mesh, "npart_perdim_percell": self.npart_perdim_percell, "check_duplicates": self.check_duplicates, "automatic_generation": self.automatic_generation}


class EntitySets():
    """Create and write to a file entity sets for nodes and particles.

    Parameters
    ----------
    mesh : :class:`~pycbg.preprocessing.Mesh` object
        Simulation's mesh. Has to be specified even if only particles sets are defined.
    particles : :class:`~pycbg.preprocessing.Particles` object
        Simulation's particles. Has to be specified even if only nodes sets are defined.
    directory : str, optional
        Directory in which the entity sets file will be saved. If the directory doesn't already exist, it will be created. It is set by default to the current working directory.

    Attributes
    ----------
    nsets : list of lists of ints
        Each element is a list of nodes' ids belonging to the same set. Its index is the id of the node set.
    psets : list of lists of ints
        Each element is a list of particles' ids belonging to the same set. Its index is the id of the particle set.
    mesh : :class:`~pycbg.preprocessing.Mesh` object
        Simulation's mesh.
    particles : :class:`~pycbg.preprocessing.Particles` object
        Simulation's particles.
    params : dict
        Dictionary containing all necessary parameters (except the `directory` parameter) to create a copy of this `EntitySets` object. For instance, one can do: `eset_copy = EntitySets(**existing_eset.params)`.

    
    Notes
    -----
     - The entity sets file is not written upon creating the object. It is thus necessary to run the `write_file` method once all entity sets are created.
     - The user has to define a function for each entity set that indicates which particle or node should be included. See `create_set` method's documention for more informations.

    Examples
    --------
    Creating a node and a particle set in a one cell mesh :

    >>> mesh = Mesh((1.,1.,1.), (1,1,1))
    >>> particles = Particles(mesh, 2)
    >>> entity_sets = EntitySets(mesh, particles)
    >>> node_set_id = entity_sets.create_set(lambda x,y,z: x==0, typ="node")
    >>> particle_set_id = entity_sets.create_set(lambda x,y,z: x<.5, typ="particle")

    Note that this example uses lambda functions to create sets with only one line. One could also use : 

    >>> def x_wall(x, y, z): return x==0
    >>> node_set_id = entity_sets.create_set(x_wall, typ="node")
    """
    def __init__(self, mesh, particles, directory=""):
        self.mesh = mesh
        self.particles = particles
        if directory == '' : directory = '/'
        if directory[-1] != '/' : directory += '/'
        if not os.path.isdir(directory): os.mkdir(directory)
        self.directory = directory

        self.psets, self.nsets = [], []
        self.__reset_params()
    
    def create_set(self, condition_function, typ="particle", kwargs={}):
        """Create a set of nodes or particles and add it to the corresponding
        list. Nodes and particles are selected using `condition_function`.

        Parameters
        ----------
        condition_function : function
            Select particles or nodes using their positions. The inputs should be at least 3 parameters `x`, `y` and `z` that correspond to the position of a node or particle, additional keyword parameters can be passed through `kwargs`. Should return `True` if the node or particle belongs to the set, `False` otherwise.
        typ : {"node", "particle"}, optional
            Type of set to be created. Default is "particle".
        kwargs : dictionary
            Contains the keyword arguments to be passed to `condition_function`
        
        Returns
        -------
        int
            Id of the set just appended.

        Examples
        --------
        Creating a node set using a mesh `mesh` and the particles `particles` previously defined :
        
        >>> mesh = Mesh((1.,1.,1.), (1,1,1))
        >>> particles = Particles(mesh, 2)
        >>> entity_sets = EntitySets(mesh, particles)
        >>> node_set_id = entity_sets.create_set(lambda x,y,z: x==0, typ="node")
        >>> particle_set_id = entity_sets.create_set(lambda x,y,z: x<.5, typ="particle")
        """
        if typ=="particle": points, set_list = self.particles.positions, self.psets
        elif typ=="node": points, set_list = self.mesh.nodes, self.nsets
        else: raise ValueError("`typ` parameter should be 'particle' or 'node'")
        
        ids = []
        for i, p in enumerate(points): 
            if condition_function(*p, **kwargs): ids.append(i)
        set_list.append(ids)
        return len(set_list)-1 if typ=="node" else len(set_list)
    
    def write_file(self, filename="entity_sets"):
        """
        Write the entity sets file formatted for CB-Geo MPM.
        
        Parameters
        ----------
        filename : str, optional
            Name of the entity set file, the extension '.json' is automatically added. Default is 'entity_sets'.
        """
        self._filename = self.directory + filename + '.json'
        main_dic = {}
        for typ, sets_tmp in zip(("particle_sets", "node_sets"), (self.psets, self.nsets)):
            if len(sets_tmp)==0: continue
            sets = []
            for i, current_set in enumerate(sets_tmp): 
                set_id = i if typ=="node_sets" else i+1
                sets.append({"id":set_id, "set":str(current_set)})
            main_dic[typ] = sets
        with open(self._filename, 'w') as fil: json.dump(main_dic, fil, sort_keys=False, indent=4)
        ## Read and rewrite the file, to erase the double quotes around the list in "set" lines (definitely ugly)
        with open(self._filename, 'r') as fil: lines = fil.readlines()
        with open(self._filename, 'w') as fil:
            for line in lines:
                if '"set"' in line: line = line[:18] + line[18:].replace('"', '')
                fil.write(line)
    
    def __reset_params(self):
        self.params = {"mesh": self.mesh, "particles": self.particles}


class Materials():
    """Create materials for particle sets.

    Parameters
    ----------
    n_dims : int, optional
        Number of dimensions in the simulation (2 for 2D, 3 for 3D). Default to 3.

    Attributes
    ----------
    materials : list of dict
        Each element is a dictionnary containing a material's parameters. The index of a material is his id.
    pset_ids : list of (ints or list of ints)
        The element i of this list is the id (or list of ids) of the particle set made of the material defined in ``materials[i]``.
    n_dims : int
        Number of dimensions in the simulation (2 for 2D, 3 for 3D).
    
    Notes
    -----
    Due to (probably) a bug in CB-Geo, materials should be created in the same order than the corresponding particle sets (so particle sets and materials have the same id). 
    
    Examples
    --------
    Creating two materials for two particle sets :

    >>> mesh = Mesh((1.,1.,1.), (1,1,1))
    >>> particles = Particles(mesh, 2)
    >>> entity_sets = EntitySets(mesh, particles)
    >>> lower_particles = entity_sets.create_set(lambda x,y,z: x<.5, typ="particle")
    >>> upper_particles = entity_sets.create_set(lambda x,y,z: x>=.5, typ="particle")
    >>> materials = Materials()
    >>> materials.create_MohrCoulomb(pset_id=lower_particles, density=750,
    ...                                                       youngs_modulus=5.26e7,
    ...                                                       poisson_ratio=.3,
    ...                                                       friction=36.,
    ...                                                       dilation=0.,
    ...                                                       cohesion=1.,
    ...                                                       tension_cutoff=1.,
    ...                                                       softening=False)
    >>> materials.create_Newtonian(pset_id=upper_particles, density=1.225, 
    ...                                                     bulk_modulus=1.42e5, 
    ...                                                     dynamic_viscosity=1.81e3)
    >>> materials.pset_ids
    [0, 1]
    """
    def __init__(self, n_dims=3): 
        self.materials = []
        self.pset_ids = []
        self.n_dims = n_dims

    def create_Newtonian(self, pset_id=0, density=1.225, 
                                          bulk_modulus=1.42e5, 
                                          dynamic_viscosity=1.81e-5):
        """Create a `Newtonian material <https://mpm.cb-geo.com/#/theory/material/newtonian>`_.

        Parameters
        ----------
        pset_id : int or list of ints
            Particle set ids that will be made of this material
        density : float, optional
            Initial density of the material (:math:`kg/m^3`). Default is 1.225 :math:`kg/m^3`.
        bulk_modulus : float, optional
            Bulk modulus of the material (:math:`Pa`). Default is 142 :math:`kPa`.
        dynamic_viscosity : float, otpional
            Dynamic viscosity of the material (:math:`Pa.s`). Default is 18.1 :math:`\mu Pa.s`

        Notes
        -----
        Defaults correspond to air's properties.
        """
        self.pset_ids.append(pset_id)
        self.materials.append({"id": len(self.materials),
                               "type": "Newtonian{:d}D".format(self.n_dims),
                               "density": density,
                               "bulk_modulus": bulk_modulus,
                               "dynamic_viscosity": dynamic_viscosity})
    
    def create_MohrCoulomb(self, pset_id=0, density=1e3,
                                            youngs_modulus=5e7,
                                            poisson_ratio=.3,
                                            friction=36.,
                                            dilation=0.,
                                            cohesion=0.,
                                            tension_cutoff=0.,
                                            softening=False,
                                            peak_pdstrain=0.,
                                            residual_pdstrain=0.,
                                            residual_friction=13.,
                                            residual_dilation=0.,
                                            residual_cohesion=0.):
        """Create a `MohrCoulomb material <https://mpm.cb-geo.com/#/theory/material/mohr-coulomb>`_.

        Parameters
        ----------
        pset_id : int or list of ints
            Particle set id that will be made of this material
        density : float
            Initial density of the material (:math:`kg/m^3`). Default is 1000 :math:`kg/m^3`.
        young_modulus : float
            Young's modulus of the material (:math:`Pa`). Default is 50 :math:`MPa`.
        poisson_ratio : float
            Poisson's ratio of the material. Default is 0.3 .
        friction : float
            Friction angle of the material (:math:`^\circ`). Default is 36 :math:`^\circ`.
        dilation : float
            Dilation angle of the material (:math:`^\circ`). Default is 0 :math:`^\circ`.
        cohesion : float
            Cohesion in the material (:math:`Pa`). Default is 0 :math:`Pa`.
        tension_cutoff : float
            Tension strength of the material (:math:`Pa`). Default is 0 :math:`Pa`.
        softening : bool, optional
            Enable softening option. If `True`, one has to set `peak_pdstrain`, `residual_pdstrain`, `residual_friction`, `residual_dilation` and `residual_cohesion`. Default is `False`.
        peak_pdstrain : float, optional
            Start point of strain softening. Default is 0.
        residual_pdstrain : float, optional
            End point of strain softening. Default is 0.
        residual_friction : float, optional
            Residual friction angle (:math:`^\circ`). Default is 13 :math:`^\circ`.
        residual_dilation : float, optional
            Residual dilation angle (:math:`^\circ`). Default is 0 :math:`^\circ`.
        residual_cohesion : float, optional
            Residual cohesion (:math:`Pa`). Default is 0 :math:`Pa`.
        """
        self.pset_ids.append(pset_id)
        self.materials.append({"id": len(self.materials),
                               "type": "MohrCoulomb{:d}D".format(self.n_dims),
                               "density": density,
                               "youngs_modulus": youngs_modulus,
                               "poisson_ratio": poisson_ratio,
                               "friction": friction,
                               "dilation": dilation,
                               "cohesion": cohesion,
                               "tension_cutoff": tension_cutoff,
                               "softening": softening,
                               "peak_pdstrain": peak_pdstrain,
                               "residual_friction": residual_friction,
                               "residual_dilation": residual_dilation,
                               "residual_cohesion": residual_cohesion,
                               "residual_pdstrain": residual_pdstrain})
    
    def create_NorSand(self, pset_id=0, density=2e3,
                                        poisson_ratio=.17,
                                        reference_pressure=1e5,
                                        friction_cs=27.,
                                        N=.3,
                                        lmbda=.3,
                                        kappa=.08,
                                        gamma=1.,
                                        chi=3.5,
                                        hardening_modulus=2e2,
                                        void_ratio_initial=.38,
                                        p_image_initial=3e6,
                                        bond_model=False,
                                        p_cohesion_initial=1.2e4,
                                        p_dilation_initial=2.4e4,
                                        m_cohesion=1e1,
                                        m_dilation=1.,
                                        m_modulus=1e2,
                                        tolerance=None):
        """Create a `NorSand material <https://mpm.cb-geo.com/#/theory/material/norsand>`_.

        Parameters
        ----------
        pset_id : int or list of ints
            Particle set id that will be made of this material
        density : float
            Initial density of the material (:math:`kg/m^3`). Default is 1000 :math:`kg/m^3`.
        poisson_ratio : float
            Poisson's ratio of the material. Default is 0.3 .
        reference_pressure : float
            Reference pressure (:math:`Pa`). Default is atmospheric pressure, 100 :math:`kPa`.
        friction_cs : float
            Critical state friction angle used to compute M (:math:`^\circ`). Default is 27 :math:`^\circ`.
        N : float
            Volumetric coupling (dilatancy) parameter. Default is 0.3 .
        lmbda : float
            Virgin compression index. Default is 0.3 .
        kappa : bool
            Swell/recompression index. Default to 0.08 .
        chi : float
            Dilatancy coefficient. Default is 3.5 .
        hardening_modulus : float
            (Minimun ?) hardening modulus (:math:`Pa`). Default is 200 :math:`Pa`.
        void_ratio_initial : float
            Initial void ratio. Default is 0.38 .
        p_image_initial : float
            Intersection between the yield surface and the critical state line (:math:`Pa`). Default is 3 :math:`GPa`.
        bond_model : bool, optional
            Enable bond model. Default is False.
        p_cohesion_initial : float, optional
            Initial cohesive pressure for the bond (:math:`Pa`). Default is 120 :math:`kPa`.
        p_dilation_initial : float, optional
            Initial dilation pressure for the bond (:math:`Pa`). Default is 200 :math:`kPa`.
        m_cohesion : float, optional
            Cohesive degradation. Default is 10 .
        m_dilation : float, optional
            Dilative degradation. Default is 1 .
        m_modulus : float, optional
            Bonded modulus effects of the cohesion and dilation. Default is 100 .
        tolerance : float, optional
            Optional tolerance value for computations such as yield condition, set default as machine epsilon .
        """
        self.pset_ids.append(pset_id)
        self.materials.append({"id": len(self.materials),
                               "type": "NorSand{:d}D".format(self.n_dims),
                               "density": density,
                               "reference_pressure": reference_pressure,
                               "poisson_ratio": poisson_ratio,
                               "friction_cs": friction_cs,
                               "N": N,
                               "lambda": lmbda,
                               "kappa": kappa,
                               "gamma": gamma,
                               "chi": chi,
                               "hardening_modulus": hardening_modulus,
                               "void_ratio_initial": void_ratio_initial,
                               "p_image_initial": p_image_initial,
                               "bond_model": bond_model
        })
        if bond_model: self.materials[-1].update({
                               "p_cohesion_initial": p_cohesion_initial,
                               "p_dilation_initial": p_dilation_initial,
                               "m_cohesion": m_cohesion,
                               "m_dilation": m_dilation,
                               "m_modulus": m_modulus,
        })
        if tolerance is not None: self.materials[-1].update({
                               "tolerance": tolerance
        })

    def create_LinearElastic(self, pset_id=0, density=1e3,
                                              youngs_modulus=5e7,
                                              poisson_ratio=.3):
        """Create a `LinearElastic material <https://mpm.cb-geo.com/#/theory/material/linear-elastic>`_.

        Parameters
        ----------
        pset_id : int or list of ints
            Particle set id that will be made of this material
        density : float
            Initial density of the material (:math:`kg/m^3`). Default is 1000 :math:`kg/m^3`.
        young_modulus : float
            Young's modulus of the material (:math:`Pa`). Default is 50 :math:`MPa`.
        poisson_ratio : float
            Poisson's ratio of the material. Default is 0.3 .
        """
        self.pset_ids.append(pset_id)
        self.materials.append({"id": len(self.materials),
                               "type": "LinearElastic{:d}D".format(self.n_dims),
                               "density": density,
                               "youngs_modulus": youngs_modulus,
                               "poisson_ratio": poisson_ratio})

    def create_CustomLaw(self, pset_id=0, density=1e3,
                                          init_state_variables=[],
                                          script_path="custom_law",
                                          function_name="custom_law",
                                          particles_ids=""
                                            ):
        """Create CustomLaw3D material. The behaviour of the material is
        computed using a user-defined python script.

        Parameters
        ----------
        pset_id : int or list of ints
            Particle set id that will be made of this material.
        density : float
            Initial density of the material (:math:`kg/m^3`). Default is 1000 :math:`kg/m^3`.
        init_state_variables : list of floats
            Contains the initial values of the states variables. The order in which they are given is their numbering among states variables : the first one is named "svars_0", the second is named "svars_1", ... Default is an empty list, for no state variables.
        script_path : str
            Path to the user-defined script that compute the material's behaviour. Note that the exentsion `.py` shouldn't be specified. Default is 'custom_law'.
        function_name : str
            Name of the function in `script_path` that compute the stress increment from the strain increment. It should take as input `6 + n_state_vars` arguments. The first 6 are the components of the engineering strain increment (in the directions `xx`, `yy`, `zz`, `xy`, `yz` and `xz` respectively), the others are the state variables. The order of the state variables in the function parameter gives their numbering in the output files (`'svars_0'`, `'svars_1'`, ...).
        """
        self.pset_ids.append(pset_id)
        material_dict = {"id": len(self.materials),
                         "type": "CustomLaw3D",
                         "density": density,
                         "script_path": script_path,
                         "function_name": function_name}
        for i, init_val in enumerate(init_state_variables): material_dict["svars_"+str(i)] = init_val
        self.materials.append(material_dict) 

    def create_PythonModel(self, pset_id=0, density=1e3,
                                            init_state_variables=[],
                                            script_path="custom_law",
                                            function_name="custom_law",
                                            particles_ids=""
                                            ):
        """Create PythonModel3D material. This material is identical to CustomLaw3D except that it is able to run in parallel. For sequential simulations, CustomLaw3D is more efficient than PythonModel3D. 

        Parameters
        ----------
        pset_id : int or list of ints
            Particle set id that will be made of this material.
        density : float
            Initial density of the material (:math:`kg/m^3`). Default is 1000 :math:`kg/m^3`.
        init_state_variables : list of floats
            Contains the initial values of the states variables. The order in which they are given is their numbering among states variables : the first one is named "svars_0", the second is named "svars_1", ... Default is an empty list, for no state variables.
        script_path : str
            Path to the user-defined script that compute the material's behaviour. Note that the exentsion `.py` shouldn't be specified. This script will be copied into pycbg's simulation directory and executed there. Default is 'custom_law'.
        function_name : str
            Name of the function in `script_path` that compute the stress increment from the strain increment. It should take as input `6 + n_state_vars` arguments. The first 6 are the components of the engineering strain increment (in the directions `xx`, `yy`, `zz`, `xy`, `yz` and `xz` respectively), the others are the state variables. The order of the state variables in the function parameter gives their numbering in the output files (`'svars_0'`, `'svars_1'`, ...).
        """
        self.pset_ids.append(pset_id)
        material_dict = {"id": len(self.materials),
                         "type": "PythonModel3D",
                         "density": density,
                         "script_path": script_path,
                         "function_name": function_name}
        for i, init_val in enumerate(init_state_variables): material_dict["svars_"+str(i)] = init_val
        self.materials.append(material_dict) 

    def _set_n_dims(self, n_dims): 
        self.n_dims = n_dims
        for mat in self.materials:
            mat["type"] = mat["type"][:-2] + "{:d}D".format(self.n_dims)


class Simulation():
    """Create a simulation.

    Parameters
    ----------
    title : str, optional
        Simulation title. Default is 'Sim_title'.
    directory : str, optional
        Path to the simulation's directory (will be created if not existent. User-indication of a final '/' is optional). Mesh, particles and entity sets files will be saved in this directory. The result folder is also set to be created by CB-Geo MPM in this directory. Default is `title`.
    input_filename : str, optional
        Name of the input file, the extension `.json` is automatically added. Default is 'input_file'.

    Attributes
    ----------
    mesh : :class:`~pycbg.preprocessing.Mesh` object
        Simulation's mesh. Created using the `create_mesh` method.
    particles : :class:`~pycbg.preprocessing.Particles` object
        Simulation's particles. Created using the `create_particles` method.
    entity_sets : :class:`~pycbg.preprocessing.EntitySets` object
        Simulation's entity sets. Created using the `init_entity_sets` method.
    materials : :class:`~pycbg.preprocessing.Materials` object
        Simulation's materials. Created upon creating the `Simulation` object. 
    init_stresses : numpy array
        Initial stresses for each particle. Noting `npart` the number of particles, its shape is ``(npart, 6)``.
    init_velocities : numpy array
        Initial velocities for each particle. Noting `npart` the number of particles, its shape is ``(npart, 3)``.
    init_volumes : numpy array
        Initial volume for each particle. Noting `npart` the number of particles, its shape is ``(npart, 1)``.
    input_filename : str
        Path to the CB-Geo MPM json input file to create, the extension '.json' is automatically added.. If `directory='.'`, the title of the simulation is automatically added before the user-specified filename. Default is 'input_file' in `directory`.
    gravity : list of floats
        Gravity vector for the simulation.
    title : str
        Simulation title.
    directory : str
        Path to the simulation's directory.
    custom_params : dict
        Dictionary containing user-defined parameters. It will be saved in the :class:`~pycbg.preprocessing.Simulation` object when the input file is written. Its element should be appended using the `add_custom_parameters` method.
    analysis_params : dict
        Dictionary containing all analysis parameters. It is possible to pass one simulation's analysis parameters to another, for instance: `sim2.set_analysis_parameters(**sim1.analysis_params)`.

    Notes
    -----
    `Mesh`, `Particles`, `EntitySets` and :class:`~pycbg.preprocessing.Materials` objects are created from the `Simulation` object.

    Examples
    --------
    Simulating under gravity a column made of two materials :
    
    >>> sim = Simulation()
    >>> sim.create_mesh(dimensions=(1.,1.,10.), ncells=(1,1,10))
    >>> sim.create_particles(npart_perdim_percell=1)
    >>> sim.init_entity_sets()
    >>> lower_particles = sim.entity_sets.create_set(lambda x,y,z: z<10, typ="particle")
    >>> upper_particles = sim.entity_sets.create_set(lambda x,y,z: z>=10, typ="particle")
    >>> sim.materials.create_MohrCoulomb(pset_id=lower_particles)
    >>> sim.materials.create_Newtonian(pset_id=upper_particles)
    >>> walls = []
    >>> walls.append([sim.entity_sets.create_set(lambda x,y,z: x==lim, typ="node") for lim in [0, sim.mesh.l0]])
    >>> walls.append([sim.entity_sets.create_set(lambda x,y,z: y==lim, typ="node") for lim in [0, sim.mesh.l1]])
    >>> walls.append([sim.entity_sets.create_set(lambda x,y,z: z==lim, typ="node") for lim in [0, sim.mesh.l2]])
    >>> for direction, sets in enumerate(walls): _ = [sim.add_velocity_condition(direction, 0., es) for es in sets]
    >>> sim.set_gravity([0,0,-9.81])
    >>> sim.set_analysis_parameters(dt=1e-3, nsteps=1.5e8, output_step_interval=7.5e6)
    >>> sim.write_input_file()
    """
    def __init__(self, title='Sim_title', directory='', input_filename='input_file'):
        if directory == '': directory = title
        if directory[-1] != '/': directory += '/'
        if not os.path.isdir(directory): os.mkdir(directory)
        self.directory = directory

        self.input_filename = input_filename + ".json"
        self.title = title
        self.materials = Materials()
        
        self.set_analysis_parameters()
        self.math_functions = []
        self.entity_sets = None
        self.init_stresses = None
        self.init_velocities = None
        self.init_volumes = None

        self.custom_params = {}

        self.__boundary_conditions = {"velocity_constraints": [],
                                      "friction_constraints": [],
                                      "particles_velocity_constraints": []}
        self.gravity = [0,0,0]
        self.__nodal_forces = []
        self.__particle_traction = []
        self.__init_stress_filename = self.directory + "particles_stresses.txt"
        self.__init_velocity_filename = self.directory + "particles_velocities.txt"
        self.__init_volumes_filename = self.directory + "particles_volumes.txt"

    def create_mesh(self, *args, **kwargs):
        """Defines the simulation's mesh, with the generation of an appropriate mesh file for CB-Geo MPM.

        Parameters
        ----------
        dimensions : tuple of floats
            Dimensions of the mesh. Its length should be 3, with `dimensions[n]` the dimension of the mesh on the axis `n`.
        ncells : tuple of ints
            Number of cells in each direction. Its length should be 3, with `ncells[n]` the number of cells on the axis `n`.
        check_duplicates : bool, optional
            See CB-Geo documentation for informations on this parameter. Default is `True`.
        cell_type : {'ED3H8', 'ED3H20', 'ED3H64G'}, optional
            Type of cell. Only 3D Hexahedrons are supported. The number of nodes can be 8, 20 or 64. Default is 'ED3H8'.
        """
        if "directory" in kwargs or len(args)>2: raise TypeError("`directory` parameter is defined by the `Simulation` object")
        self.mesh = Mesh(*args, directory=self.directory, **kwargs)
        self.mesh.write_file(filename="mesh")
        self.materials._set_n_dims(self.mesh.n_dims)

    def create_particles(self, *args, **kwargs):
        """Defines the simulation's particles, with the generation of an appropriate particles file for CB-Geo MPM.

        Parameters
        ----------
        npart_perdim_percell : int, optional
            Number of particles for each dimension in one cell. All cells will contain ``npart_perdim_percell**3`` equally spaced particles. Note that particles are equally spaced within a cell, not between cells. Default is 1 .
        check_duplicates : bool, optional
            See CB-Geo documentation for informations on this parameter. Default is `True`.
        """
        if "mesh" in kwargs: raise TypeError("`mesh` parameter is defined by the `Simulation` object")
        if "directory" in kwargs or len(args)>1: raise TypeError("`directory` parameter is defined by the `Simulation` object")
        self.particles = Particles(mesh=self.mesh, *args, directory=self.directory, **kwargs)
        self.particles.write_file(filename="particles")

    def init_entity_sets(self): 
        """Create the simulation's :class:`~pycbg.preprocessing.EntitySets`
        object.

        Has to be called after mesh and particles creation.
        """
        self.entity_sets = EntitySets(mesh=self.mesh, particles=self.particles, directory=self.directory)

    def add_velocity_condition(self, dir, vel_value, entity_set, typ="node", math_function_id=None):
        """Add a velocity condition on a node or particle set.

        Parameters
        ----------
        dir : {0, 1, 2}
            Axis on which the velocity is imposed.
        vel_value : float
            Imposed velocity's value (:math:`m.s^{-1}`).
        entity_set : int
            Id of the entity set on which the velocity is imposed.
        typ : {"node", "particle"}, optional
            Type of set on which the velocity is imposed. Default is "particle".
        math_function_id : int, optional
            Id of the math function to use. Default value is `None` (the velocity is then constant).
        """
        if typ=="particle": list_name, key_name = "particles_velocity_constraints", "pset_id"
        elif typ=="node": list_name, key_name = "velocity_constraints", "nset_id"
        else: raise ValueError("`typ` parameter should be 'particle' or 'node'")

        self.__boundary_conditions[list_name].append({key_name: entity_set,
                                                      "dir": dir,
                                                      "velocity": vel_value})
        if math_function_id != None: self.__boundary_conditions[list_name][-1]["math_function_id"] = math_function_id
        
    def add_friction_condition(self, dir, sgn_n, frict_value, node_set):
        """Add a friction condition on a node set.

        Parameters
        ----------
        dir : {0, 1, 2}
            Axis of the normal vector to the plane where friction is acting.
        sgn_n : {-1, 1}
            Sign of the normal vector to the plane where friction is acting.
        vel_value : float
            Imposed friction coefficient's value.
        entity_set : int
            Id of the entity set on which the velocity is imposed.
        typ : {"node", "particle"}, optional
            Type of set on which the velocity is imposed. Default is "particle".
        """
        self.__boundary_conditions["friction_constraints"].append({"nset_id": node_set,
                                                                   "dir": dir,
                                                                   "sign_n": sgn_n,
                                                                   "friction": frict_value})
    
    def add_math_function(self, times, values):
        """Add a math function to the simulation. The function can only be
        piecewise-linear.

        Parameters
        ----------
        times : list of floats
            Contains the times at which the values of the math function are given. The first element should always be `0.` and the last should always be `nsteps*dt`.
        values : list of floats
            Contains the values of the math function for each time given in `times`.

        Returns
        -------
        int
            Id of the math function just appended.
        """
        fct_id = len(self.math_functions)
        self.math_functions.append({"id":fct_id, "type": 'Linear', "xvalues": str(times), "fxvalues": str(values)})

        return fct_id

    def add_force(self, dir, force, entity_set, typ="node", math_function_id=None):
        """Add a force on all the element in a entity set.

        Parameters
        ----------
        dir : {0, 1, 2}
            Axis on which the force is imposed.
        force : float
            Imposed force's value (:math:`N`).
        entity_set : int
            Id of the entity set on which the velocity is imposed.
        typ : {"node", "particle"}, optional
            Type of set on which the velocity is imposed. Default is "particle".
        math_function_id : int, optional
            Id of the math function to use. Default value is `None` (the load is then static).
        """
        if typ=="node":
            self.__nodal_forces.append({"nset_id": entity_set,
                                        "dir": dir,
                                        "force": force})
            if math_function_id != None: self.__nodal_forces[-1]["math_function_id"] = math_function_id
        elif typ=="particle":
            self.__particle_traction.append({"pset_id": entity_set,
                                             "dir": dir,
                                             "traction": force})
            if math_function_id != None: self.__particle_traction[-1]["math_function_id"] = math_function_id

    
    def set_initial_particles_stresses(self, init_stresses):
        """Set the initial stresses for each particle.

        Parameters
        ----------
        init_stresses : numpy array
            Initial stresses for each particle. Noting `npart` the number of particles, it should have the shape ``(npart, 6)``.
        """
        self.init_stresses = init_stresses

        psfile = open(self.__init_stress_filename, "w") 
        psfile.write("{:d}\n".format(len(self.particles.positions)))   
        for ps in init_stresses: psfile.write("{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\n".format(*ps)) 

    def set_initial_particles_velocities(self, init_velocities):
        """Set the initial velocities for each particle.

        Parameters
        ----------
        init_velocities : numpy array
            Initial velocities for each particle. Noting `npart` the number of particles, it should have the shape ``(npart, 3)``.
        """
        self.init_velocities = init_velocities

        psfile = open(self.__init_velocity_filename, "w") 
        psfile.write("{:d}\n".format(len(self.particles.positions)))   
        for ps in init_velocities: psfile.write("\t".join(["{:e}"]*self.mesh.n_dims).format(*ps)+"\n") 
    
    def set_initial_particles_volumes(self, init_volumes):
        """Set the initial volume for each particle.

        Parameters
        ----------
        init_volumes : numpy array
            Initial volumes for each particle. Noting `npart` the number of particles, it should have the shape ``(npart, 1)``.
        """
        self.init_volumes = init_volumes

        psfile = open(self.__init_volumes_filename, "w") 
        for i, ps in enumerate(init_volumes): psfile.write("{:d}\t{:e}\n".format(i, ps)) 

    
    def set_gravity(self, gravity): 
        """Set the value of gravity. If this method isn't called, gravity is
        `[0,0,0]`.

        Parameters
        ----------
        gravity : list of floats
            Gravity's value on each axis (:math:`m/s^2`).
        """
        self.gravity = gravity

    def set_analysis_parameters(self, type="MPMExplicit3D", mpm_scheme="usf", damping=0.05, locate_particles=False, dt=1e-05, velocity_update='flip', nsteps=2000, verbosity=100, output_step_interval=100):
        """Set the analysis parameters. Has to be called before
        `write_input_file`.

        Parameters
        ----------
        type : {'MPMExplicit2D', 'MPMExplicit3D'}, optional
            Analysis type. Default is 'MPMExplicit3D'.
        mpm_scheme : {'usf', 'usl', 'musl'}, optional
            MPM scheme for the stress update. The scheme can be "Update Stress First" ('usf'), "Update Stress Last" ('usl') or "Modified Update Stress Last" ('musl').
        damping : float, optional
            Cundall's damping. Should verify : ``0 <= damping < 1``. Default is 0.05 .
        locate_particles : bool, optional
            Stops the simulation when particles go outside the mesh if `True`. Default is `False`.
        dt : float, optional
            Time step (:math:`s`). Default is `1e-5` :math:`s`.
        velocity_update : {'pic', 'flip', 'flipX', 'apic', 'nflip', 'nflipX'}, optional
            How to compute particle's velocity. If `'pic'` nodal velocity is directly interpolated to particles. If `'flip'` nodal velocity is computed from the acceleration interpolated to particles. If `'flipX'`, the nodal velocity will be a proportion between the `'flip'` velocity and the `'pic'` velocity: `X*velocity('flip') + (1-X)*velocity('pic')`. If `'apic'`, momentum is interpolated according to *Jiang, Chenfanfu, et al., 2015*. Default is `'flip'`.
            If `'nflip'`, the velocity will be computed using `'flip'` but the material points will be moved according to this velocity (not the nodal one). The `'nflipX'` value combines `'flipX'` and `'nflip'`.
        nsteps : int, optional
            Number of steps to be performed. Default is 2000.
        verbosity : int, optional
            Number of info lines (with current step number) to be printed in the console during MPM execution. Default is 100.
        output_step_interval : int, optional
            Number of steps between two data points. Default is 100.
        """
        try: detected_type = 'MPMExplicit{:d}D'.format(self.mesh.n_dims)
        except: detected_type = type
        self.analysis_params = {"type": detected_type,
                           "mpm_scheme": mpm_scheme,
                           "locate_particles": locate_particles,
                           "dt": dt,
                           "damping": {"type": "Cundall", "damping_factor": damping},
                           "velocity_update": velocity_update,
                           "nsteps": int(nsteps),
                           "verbosity": int(verbosity)}
        self.post_processing = {"path": self.directory + "results/",
                                "output_steps": int(output_step_interval)}
    
    def write_input_file(self):
        """Write the input file."""
        mesh_dic = {"mesh": self.mesh._filename,
                    "boundary_conditions": self.__boundary_conditions,
                    "isoparametric": self.mesh._isoparametric,
                    "check_duplicates": self.mesh.check_duplicates,
                    "cell_type": self.mesh.cell_type,
                    "io_type": self.mesh._io_type,
                    "node_type": self.mesh._node_type}
        if self.entity_sets is not None: 
            self.entity_sets.write_file(filename="entity_sets")
            mesh_dic["entity_sets"] = self.entity_sets._filename
        if self.init_stresses is not None: 
            mesh_dic["particles_stresses"] = self.__init_stress_filename
        if self.init_velocities is not None: 
            mesh_dic["particles_velocities"] = self.__init_velocity_filename
        if self.init_volumes is not None: 
            mesh_dic["particles_volumes"] = self.__init_volumes_filename
        if self.particles.type == "file":
            particles_list = [{"generator": {"check_duplicates": self.particles.check_duplicates,
                                            "location": self.particles._filename,
                                            "io_type": self.particles._io_type,
                                            "pset_id": 0,
                                            "particle_type": self.particles._particle_type,
                                            "material_id": 0,
                                            "type": self.particles.type}}]
        elif self.particles.type == "gauss":
            particles_list = [{"generator": {"check_duplicates": self.particles.check_duplicates,
                                            "pset_id": 0,
                                            "cset_id": self.particles.cset_id, 
                                            "nparticles_per_dir": self.particles.nparticles_per_dir,
                                            "particle_type": self.particles._particle_type,
                                            "material_id": 0,
                                            "type": self.particles.type}}]

        material_sets_list = [] 
        for m_id, ps_id in enumerate(self.materials.pset_ids):
            if type(ps_id)==list: 
                for ps_id_p in ps_id: material_sets_list.append({"material_id": m_id, "pset_id": ps_id_p})
            else: material_sets_list.append({"material_id": m_id, "pset_id": ps_id})

        if self.gravity==[0,0,0] and self.mesh.n_dims==2: self.gravity = [0,0]
        external_loading_conditions_dic = {"gravity": self.gravity}
        if len(self.__nodal_forces) != 0: external_loading_conditions_dic["concentrated_nodal_forces"] = self.__nodal_forces
        if len(self.__particle_traction) != 0: external_loading_conditions_dic["particle_surface_traction"] = self.__particle_traction
    
        dic = {"title": self.title,
               "mesh": mesh_dic,
               "particles": particles_list,
               "materials": self.materials.materials,
               "material_sets": material_sets_list,
               "external_loading_conditions": external_loading_conditions_dic,
               "analysis": dict(self.analysis_params, **{"uuid": self.title}),
               "post_processing": self.post_processing}
        
        if len(self.math_functions) != 0: dic["math_functions"] = self.math_functions

        with open(self.directory + self.input_filename, 'w') as fil: json.dump(dic, fil, sort_keys=False, indent=4)

        ## Read and rewrite the file, to erase the double quotes around the list in "xvalues" and "fxvalues" lines (definitely ugly)
        with open(self.directory + self.input_filename, 'r') as fil: lines = fil.readlines()
        with open(self.directory + self.input_filename, 'w') as fil:
            for line in lines:
                if '"xvalues"' in line: line = line[:22] + line[22:].replace('"', '')
                if '"fxvalues"' in line: line = line[:23] + line[23:].replace('"', '')
                fil.write(line)

        save_name = self.directory + self.title + ".Simulation"
        with open(save_name, 'wb') as fil : pickle.dump(self, fil)

        ## If PythonModel3D or CustomLaw3D is used, copy the rve script into the simulation's directory
        mat_types = [mat["type"] for mat in self.materials.materials]
        if "PythonModel3D" in mat_types or "CustomLaw3D" in mat_types:
            try: mat = self.materials.materials[mat_types.index("PythonModel3D")]
            except ValueError: mat = mat_types[mat_types.index("CustomLaw3D")]

            shutil.copy(mat["script_path"]+".py", self.directory)
    
    def add_custom_parameters(self, dic):
        """Add `dict` content in `custom_params`.

        Parameters
        ----------
        dic : dict
            Dictionary containing the parameters to be appended in `custom_params`.
        """
        for key, val in dic.items(): self.custom_params[key] = val


def setup_batch(script_path, params, directory='', cbgeo_executable=None, ncores="max", ncores_perjob="max"):
    """Setup a serie of simulations. The simulations are based on the script
    which path is `script_path` but additional variables are defined, each
    simulation can have different values for these variables (defined in
    `params`). One should be careful that those additional parameters do not
    get overwritten in the script. Two additional variables will always be
    insert into the script :

     - `sim_dir`, the `str` that contains the path of the simulation's directory
     - `sim_title`, the `str` that contains the title of the simulation
    These two parameters should be used by the base script.

    Parameters
    ----------
    script_path : str
        Path of the base script.
    params : list of dict, dict
        Batch's parameters, it contains the different values for the additional parameters. If it is a list of dict each element should be the dictionary containing a parameter set, the keys being the names of the variables. If it is a dict all keys should be the variables' names and their values should be lists, a list of parameter sets containing all the combinations between the parameters will then be used.
    directory : str, optional
        Path to the batch's directory. The directory of each simulation will be inside. Default is `''`.
    cbgeo_executable : str, optional
        Path to a cbgeo executable. It will be used to generate a bash script that launches the batch. Default is `None`.
    ncores : int or "max", optional
        Number of core to use for the batch. It will only affect the bash script that launches the batch. Default is `'max'`.
    ncores_perjob : int or "max", optional
        Number of core to use with each simulation. It will only affect the bash script that launches the batch. Default is `'max'`.

    Notes
    -----
     - If `ncores` and `ncores_perjob` ar set to `max`, all CPUs will be used for each simulation. Simulations will then be executed one at a time.
     - Be careful if you don't execute the pycbg script on the same machine used to run the simulation: the maximum number of cpu will be the one of the first machine.
    """
    if directory == '' : directory = '/'
    if directory[-1] != '/' : directory += '/'
    if not os.path.isdir(directory): os.mkdir(directory)

    set_executable = cbgeo_executable is not None

    if type(params) == list : param_sets = params
    elif type(params) == dict :
        all_combinations = it.product(*params.values())
        param_sets = [{key:val for key, val in zip(params.keys(), val_set)} for val_set in all_combinations]

    with open(script_path, 'r') as fil: script = fil.readlines()

    script_a, script_b, insert_flag = [], [], True
    for line in script:
        if "BATCH PARAMETERS INSERTION" in line:
            insert_flag = False
        else:
            if insert_flag: script_a.append(line)
            else: script_b.append(line)
    if insert_flag: script_b, script_a = script_a, []
    
    table_file = open(directory + "parameters_sets.table", "w")
    header = "sim_id"
    if type(params) == list : 
        for key in params[0].keys(): header += "\t" + key
    elif type(params) == dict :
        for key in params.keys(): header += "\t" + key
    table_file.write(header + "\n")

    if set_executable: 
        if ncores == 'max': ncores = multiprocessing.cpu_count()
        if ncores_perjob == 'max': ncores_perjob = multiprocessing.cpu_count()
        if ncores_perjob > ncores: ncores = ncores_perjob
        cores_str = " -p {:d}".format(ncores_perjob)
        n_simultaneous_sim = ncores//ncores_perjob

        batch_launcher_file = open(directory + "start_batch.sh", "w")
        batch_launcher_file.write("#!/bin/bash\n\n")
        batch_launcher_file.write("""n_sim_slot_available="{:d}"\necho "n_sim_slot_available=$n_sim_slot_available" > /tmp/n_sim_slot_available.sh\n\n""".format(n_simultaneous_sim))

    
    for sim_id, param_set in enumerate(param_sets):
        sim_title = "sim{:d}".format(sim_id)
        sim_dir = directory + sim_title + "/"
        if not os.path.isdir(sim_dir): os.mkdir(sim_dir)

        affectation_lines  = ["\n# Batch parameters:\n"]
        affectation_lines += [key + " = " + str(val) + "\n" for key, val in param_set.items()]
        affectation_lines += ["sim_title = '{:}'\n".format(sim_title)]
        affectation_lines += ["sim_dir = '{:}'\n\n".format(sim_dir)]
        affectation_lines += ["# Base script:\n"]

        out_script = script_a + affectation_lines + script_b

        out_script_path = sim_dir + "pycbg_script.py"
        with open(out_script_path, "w") as fil: 
            for line in out_script: fil.write(line)
        runpy.run_path(out_script_path)
        
        param_line = str(sim_id)
        for val in param_set.values(): param_line += "\t" + str(val)
        table_file.write(param_line + "\n")

        if set_executable: 
            batch_launcher_file.write("""until [ "$n_sim_slot_available" -gt "0" ] \ndo\n\tsleep 0.1\n\t. /tmp/n_sim_slot_available.sh\ndone\n""")
            batch_launcher_file.write("""n_sim_slot_available=$((--n_sim_slot_available)); echo "n_sim_slot_available=$n_sim_slot_available" > /tmp/n_sim_slot_available.sh \n""")
            batch_launcher_file.write("""( {:}{:} -f "$(pwd)/" -i {:}input_file.json >> {:}cbgeo.log; . /tmp/n_sim_slot_available.sh; n_sim_slot_available=$((++n_sim_slot_available)); echo "n_sim_slot_available=$n_sim_slot_available" > /tmp/n_sim_slot_available.sh ) &\n""".format(cbgeo_executable, cores_str, sim_dir, sim_dir))

    table_file.close()
    if set_executable: batch_launcher_file.close()
    sys.exit()


    


