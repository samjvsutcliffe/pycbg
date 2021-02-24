import gmsh, os, json, pickle, csv, runpy, sys
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
    directory : str, optional
        Directory in which the mesh file will be saved. If the directory doesn't already exist, it will be created. It is set by default to the current working directory.
    check_duplicates : bool, optional  
        See CB-Geo documentation for informations on this parameter. Default is `True`.
    cell_type : {'ED3H8', 'ED3H20', 'ED3H64'}, optional
        Type of cell. Only 3D Hexahedrons are supported. The number of nodes can be 8, 20 or 64. Default is 'ED3H8'.

    Attributes
    ----------
    nodes : numpy array
        Positions of all nodes in the mesh. The id of a node is the index of its line in this array. Noting `nnodes` the number of nodes, the shape of `nodes` is ``(nnodes,3)``.
    cells : numpy array
        Connections between cells and nodes. Each line corresponds to a cell, its index is the cell's id. The columns correspond to the ids of the nodes composing a cell. Noting `nnode_pcell` the number of nodes per cell (8, 20 or 64), the shape of `cells` is ``(ncells,nnode_pcell)``.
    filename : str
        Path to the mesh file from the current working directory.
    dimensions : tuple of floats
        Dimensions of the mesh.
    l0, l1, l2 : floats
        Dimensions of the mesh (``self.l0, self.l1, self.l2 = self.dimensions``).
    ncells : tuple of ints
        Number of cells in each direction.
    nc1, nc2, nc3 : ints
        Number of cells in each direction (``self.nc0, self.nc1, self.nc2 = self.ncells``).
    directory : str
        Directory in which the mesh file will be saved.
    check_duplicates : bool
        See CB-Geo documentation.
    cell_type : {'ED3H8', 'ED3H20', 'ED3H64'}
        Type of cell. 

    Notes
    -----
     - The mesh file is written upon creating the object.
     - The nodes' coordinates are always positive.
     - The point `(0,0,0)` will always be a node of the mesh. 
     - The maximum coordinates are ``self.l0``, ``self.l1`` and ``self.l2`` on the axis 0, 1 and 2 respectively.

    Examples
    --------
    Creating a cubic mesh of 1000 cells : 

    >>> mesh = Mesh((1.,1.,1.), (10,10,10))
    >>> mesh.nc0 * mesh.nc1 * mesh.nc2
    1000
    """
    ## TODO: - Test 'ED3H20' and 'ED3H64'
    ##       - Avoid having to write the mesh file from gmsh for rewritting it again
    ##       - Make crete_mesh usable by the user 

    def __init__(self, dimensions, ncells, directory="", check_duplicates=True, cell_type="ED3H8"):
        self.set_parameters(dimensions, ncells)
        if not os.path.isdir(directory) and directory!='' : os.mkdir(directory)
        self.filename = directory + "mesh.msh"
        
        self.write_file()
        self.cells, self.nodes = np.array(self.cells), np.array(self.nodes)

        self.check_duplicates = check_duplicates
        self.cell_type = cell_type
        self._isoparametric = False # Shouldn't have to be set to another value
        self._io_type = "Ascii3D" # Shouldn't have to be set to another value
        self._node_type = "N3D" # Shouldn't have to be set to another value
    def set_parameters(self, dimensions, ncells):
        """Set the dimensions and number of cells of the mesh.

        Parameters
        ----------
        dimensions : tuple of floats
            Dimensions of the mesh. Its length should be 3, with `dimensions[n]` the dimension of the mesh on the axis `n`.
        ncells : tuple of ints
            Number of cells in each direction. Its length should be 3, with `ncells[n]` the number of cells on the axis `n`.
        """
        self.l0, self.l1, self.l2 = dimensions
        self.nc0, self.nc1, self.nc2 = ncells

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
        
        p = gmsh.model.geo.addPoint(0, 0, 0)
        l = gmsh.model.geo.extrude([(0, p)], self.l0, 0, 0, [self.nc0], [1])
        s = gmsh.model.geo.extrude([l[1]], 0, self.l1, 0, [self.nc1], [1], recombine=True)
        v = gmsh.model.geo.extrude([s[1]], 0, 0, self.l2, [self.nc2], [1], recombine=True)
        gmsh.model.geo.synchronize()
        gmsh.model.addPhysicalGroup(3, [v[1][1]])
        gmsh.model.mesh.generate(3)

    def write_file(self):
        """Write the mesh file formated for CB-Geo."""
        self.create_mesh()

        gmsh.write(self.filename)
        gmsh.finalize()

        self.__reformat_from_gmsh()
    
    def __reformat_from_gmsh(self): # Not meant for the user
        """Reads the mesh file generated by gmsh and reformats it for CB-
        Geo."""
        with open(self.filename, 'r') as fil: lines = fil.readlines()
        for i, line in enumerate(lines):
            if "$Nodes" in line : nn, start_nodes = int(float(lines[i+1])), i+2
            elif "$EndNodes" in line : end_nodes = i
            elif "$Elements" in line : ne, start_ele = int(float(lines[i+1])), i+2
            elif "$EndElements" in line : end_ele = i
        
        self.cells, self.nodes = [], []
        with open(self.filename, 'w') as fil:
            fil.write("{:d}\t{:d}\n".format(nn, ne))
            for line in lines[start_nodes:end_nodes]: 
                sl = line.split(' ')
                self.nodes.append([float(c) for c in sl[-3:]])
                out_line = ""
                for node in sl[-3:-1]: out_line += node + " "
                out_line += sl[-1]
                fil.write(out_line)
            for line in lines[start_ele:end_ele]: 
                sl = line.split(' ')
                self.cells.append([int(float(node)-1) for node in sl[-8:]])
                out_line = ""
                for node in sl[-8:-1]: out_line += str(int(float(node)-1)) + " "
                out_line += str(int(float(sl[-1])-1)) + "\n"
                fil.write(out_line)
        

class Particles():
    """Create and write to a file particles from a :class:`~pycbg.preprocessing.Mesh` object.

    Parameters
    ----------
    mesh : :class:`~pycbg.preprocessing.Mesh` object
        Mesh in which the particles will be generated.
    npart_perdim_percell : int, optional
        Number of particles for each dimensions in one cell. All cells will contain ``npart_perdim_percell**3`` equally spaced particles. Note that particles are equally spaced within a cell, not between cells. Default is 1 .
    directory : str, optional
        Directory in which the particles file will be saved. If the directory doesn't already exist, it will be created. It is set by default to the current working directory.
    check_duplicates : bool, optional  
        See CB-Geo documentation for informations on this parameter. Default is `True`.

    Attributes
    ----------
    particles : numpy array
        Positions of all particles written into the particles file. The id of a particle is the index of its line in this array. Noting `npart` the number of particles, the shape of `particles` is ``(npart,3)``.
    filename : str
        Path to the particles file from the current working directory.
    npart_perdim_percell : int
        Number of particles for each dimensions in one cell.
    directory : str
        Directory in which the particles file will be saved.
    check_duplicates : bool 
        See CB-Geo documentation.

    Notes
    -----
     - The particles file is written upon creating the object.
     - One can manually generate the particles by directly setting the `particles` attribute. It is then necessary to rewrite the particles file using the `write_file` method.

    Examples
    --------
    Generating 8 particles in a one cell mesh :

    >>> mesh = Mesh((1.,1.,1.), (1,1,1))
    >>> particles = Particles(mesh, 2)
    >>> particles.particles
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
    >>> particles.particles = np.array([[.02, .02, .02],
    ...                                 [.05, .02, .02],
    ...                                 [.02, .05, .02],
    ...                                 [.02, .02, .05]])
    >>> particles.write_file()

    Note that a mesh has to be specified even if it isn't used.
    """
    ## TODO: Make the empty initialisation of Particles object possible (without specifying a mesh)

    def __init__(self, mesh, npart_perdim_percell=1, directory="", check_duplicates=True):
        if not os.path.isdir(directory) and directory!='' : os.mkdir(directory)
        self.filename = directory + "particles.txt"
        self.particles = []
        self.create_particles(mesh, npart_perdim_percell)
        
        self.check_duplicates = check_duplicates
        self._io_type = "Ascii3D" # Shouldn't have to be set to another value
        self._particle_type = "P3D" # Shouldn't have to be set to another value
        self._type = "file" # Shouldn't have to be set to another value

    def create_particles(self, mesh, npart_perdim_percell=1):
        """Create the particles using the given mesh.

        Parameters
        ----------
        mesh : :class:`~pycbg.preprocessing.Mesh` object
            Mesh in which the particles will be generated.
        npart_perdim_percell : int, optional
            Number of particles for each dimensions in one cell. All cells will contain ``npart_perdim_percell**3`` equally spaced particles. Note that particles are equally spaced within a cell, not between cells. Default is 1 .
        """
        for ie, e in enumerate(mesh.cells):
            coors = np.array([mesh.nodes[i] for i in e])
            mins, maxs = coors.min(axis=0), coors.max(axis=0)
            steps = (maxs-mins)/(npart_perdim_percell+1)
            poss = [mins + steps*i for i in range(1, npart_perdim_percell+1)]
            xs, ys, zs = [p[0] for p in poss], [p[1] for p in poss], [p[2] for p in poss]
            for x in xs:
                for y in ys:
                    for z in zs:
                        self.particles.append([x, y, z])
        self.particles = np.array(self.particles)

    def write_file(self):
        """Write the particles file formatted for CB-Geo."""
        pfile = open(self.filename, "w") 
        pfile.write("{:d}\n".format(len(self.particles)))   
        for p in self.particles: pfile.write("{:e}\t{:e}\t{:e}\n".format(*p)) 

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
    filename : str
        Path to the entity sets file from the current working directory.
    mesh : :class:`~pycbg.preprocessing.Mesh` object
        Simulation's mesh.
    particles : :class:`~pycbg.preprocessing.Particles` object
        Simulation's particles.
    
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
        if not os.path.isdir(directory) and directory!='' : os.mkdir(directory)
        self.filename = directory + "entity_sets.txt"

        self.psets, self.nsets = [], []
    
    def create_set(self, condition_function, typ="particle"):
        """Create a set of nodes or particles and add it to the corresponding
        list. Nodes and particles are selected using `condition_function`.

        Parameters
        ----------
        condition_function : function
            Select particles or nodes using their positions. The inputs should be 3 parameters `x`, `y` and `z` that correspond to the position of a node or particle. Should return `True` if the node or particle belongs to the set, `False` otherwise.
        typ : {"node", "particle"}, optional
            Type of set to be created. Default is "particle".
        
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
        if typ=="particle": points, set_list = self.particles.particles, self.psets
        elif typ=="node": points, set_list = self.mesh.nodes, self.nsets
        else: raise ValueError("`typ` parameter should be 'particle' or 'node'")
        
        ids = []
        for i, p in enumerate(points): 
            if condition_function(*p): ids.append(i)
        set_list.append(ids)
        return len(set_list)-1
    
    def write_file(self):
        """Write the entity sets file formatted for CB-Geo."""
        main_dic = {}
        for typ, sets_tmp in zip(("particle_sets", "node_sets"), (self.psets, self.nsets)):
            if len(sets_tmp)==0: continue
            sets = []
            for i, current_set in enumerate(sets_tmp): sets.append({"id":i, "set":str(current_set)})
            main_dic[typ] = sets
        with open(self.filename, 'w') as fil: json.dump(main_dic, fil, sort_keys=False, indent=4)
        ## Read and rewrite the file, to erase the double quotes around the list in "set" lines (definitely ugly)
        with open(self.filename, 'r') as fil: lines = fil.readlines()
        with open(self.filename, 'w') as fil:
            for line in lines:
                if '"set"' in line: line = line[:18] + line[18:].replace('"', '')
                fil.write(line)

class Materials():
    """Create materials for particle sets.

    Attributes
    ----------
    materials : list of dict
        Each element is a dictionnary containing a material's parameters. The index of a material is his id.
    psets_ids : list of ints
        The element i of this list is the id of the particle set made of the material defined in ``materials[i]``.
    
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
    >>> materials.create_MohrCoulomb3D(pset_id=lower_particles, density=750,
    ...                                                         youngs_modulus=5.26e7,
    ...                                                         poisson_ratio=.3,
    ...                                                         friction=36.,
    ...                                                         dilation=0.,
    ...                                                         cohesion=1.,
    ...                                                         tension_cutoff=1.,
    ...                                                         softening=False)
    >>> materials.create_Newtonian3D(pset_id=upper_particles, density=1.225, 
    ...                                                       bulk_modulus=1.42e5, 
    ...                                                       dynamic_viscosity=1.81e3)
    >>> materials.pset_ids
    [0, 1]
    """
    def __init__(self): 
        self.materials = []
        self.pset_ids = []

    def create_Newtonian3D(self, pset_id=0, density=1.225, 
                                            bulk_modulus=1.42e5, 
                                            dynamic_viscosity=1.81e-5):
        """Create Newtonian3D material, as specified by CB-Geo documentation.

        Parameters
        ----------
        pset_id : int
            Particle set id that will be made off this material
        density : float, optional
            Density of the material (:math:`kg/m^3`). Default is 1.225 :math:`kg/m^3`.
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
                               "type": "Newtonian3D",
                               "density": density,
                               "bulk_modulus": bulk_modulus,
                               "dynamic_viscosity": dynamic_viscosity})
    
    def create_MohrCoulomb3D(self, pset_id=0, density=1e3,
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
        """Create MohrCoulomb3D material, as specified by CB-Geo documentation.

        Parameters
        ----------
        pset_id : int
            Particle set id that will be made off this material
        density : float, optional
            Density of the material (:math:`kg/m^3`). Default is 1.225 :math:`kg/m^3`.
        young_modulus : float, optional
            Young's modulus of the material (:math:`Pa`). Default is 50 :math:`GPa`.
        poisson_ratio : float, otpional
            Poisson's ratio of the material. Default is 0.3 .
        friction : float, optional
            Friction angle of the material (:math:`^\circ`). Default is 36 :math:`^\circ`.
        dilation : float, optional
            Dilation angle of the material (:math:`^\circ`). Default is 0 :math:`^\circ`.
        cohesion : float, optional
            Cohesion in the material (:math:`Pa`). Default is 0 :math:`Pa`.
        tension_cutoff : float, optional
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
                               "type": "MohrCoulomb3D",
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

    def create_LinearElastic3D(self, pset_id=0, density=1e3,
                                              youngs_modulus=5e7,
                                              poisson_ratio=.3):
        """Create LinearElastic3D material.

        Parameters
        ----------
        pset_id : int
            Particle set id that will be made off this material
        density : float, optional
            Density of the material (:math:`kg/m^3`). Default is 1.225 :math:`kg/m^3`.
        young_modulus : float, optional
            Young's modulus of the material (:math:`Pa`). Default is 50 :math:`GPa`.
        poisson_ratio : float, otpional
            Poisson's ratio of the material. Default is 0.3 .
        """
        self.pset_ids.append(pset_id)
        self.materials.append({"id": len(self.materials),
                               "type": "LinearElastic3D",
                               "density": density,
                               "youngs_modulus": youngs_modulus,
                               "poisson_ratio": poisson_ratio})

class Simulation():
    """Create a simulation.

    Parameters
    ----------
    title : str, optional
        Simulation title. Default is 'Sim_title'.
    directory : str, optional
        Path to the simulation's directory. Mesh, particles and entity sets files will be saved in this directory. The result folder is also set to be created by CB-Geo in this directory. Default is `title`.

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
        Initial stresses for each particle. Noting `npart` the number of particles, its shape is ``(npart, 3)``.
    input_filename : str
        Path to the input file.
    title : str
        Simulation title.
    directory : str
        Path to the simulation's directory.

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
    >>> sim.materials.create_MohrCoulomb3D(pset_id=lower_particles)
    >>> sim.materials.create_Newtonian3D(pset_id=upper_particles)
    >>> walls = []
    >>> walls.append([sim.entity_sets.create_set(lambda x,y,z: x==lim, typ="node") for lim in [0, sim.mesh.l0]])
    >>> walls.append([sim.entity_sets.create_set(lambda x,y,z: y==lim, typ="node") for lim in [0, sim.mesh.l1]])
    >>> walls.append([sim.entity_sets.create_set(lambda x,y,z: z==lim, typ="node") for lim in [0, sim.mesh.l2]])
    >>> for direction, sets in enumerate(walls): _ = [sim.add_velocity_condition(direction, 0., es) for es in sets]
    >>> sim.set_gravity([0,0,-9.81])
    >>> sim.set_analysis_parameters(dt=1e-3, nsteps=1.5e8, output_step_interval=7.5e6)
    >>> sim.write_input_file()
    """
    def __init__(self, title='Sim_title', directory=''):
        if directory == '' : directory = title
        if directory[-1] != '/' : directory += '/'
        if not os.path.isdir(directory): os.mkdir(directory)
        self.directory = directory
        self.input_filename = directory + "input_file.json"
        self.title = title
        self.materials = Materials()
        
        self.set_analysis_parameters()
        self.math_functions = []
        self.entity_sets = None
        self.init_stresses = None


        self.__boundary_conditions = {"velocity_constraints": [],
                                      "friction_constraints": [],
                                      "particles_velocity_constraints": []}
        self.__gravity = [0,0,0]
        self.__nodal_forces = []
        self.__init_stress_filename = self.directory + "particles_stresses.txt"

    def create_mesh(self, *args, **kwargs):
        """Create the simulation's mesh.

        Parameters
        ----------
        dimensions : tuple of floats
            Dimensions of the mesh. Its length should be 3, with `dimensions[n]` the dimension of the mesh on the axis `n`.
        ncells : tuple of ints
            Number of cells in each direction. Its length should be 3, with `ncells[n]` the number of cells on the axis `n`.
        check_duplicates : bool, optional
            See CB-Geo documentation for informations on this parameter. Default is `True`.
        cell_type : {'ED3H8', 'ED3H20', 'ED3H64'}, optional
            Type of cell. Only 3D Hexahedrons are supported. The number of nodes can be 8, 20 or 64. Default is 'ED3H8'.
        """
        if "directory" in kwargs or len(args)>2: raise TypeError("`directory` parameter is defined by the `Simulation` object")
        self.mesh = Mesh(*args, directory=self.directory, **kwargs)
        self.mesh.write_file()

    def create_particles(self, *args, **kwargs):
        """Create the simultation's particles.

        Parameters
        ----------
        npart_perdim_percell : int, optional
            Number of particles for each dimensions in one cell. All cells will contain ``npart_perdim_percell**3`` equally spaced particles. Note that particles are equally spaced within a cell, not between cells. Default is 1 .
        check_duplicates : bool, optional
            See CB-Geo documentation for informations on this parameter. Default is `True`.
        """
        if "mesh" in kwargs: raise TypeError("`mesh` parameter is defined by the `Simulation` object")
        if "directory" in kwargs or len(args)>1: raise TypeError("`directory` parameter is defined by the `Simulation` object")
        self.particles = Particles(mesh=self.mesh, *args, directory=self.directory, **kwargs)
        self.particles.write_file()

    def init_entity_sets(self): 
        """Create the simulation's :class:`~pycbg.preprocessing.EntitySets` object.

        Has to be called after mesh and particles creation.
        """
        self.entity_sets = EntitySets(mesh=self.mesh, particles=self.particles, directory=self.directory)

    def add_velocity_condition(self, dir, vel_value, entity_set, typ="node"):
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
        """
        if typ=="particle": list_name, key_name = "particles_velocity_constraints", "pset_id"
        elif typ=="node": list_name, key_name = "velocity_constraints", "nset_id"
        else: raise ValueError("`typ` parameter should be 'particle' or 'node'")

        self.__boundary_conditions[list_name].append({key_name: entity_set,
                                                      "dir": dir,
                                                      "velocity": vel_value})
        
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
                                                                   "velocity": frict_value})
    
    def add_math_function(self, function_wrt_time):
        """Add a math function to the simulation.

        Parameters
        ----------
        function_wrt_time : function
            Function that take the time `t` as input and returns the value of the math function.

        Notes
        -----
        Should be called after analysis parameters are set (since the time values are computed from the time step and number of steps).

        Returns
        -------
        int
            Id of the math function just appended. 
        """
        ft_values = [function_wrt_time(t) for t in (self.__analysis["dt"] * i for i in np.array(range(self.__analysis["nsteps"])))]
        fct_id = len(self.math_functions)
        
        self.math_functions.append({"id":fct_id, "type": 'Linear', "xvalues": str(list(self.__times)), "fxvalues": str(ft_values)})

        return fct_id

    def add_nodal_force(self, dir, force, node_set, math_function_id=None):
        """Add a force on all the node in a node set.

        Parameters
        ----------
        dir : {0, 1, 2}
            Axis on which the force is imposed.
        force : float
            Imposed force's value (:math:`N`).
        node_set : int
            Id of the node set on which the force is imposed.
        math_function_id : int, optional
            Id of the math function to use. Default value is `None` (the load is then static).
        """
        self.__nodal_forces.append({"nset_id": node_set,
                                    "dir": dir,
                                    "force": force})
        if math_function_id != None: self.__nodal_forces[-1]["math_function_id"] = math_function_id
    
    def set_initial_particles_stresses(self, init_stresses):
        """Set the initial stresses for each particle.

        Parameters
        ----------
        init_stresses : numpy array
            Initial stresses for each particle. Noting `npart` the number of particles, it should have the shape ``(npart, 3)``.
        """
        self.init_stresses = init_stresses

        psfile = open(self.__init_stress_filename, "w") 
        psfile.write("{:d}\n".format(len(self.particles.particles)))   
        for ps in init_stresses: psfile.write("{:e}\t{:e}\t{:e}\n".format(*ps)) 

    
    def set_gravity(self, gravity): 
        """Set the value of gravity. If this method isn't called, gravity is
        `[0,0,0]`.

        Parameters
        ----------
        gravity : list of floats
            Gravity's value on each axis (:math:`m/s^2`).
        """
        self.__gravity = gravity

    def set_analysis_parameters(self, type="MPMExplicit3D", mpm_scheme="usl", damping=0.05, locate_particles=False, dt=1e-05, nsteps=2000, output_step_interval=100):
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
        nsteps : int, optional
            Number of steps to be performed. Default is 2000.
        output_step_interval : int, optional
            Number of steps between two data points. Default is 100.
        """
        self.__analysis = {"type": type,
                           "mpm_scheme": mpm_scheme,
                           "locate_particles": locate_particles,
                           "dt": dt,
                           "damping": {"type": "Cundall", "damping_factor": damping},
                           "nsteps": int(nsteps),
                           "uuid": self.title}
        self.post_processing = {"path": self.directory + "results/",
                                "output_steps": int(output_step_interval)}
    
    def write_input_file(self):
        """Write the input file."""
        mesh_dic = {"mesh": self.mesh.filename,
                    "boundary_conditions": self.__boundary_conditions,
                    "isoparametric": self.mesh._isoparametric,
                    "check_duplicates": self.mesh.check_duplicates,
                    "cell_type": self.mesh.cell_type,
                    "io_type": self.mesh._io_type,
                    "node_type": self.mesh._node_type}
        if self.entity_sets is not None: 
            self.entity_sets.write_file()
            mesh_dic["entity_sets"] = self.entity_sets.filename
        if self.init_stresses is not None: 
            mesh_dic["particles_stresses"] = self.__init_stress_filename
        
        particles_list = [{"generator": {"check_duplicates": self.particles.check_duplicates,
                                         "location": self.particles.filename,
                                         "io_type": self.particles._io_type,
                                         "pset_id": 0,
                                         "particle_type": self.particles._particle_type,
                                         "material_id": 0,
                                         "type": self.particles._type}}]

        material_sets_list = [{"material_id": m_id, "pset_id": ps_id} for m_id, ps_id in enumerate(self.materials.pset_ids)]

        external_loading_conditions_dic = {"gravity": self.__gravity}
        if len(self.__nodal_forces) != 0: external_loading_conditions_dic["concentrated_nodal_forces"] = self.__nodal_forces
    
        dic = {"title": self.title,
               "mesh": mesh_dic,
               "particles": particles_list,
               "materials": self.materials.materials,
               "material_sets": material_sets_list,
               "external_loading_conditions": external_loading_conditions_dic,
               "analysis": self.__analysis,
               "post_processing": self.post_processing}
        
        if len(self.math_functions) != 0: dic["math_functions"] = self.math_functions

        with open(self.input_filename, 'w') as fil: json.dump(dic, fil, sort_keys=False, indent=4)

        ## Read and rewrite the file, to erase the double quotes around the list in "xvalues" and "fxvalues" lines (definitely ugly)
        with open(self.input_filename, 'r') as fil: lines = fil.readlines()
        with open(self.input_filename, 'w') as fil:
            for line in lines:
                if '"xvalues"' in line: line = line[:22] + line[22:].replace('"', '')
                if '"fxvalues"' in line: line = line[:23] + line[23:].replace('"', '')
                fil.write(line)

        save_name = self.directory + self.title + ".Simulation"
        with open(save_name, 'wb') as fil : pickle.dump(self, fil)

def setup_batch(params, directory='', cbgeo_executable=None, ncores="max"):
    if directory == '' : directory = '/'
    if directory[-1] != '/' : directory += '/'
    if not os.path.isdir(directory): os.mkdir(directory)

    set_executable = cbgeo_executable is not None

    if type(params) == list : param_sets = params
    elif type(params) == dict :
        all_combinations = it.product(*params.values())
        param_sets = [{key:val for key, val in zip(params.keys(), val_set)} for val_set in all_combinations]

    with open(os.path.basename(main.__file__), 'r') as fil: script = fil.readlines()
    for i, line in enumerate(script):
        if "setup_batch(" in line: insert_line = i
    
    
    table_file = open(directory + "parameters_sets.table", "w")
    header = "sim_id"
    for key in params.keys(): header += "\t" + key
    table_file.write(header + "\n")

    if set_executable: batch_launcher_file = open(directory + "start_batch.sh", "w")
    
    for sim_id, param_set in enumerate(param_sets):
        sim_dir = directory + "sim{:d}/".format(sim_id)
        if not os.path.isdir(sim_dir): os.mkdir(sim_dir)

        affectation_lines = [key + " = " + str(val) + "\n" for key, val in param_set.items()]
        affectation_lines += ["sim_dir = '{:}'\n".format(sim_dir)]
        out_script = script[:insert_line] + affectation_lines + script[insert_line+1:]

        out_script_path = sim_dir + "pycbg_script.py"
        with open(out_script_path, "w") as fil: 
            for line in out_script: fil.write(line)
        runpy.run_path(out_script_path)
        #exec(open(out_script_path).read(), globals(), locals())
        
        param_line = str(sim_id)
        for val in param_set.values(): param_line += "\t" + str(val)
        table_file.write(param_line + "\n")

        if set_executable: 
            if ncores!="max": cores_str = " -p {:d}".format(ncores)
            else : cores_str = ""
            batch_launcher_file.write("""{:}{:} -f "$(pwd)/" -i {:}input_file.json >> {:}cbgeo.log\n""".format(cbgeo_executable, cores_str, sim_dir, sim_dir))
    
    table_file.close()
    if set_executable: batch_launcher_file.close()
    sys.exit()


    

