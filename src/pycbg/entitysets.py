import gmsh, os, json, pickle, csv, runpy, sys, shutil
import multiprocessing
import numpy as np
import itertools as it  

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

        self.psets, self.nsets, self.csets = [], [], []
        self.__reset_params()
    
    def create_set(self, condition_function, typ="particle", kwargs={}):
        """Create a set of nodes or particles and add it to the corresponding
        list. Nodes and particles are selected using `condition_function`.

        Parameters
        ----------
        condition_function : function
            Function to select particles or nodes using their positions. It should take as inputs at least 3 parameters `x`, `y` and `z` (2 parameters `x` and `y` for 2D analyses) that correspond to the position of a node or particle, additional keyword parameters can be passed through `kwargs`. It should return `True` if the node or particle belongs to the set, `False` otherwise.
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
        elif typ=="cell": points, set_list = self.mesh.cells, self.csets
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
        for typ, sets_tmp in zip(("particle_sets", "node_sets", "cell_sets"), (self.psets, self.nsets,self.csets)):
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

