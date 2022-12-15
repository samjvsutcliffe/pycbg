import gmsh, os, json, pickle, csv, runpy, sys, shutil
import multiprocessing
import numpy as np
import itertools as it  

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
                if mesh.cell_type=="ED2Q36BS": cell_nodes = [e[i] for i in [14, 15, 20, 21]]
                else: cell_nodes = e
                coors = np.array([mesh.nodes[i] for i in cell_nodes])
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


