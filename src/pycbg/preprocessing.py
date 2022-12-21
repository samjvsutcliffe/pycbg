import gmsh, os, json, pickle, csv, runpy, sys, shutil
import multiprocessing
import numpy as np
import itertools as it  
import __main__ as main
from pycbg import __version__ as pycbg_version
from pycbg.mesh import Mesh
from pycbg.materials import Materials 
from pycbg.particles import Particles
from pycbg.entitysets import EntitySets 

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
        Initial volume for each particle. Noting `npart` the number of particles, its shape is ``(npart,)``.
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
            Type of set on which the velocity is imposed. Default is "node".
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
            Type of set on which the force is imposed. Default is "node".
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
            Initial volumes for each particle. Noting `npart` the number of particles, it should have a ``(npart, 1)`` or ``(npart,)`` shape.
        """
        if len(init_volumes.shape)>1:
            if len(init_volumes.shape)>2 or init_volumes.shape[1]>1: raise ValueError("Incorrect shape for given volumes array, please check your inputs with respect to documentation")
            else: self.init_volumes = init_volumes.reshape(init_volumes.shape[0]) 
        else: self.init_volumes = init_volumes

        psfile = open(self.__init_volumes_filename, "w") 
        for i, ps in enumerate(self.init_volumes): psfile.write("{:d}\t{:e}\n".format(i, ps)) 

    
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
            Cundall's damping. Should verify : ``0. <= damping < 1``. Default is 0.05 .
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
        if _type(damping) == float or _type(damping) == int: damping_param = {"type": "Cundall", "damping_factor": damping}
        elif _type(damping) == dict: damping_param = damping
        else: raise ValueError("`damping` parameter wasn't correctly set: it has to be a float or a dictionnary, please check your script")

        self.analysis_params = {"type": detected_type,
                           "mpm_scheme": mpm_scheme,
                           "locate_particles": locate_particles,
                           "dt": dt,
                           "damping": damping_param,
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
        elif self.particles.type == "inject":
            particles_list = [{"generator": {"check_duplicates": self.particles.check_duplicates,
                                            "pset_id": 0,
                                            "cset_id": self.particles.cset_id, 
                                            "nparticles_per_dir": self.particles.nparticles_per_dir,
                                            "particle_type": self.particles._particle_type,
                                            "velocity" : self.particles._particle_velocity,
                                            "duration" : self.particles._particle_duration,
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

        batch_launcher_file.write("""echo -e "\e[90m$(date +'[%Y-%m-%d %T.%N]')\e[0m\tStarting batch from directory {:}"\n\n""".format(directory))
        batch_launcher_file.write("""echo -e "\tPyCBG version: {:}\n\tNumber of CPUs to be used: {:d}\n\tNumber of CPUs per job: {:d}\n"\n\n""".format(pycbg_version, ncores, ncores_perjob))
        batch_launcher_file.write("""startT=$(date +%s)\n\n""")

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
            batch_launcher_file.write("""( startT{:}=$(date +%s); echo -e "\e[90m$(date +'[%Y-%m-%d %T.%N]')\e[0m\tStarting {:};\tlog file: {:}"; {:}{:} -f "$(pwd)/" -i {:}input_file.json >> {:}cbgeo.log; t{:}=$((( $(date +%s) - $startT{:} ))); if test $t{:} -gt $(((3600*24))); then days{:}="%jd-"; t{:}=$((($t{:}-3600*24))); else days{:}=""; fi; echo -e "\e[90m$(date +'[%Y-%m-%d %T.%N]')\e[0m\tFinished {:};\tduration: $(date -u --date @$t{:} +$days{:}%H:%M:%S)"; . /tmp/n_sim_slot_available.sh; n_sim_slot_available=$((++n_sim_slot_available)); echo "n_sim_slot_available=$n_sim_slot_available" > /tmp/n_sim_slot_available.sh ) &\n""".format(sim_title, sim_title, sim_dir+"cbgeo.log", cbgeo_executable, cores_str, sim_dir, sim_dir, sim_title, sim_title, sim_title, sim_title, sim_title, sim_title, sim_title, sim_title, sim_title, sim_title))

    table_file.close()
    if set_executable: 
        batch_launcher_file.write("""for job in `jobs -p`; do wait $job; done; tall=$((( $(date +%s) - $startT ))); if test $tall -gt $(((3600*24))); then daysall="%jd-"; tall=$((($tall-3600*24))); else daysall=""; fi; echo -ne "\n\e[90m$(date +'[%Y-%m-%d %T.%N]')\e[0m\tBatch over, total duration: $(date -u --date @$tall +$daysall%H:%M:%S)\n" """)
        batch_launcher_file.close()
    sys.exit()

def _type(*args, **kwargs): # Just so the variable name "type" can be used as function argument, and still keep access to the "type" built-in function
    return type(*args, **kwargs)
