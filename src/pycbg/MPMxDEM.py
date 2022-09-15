import sys, os, warnings, inspect
import glob as gb, pickle, numpy as np, itertools as it
from matplotlib import use
import pycbg

## Beware, anything defined globally in this module (except variable whose names are in the no_auto_import list) is also imported in the main script (the one importing this module) upon calling __update_imports (which is called by several functions of this module)

no_auto_import = ["glob", "rve_directory", "pycbg_sim", "yade_sha1"]

#script_dir = os.path.dirname(os.path.realpath(__file__))
#os.chdir(inspect.stack()[-1].filename.rsplit("/", 1)[0]) # not sure about portability to other systemsz

# Initialise glob dictionary
glob = set(globals())

def __update_imports(custom_dict={}):
    global glob

    # Create a dictionary wwith all variables to be imported by main module
    all_vars = {**globals(), **custom_dict}

    # Get new keys with respect to last call
    new_keys = list(glob ^ set(all_vars))

    # Import all non-system variable into main script
    for key in new_keys: 
        if key not in no_auto_import: sys.modules['builtins'].__dict__[key] = all_vars[key]

    # Update glob dictionary, so at next call only new items will be processed
    glob = set(globals())

def __on_yade_setup():
    global rve_directory, pycbg_sim, yade_sha1
    # Get PyCBG simulation object
    with open(gb.glob("*.Simulation")[0], 'rb') as fil : pycbg_sim = pickle.load(fil)

    # Create rve_data directory
    rve_directory = "rve_data/"
    if not os.path.isdir(rve_directory): os.mkdir(rve_directory)

    # Update variables in main script
    loc = locals().copy()
    __update_imports(loc)

def setup_yade(yade_exec="/usr/bin/yade"):
    """Import YADE in the current scope and create a directory for all RVEs data (samples, vtk files, ...). The YADE simulation should be set as periodic by the user (non periodic simulations are not supported).

    Parameters
    ----------
    yade_exec : str
        Full path to the YADE executable. 
    """
    global yade_sha1

    # Load YADE
    exec_path, exec_name = yade_exec.rsplit("/", 1)

        ## Add exec_path to path
    sys.path.append(exec_path)

        ## Perform an 'import *' on yade module. The following is NOT a good practice (tempering with the vars dict), I should try to improve this (thanks https://stackoverflow.com/a/11007138/16796697)
    for key, val in vars(__import__(exec_name)).items():
        if key.startswith('__') and key.endswith('__'): continue
        vars()[key] = val
        globals()[key] = val

        ## Get git's SHA1
    try: yade_sha1 = version.split("-")[-1]
    except NameError: yade_sha1 = "release-"

        ## Print all versions to a file
    if not os.path.isfile('yade_all_versions.txt'):
        original_stdout = sys.stdout 
        with open('yade_all_versions.txt', 'w') as f:
            sys.stdout = f 
            printAllVersions()
            sys.stdout = original_stdout
    
    # Print pycbg version to a file
    if not os.path.isfile('pycbg_version.txt'):
        original_stdout = sys.stdout 
        with open('pycbg_version.txt', 'w') as f:
            sys.stdout = f 
            print(pycbg.__version__)
            sys.stdout = original_stdout

    try: __on_yade_setup()
    except: warnings.warn("Extra setup steps coudln't be performed, the current session is then a simple YADE session")

    # Update variables in main script
    loc = locals().copy()
    for loc_var in ["yade_exec", "exec_path", "exec_name"]: del loc[loc_var]
    __update_imports(loc)

class DefineCallable():
    """Callable object to be called at each MPM step, with CB-Geo's required signature. The YADE periodic simulation defined in the script creating this callable object will be deformed at each MPM iteration using `O.cell.velGrad`. The velocity gradient is computed using the strain increment provided by CB-Geo and the `dem_strain_rate` parameter provided by the user: `O.cell.velGrad = strain_increment_matrix / max(strain_increment_matrix) * dem_strain_rate`.

    Parameters
    ----------
    dem_strain_rate : float or None
        Strain rate applied to the RVE. If `None`, the MPM strain rate is used.
    run_on_setup : callable or None
        Name of the function to be run on RVE setup, if not None. This function is called after `rve_id` is defined, `run_on_setup` can thus refer to it.
    vtk_period : int
        `iterPeriod` for YADE's `VTKRecorder` engine. Default is 0 (no VTK file is saved).
    state_vars : list of str
        List of python expressions that should return a scalar, to be save as state variable in the CB-Geo simulation. `state_vars` should have at most 19 elements, the first element being called `svars_1` in CB-Geo, the second `svars_2`, ... . Default to `["O.iter, O.time, O.dt"]`.
    svars_dic : dict
        Dictionary in which all elements of `state_vars`are evaluated. Most users need to set `svars_dic=globals()`.
    save_final_state : bool
        Wether or not to save the RVE final state in a ".{SHA1}yade.bz2" file, where "{SHA1}" is git's last commit SHA1 of YADE. Default is `False`.
    flip_cell_period : int
        YADE's `flipCell` is called every `flip_cell_period`, which flips the RVE toward more axis-aligned base cell vectors if it is possible.
    use_gravity : bool
        Wether or not to compute the DEM cell's global stress considering gravity

    Attributes
    ----------
    dem_strain_rate : float
        Strain rate applied to the RVE.
    run_on_setup : callable or None
        Name of the function to be run on RVE setup, if not None. This function is called after `rve_id` is defined, `run_on_setup` can thus refer to it.
    vtk_period : int
        `iterPeriod` for YADE's `VTKRecorder` engine. Default is 0 (no VTK file is saved).
    state_variables : list of str
        List of python expressions that should return a scalar, to be save as state variable in the CB-Geo simulation. `state_vars` should have at most 19 elements, the first element being called `svars_1` in CB-Geo, the second `svars_2`, ... . Default to `["O.iter, O.time, O.dt"]`.
    save_final_state : bool
        Wether or not to save the RVE final state in a ".{SHA1}yade.bz2" file, where "{SHA1}" is git's last commit SHA1 of YADE. Default is `False`.
    rve_id : int
        The 'particle_id' of the current RVE, as numbered by CB-Geo. Before the first call of the object by CB-Geo, `rve_id=nan`, it is set to the actual particle id right before calling `run_on_setup`.
    rve_directory : str
        Path to the directory containing all RVEs data (samples, vtk files, ...), which is a subdirectory of PyCBG's simulation directory.
    pycbg_sim : :class:`~pycbg.preprocessing.Simulation` object
        PyCBG's simulation object used to create the input files.
    yade_sha1 : str
        Partial SHA1 of YADE's version
    flip_cell_period : int
        YADE's `flipCell` is called every `flip_cell_period`, which flips the RVE toward more axis-aligned base cell vectors if it is possible.
    flip_count : int
        Number of time the DEM cell has been flipped
    use_gravity : bool
        Wether or not to compute the DEM cell's global stress considering gravity
    mpm_dt : float
        The MPM time step for the current simulation
    mpm_iter : int
        The current MPM iteration
    dstrain : numpy array of shape (3,3)
        Strain increment for the current MPM iteration
    dstress : numpy array of shape (3,3)
        Stress increment for the current MPM iteration
    """

    def __init__(self, dem_strain_rate, run_on_setup=None, vtk_period=0, state_vars=["O.iter, O.time, O.dt"], svars_dic={}, save_final_state=False, flip_cell_period=0, use_gravity=False): 
        self.dem_strain_rate = dem_strain_rate
        self.run_on_setup = run_on_setup
        self.vtk_period = vtk_period
        self.state_variables = state_vars
        self.svars_dic = svars_dic
        self.save_final_state = save_final_state
        self.rve_directory = rve_directory
        self.pycbg_sim = pycbg_sim
        self.yade_sha1 = yade_sha1
        self.rve_id = np.nan
        self.flip_cell_period = flip_cell_period
        self.flip_count = 0
        self.use_gravity = use_gravity
        self.mpm_iter = 0
        self.mpm_dt = pycbg_sim.analysis_params["dt"]
        self.dstrain = np.zeros((3,3))
        self.dstress = np.zeros((3,3))

    def __call__(self, rid, de_xx, de_yy, de_zz, de_xy, de_yz, de_xz, mpm_iteration, *state_vars):
        global sigma0

        # Update mpm_iter attribute
        self.mpm_iter = mpm_iteration

        # Use usual strain, not the engineering one computed by CB-Geo
        de_xy, de_yz, de_xz = .5*de_xy, .5*de_yz, .5*de_xz

        # Set the dstrain attribute
        self.dstrain = np.array([[de_xx, de_xy, de_xz], [de_xy, de_yy, de_xz], [de_xz, de_yz, de_zz]])

        # If this function is called for the first time
        if mpm_iteration==0:
            ## Keep rve_id
            self.rve_id = int(rid)

            ## Run user's setup function
            if self.run_on_setup is not None: self.run_on_setup()

            ## Create RVE directory
            vtk_dir = rve_directory + "RVE_{:}/".format(self.rve_id) 
            if not os.path.isdir(vtk_dir): os.mkdir(vtk_dir)

            ## Add VTKRecorder to engines
            if self.vtk_period!=0: O.engines += [VTKRecorder(fileName=vtk_dir, recorders=["all"], iterPeriod=self.vtk_period)]

            ## Measure initial stress
            if not self.use_gravity: sigma0 = getStress(O.cell.volume)
            else: 
                ### If gravity is used, the sample global stress has to be computed manually, a list of particles and walls are thus
                _get_bodies_walls()
                sigma0 = _getStress_gravity()

        # Shaping dstrain increment matrix
        dstrain_matrix = Matrix3((de_xx, de_xy, de_xz,
                                    de_xy, de_yy, de_yz,
                                    de_xz, de_yz, de_zz))
        
        # Compute the DEM deformation time to keep the simulation quasistatic 
        max_deps = max([abs(i) for i in [de_xx, de_yy, de_zz, de_xy, de_yz, de_xz]])
        deformation_time = max_deps / self.dem_strain_rate if self.dem_strain_rate is not None else self.mpm_dt

        # Compute the number of DEM iterations
        n_dem_iter = -(-deformation_time//O.dt)
        deformation_time = n_dem_iter * O.dt # increase a little deformation_time so the number of iteration is exactly an integer

        # Compute the velocity gradient, assuming no rotation
        O.cell.velGrad = dstrain_matrix / deformation_time

        # Run DEM steps
        for i in range(int(n_dem_iter)): 
            if self.flip_cell_period>0: 
                if O.iter % self.flip_cell_period == 0:
                    if flipCell()!=Matrix3(0,0,0, 0,0,0, 0,0,0): self.flip_count += 1
            O.step()
        
        # Finnish the MPM iteration
        mpm_iteration += 1
        if not self.use_gravity: new_stress = getStress(O.cell.volume)
        else: new_stress = _getStress_gravity()
        dsigma = new_stress - sigma0
        sigma0 = new_stress

        # Set the dstress attribute
        self.dstress = np.array(dsigma)

        # Update state variables
        state_vars = [eval(var, self.svars_dic) for var in self.state_variables]

        # Save final state
        if mpm_iteration == pycbg_sim.analysis_params["nsteps"] and self.save_final_state:
            O.save(rve_directory + "RVE_{:}/".format(self.rve_id) + "rve{:d}_final_state.{:}yade.bz2".format(self.rve_id, self.yade_sha1))

        return (dsigma[0,0], dsigma[1,1], dsigma[2,2], dsigma[0,1], dsigma[1,2], dsigma[0,2], mpm_iteration) + tuple(state_vars)

def _get_bodies_walls():
    global bodies_id, walls_id
    bodies_id, walls_id = [], []
    for b in O.bodies:
        if type(b.shape) in [Facet, Box, Wall]: walls_id.append(b.id)
        else: bodies_id.append(b.id)

def _getStress_gravity():
    volume = O.cell.volume

    for e in O.engines:
        if type(e)==NewtonIntegrator: gravity = e.gravity
     
    ## Interaction contribution
    sigma_a = np.zeros((3,3))
    for inter in O.interactions:
        if inter.id1 in walls_id or inter.id2 in walls_id: # If one of the the bodies is a wall
            if inter.id1 in walls_id and inter.id2 not in walls_id: sgn = 1 # If only body 1 is a wall
            elif inter.id2 in walls_id and inter.id1 not in walls_id: sgn = -1 # If only body 2 is a wall
            
            f, x = sgn*(inter.phys.normalForce + inter.phys.shearForce), inter.geom.contactPoint
            for i, fi in enumerate(f):
                for j, xj in enumerate(x): sigma_a[i,j] += fi*xj/volume
            
    ## Gravity contribution
    sigma_b = np.zeros((3,3))
    for b_id in bodies_id:
        b = O.bodies[int(b_id)]
        m, x = b.state.mass, b.state.pos
        for i, gi in enumerate(gravity):
            for j, xj in enumerate(x): sigma_b[i,j] += m*gi*xj/volume

    return sigma_a + sigma_b
