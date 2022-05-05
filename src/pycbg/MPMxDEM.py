import sys, os
import glob as gb, pickle, numpy as np

## Beware, anything defined globally in this module (except variable whose names are in the no_auto_import list) is also imported in the main script (the one importing this module) upon calling __update_imports (which is called by several functions of this module)

no_auto_import = ["glob", "rve_directory", "pycbg_sim", "yade_sha1"]

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

def setup_yade(yade_exec="/usr/bin/yade"):
    """Import YADE in the current scope and create a directory for all RVEs data (samples, vtk files, ...). The YADE simulation should be set as periodic by the user (non periodic simulations are not supported).

    Parameters
    ----------
    yade_exec : str
        Full path to the YADE executable. 
    """
    global rve_directory, pycbg_sim, yade_sha1

    # Get PyCBG simulation object
    with open(gb.glob("*.Simulation")[0], 'rb') as fil : pycbg_sim = pickle.load(fil)

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
    yade_sha1 = version.split("-")[-1]

        ## Print all versions to a file
    if not os.path.isfile('yade_all_versions.txt'):
        original_stdout = sys.stdout 
        with open('yade_all_versions.txt', 'w') as f:
            sys.stdout = f 
            printAllVersions()
            sys.stdout = original_stdout

    # Create rve_data directory
    rve_directory = "rve_data/"
    if not os.path.isdir(rve_directory): os.mkdir(rve_directory)

    # Update variables in main script
    loc = locals().copy()
    for loc_var in ["yade_exec", "exec_path", "exec_name"]: del loc[loc_var]
    __update_imports(loc)

class DefineCallable():
    """Callable object to be called at each MPM step, with CB-Geo's required signature. The YADE periodic simulation defined in the script creating this callable object will be deformed at each MPM iteration using `O.cell.velGrad`. The velocity gradient is computed using the strain increment provided by CB-Geo and the `dem_strain_rate` parameter provided by the user: `O.cell.velGrad = strain_increment_matrix / max(strain_increment_matrix) * dem_strain_rate`.

    Parameters
    ----------
    dem_strain_rate : float
        Strain rate applied to the RVE.
    run_on_setup : callable or None
        Name of the function to be run on RVE setup, if not None. This function is called after `rve_id` is defined, `run_on_setup` can thus refer to it.
    vtk_period : int
        `iterPeriod` for YADE's `VTKRecorder` engine. Default is 0 (no VTK file is saved).
    state_vars : list of str
        List of python expressions that should return a scalar, to be save as state variable in the CB-Geo simulation. `state_vars` should have at most 19 elements, the first element being called `svars_1` in CB-Geo, the second `svars_2`, ... . Default to `["O.iter, O.time, O.dt"]`.
    save_final_state : bool
        Wether or not to save the RVE final state in a ".{SHA1}yade.bz2" file, where "{SHA1}" is git's last commit SHA1 of YADE. Default is `False`.


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
        The 'particle_id' of the current RVE, as numbered by CB-Geo. Before the first call of the object by CB-Geo, `rve_id=nan`, it is set to the actual particle id right before callin `run_on_setup`.
    rve_directory : str
        Path to the directory containing all RVEs data (samples, vtk files, ...), which is a subdirectory of PyCBG's simulation directory.
    pycbg_sim : :class:`~pycbg.preprocessing.Simulation` object
        PyCBG's simulation object used to create the input files.
    yade_sha1 : str
        Partial SHA1 of YADE's version
    """

    def __init__(self, dem_strain_rate, run_on_setup=None, vtk_period=0, state_vars=["O.iter, O.time, O.dt"], svars_dic={}, save_final_state=False): 
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

    def __call__(self, rid, de_xx, de_yy, de_zz, de_xy, de_yz, de_xz, mpm_iteration, *state_vars):

        # Use usual strain, not the engineering one computed by CB-Geo
        de_xy, de_yz, de_xz = .5*de_xy, .5*de_yz, .5*de_xz

        # If this function is called for the first time
        if mpm_iteration==0:
            ## Keep rve_id
            self.rve_id = int(rid)

            ## Run user's setup function
            self.run_on_setup()

            ## Create RVE directory
            vtk_dir = rve_directory + "RVE_{:}/".format(self.rve_id) 
            if not os.path.isdir(vtk_dir): os.mkdir(vtk_dir)

            ## Add VTKRecorder to engines
            O.engines += [VTKRecorder(fileName=vtk_dir, recorders=["all"], iterPeriod=self.vtk_period)]

        # Shaping dstrain increment matrix
        dstrain_matrix = Matrix3((de_xx, de_xy, de_xz,
                                    de_xy, de_yy, de_yz,
                                    de_xz, de_yz, de_zz))
        
        # Compute the DEM deformation time to keep the simulation quasistatic 
        max_deps = max([abs(i) for i in [de_xx, de_yy, de_zz, de_xy, de_yz, de_xz]])
        deformation_time = max_deps / self.dem_strain_rate

        # Compute the number of DEM iterations
        n_dem_iter = -(-deformation_time//O.dt)
        deformation_time = n_dem_iter * O.dt # increase a little deformation_time so the number of iteration is exactly an integer

        # Compute the velocity gradient, assuming no rotation
        O.cell.velGrad = dstrain_matrix / deformation_time

        # Measure initial stress
        sigma0 = getStress(O.cell.volume)

        # Run DEM steps
        for i in range(int(n_dem_iter)): O.step()
        
        # Finnish the MPM iteration
        mpm_iteration += 1
        dsigma = getStress(O.cell.volume)-sigma0

        # Update state variables
        state_vars = [eval(var, self.svars_dic) for var in self.state_variables]

        # Save final state
        if mpm_iteration == pycbg_sim.analysis_params["nsteps"] and self.save_fstate:
            O.save(rve_directory + "rve{:d}_final_state.{:}yade.bz2".format(self.rve_id, self.yade_sha1))

        return (dsigma[0,0], dsigma[1,1], dsigma[2,2], dsigma[0,1], dsigma[1,2], dsigma[0,2], mpm_iteration) + tuple(state_vars)