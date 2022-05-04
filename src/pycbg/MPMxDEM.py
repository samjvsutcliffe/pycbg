import sys, os
import glob as gb, pickle

## Beware, anything defined globally in this module (except variable whose names are in the no_auto_import list) is also imported in the main script (the one importing this module) upon calling __update_imports (which is called by several functions of this module)

no_auto_import = ["glob"]

# Default value for optional parameters
state_variables = []
vtk_period = 0
state_vars = ["O.iter, O.time, O.dt"]
save_final_state = False

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
    """Import YADE in the current scope and create a directory for all RVEs data (samples, vtk files, ...). This directory's path is stored in `rve_directory`, automatically imported in the main script. YADE's version is stored as the last commit SHA1 in `yade_sha1`, also automatically imported in the main script. The YADE simulation should be set as periodic by the user (non periodic simulations are not supported).

    Parameters
    ----------
    yade_exec : str
        Full path to the YADE executable. 
    """
    global rve_directory, pycbg_sim, yade_sha1

    # Get PyCBG simulation object
    pycbg_sim = pickle.load(gb.glob(os.environ["PWD"] + "*.Simulation")[0])

    # Load YADE
    exec_path, exec_name = yade_exec.rstring("/", 1)

        ## Add exec_path to path
    sys.path.append(exec_path)

        ## Perform an 'import *' on yade module. The following is NOT a good practice (tempering with the vars dict), I should try to improve this (thanks https://stackoverflow.com/a/11007138/16796697)
    for key, val in vars(__import__(exec_name)).items():
        if key.startswith('__') and key.endswith('__'): continue
        vars()[key] = val

        ## Get git's SHA1
    yade_sha1 = version.split("-")[-1]

        ## Print all versions to a file
    if not os.path.isfile(pycbg_sim.directory + 'yade_all_versions.txt'):
        original_stdout = sys.stdout 
        with open(pycbg_sim.directory + 'yade_all_versions.txt', 'w') as f:
            sys.stdout = f 
            printAllVersions()
            sys.stdout = original_stdout

    # Create rve_data directory
    rve_directory = pycbg_sim.directory + "rve_data/"
    if not os.path.isdir(rve_directory): os.mkdir(rve_directory)

    # Update variables in main script
    loc = locals().copy()
    for loc_var in ["yade_exec", "exec_path", "exec_name", "pycbg_sim"]: del loc[loc_var]
    __update_imports(loc)


def set_opt_params(vtk_period=0, state_vars=["O.iter, O.time, O.dt"], save_final_state=False):
    """Set MPMxDEM optional parameters. `MPMxDEM.set_opt_params` should be called before `MPMxDEM.define_compute_stress`.

    Parameters
    ----------
    vtk_period : int
        `iterPeriod` for YADE's `VTKRecorder` engine. Its value is stored in `vtk_p`, automatically imported in the main script.
    state_vars : list of str
        List of python expressions that should return a scalar, to be save as state variable in the CB-Geo simulation. `state_vars` should have at most 19 elements, the first element being called `svars_1` in CB-Geo, the second `svars_2`, ... . Default to `["O.iter, O.time, O.dt"]`
    save_final_state : bool
        Wether or not to save the RVE final state in a ".{SHA1}yade.bz2" file, where "{SHA1}" is git's last commit SHA1 of YADE. Its value is stored in `save_fstate`, automatically imported in the main script.

    """
    global state_variables, vtk_p, save_fstate
    state_variables, vtk_p, save_fstate = state_vars, vtk_period, save_final_state


def define_compute_stress(dem_strain_rate, function_name="compute_stress", run_on_setup=None):
    """Define the function to be called at each MPM step, with CB-Geo's required signature. This function sets the global variable `rve_id`, which corresponds to the RVE's particle id in CB-Geo. User has to `MPMxDEM.set_opt_params` should be called before `MPMxDEM.define_compute_stress`.

    Parameters
    ----------
    dem_strain_rate : float
        Strain rate applied to the RVE.
    function_name : str
        Name of the function to be defined (also the one to be passed to PyCBG's PythonModel3D material). Default is `"compute_stress"`.
    run_on_setup : str or None
        Name of the function to be run on RVE setup, if not None. This function is called after `rve_id` is defined, `run_on_setup` can thus refer to it.
    """
            
    def fct(rid, de_xx, de_yy, de_zz, de_xy, de_yz, de_xz, mpm_iteration, *state_vars):
        global rve_id, state_variables

        # Use usual strain, not the engineering one computed by CB-Geo
        de_xy, de_yz, de_xz = .5*de_xy, .5*de_yz, .5*de_xz

        # If this function is called for the first time
        if mpm_iteration==0:
            ## Keep rve_id
            rve_id = rid

            ## Run user's setup function
            globals()[run_on_setup]()

            ## Create RVE directory
            vtk_dir = rve_directory + "RVE_{:}/".format(rve_id) 
            if not os.path.isdir(vtk_dir): os.mkdir(vtk_dir)

            ## Add VTKRecorder to engines
            O.engines += [VTKRecorder(fileName=vtk_dir, recorders=["all"], iterPeriod=vtk_p)]

        # Shaping dstrain increment matrix
        dstrain_matrix = Matrix3((de_xx, de_xy, de_xz,
                                    de_xy, de_yy, de_yz,
                                    de_xz, de_yz, de_zz))
        
        # Compute the DEM deformation time to keep the simulation quasistatic 
        max_deps = max([abs(i) for i in [de_xx, de_yy, de_zz, de_xy, de_yz, de_xz]])
        deformation_time = max_deps / dem_strain_rate

        # Compute the number of DEM iterations
        n_dem_iter = -(-deformation_time//O.dt)
        deformation_time = n_dem_iter * O.dt # increase a little deformation_time so the number of iteration is exactly an integer

        # Compute the velocity gradient, assuming no rotation
        O.cell.velGrad = dstrain_matrix / deformation_time

        # Measure initial stress
        sigma0 = getStress(O.cell.volume)

        # Run DEM steps
        O.run(n_dem_iter)
        
        # Finnish the MPM iteration
        mpm_iteration += 1
        dsigma = getStress(O.cell.volume)-sigma0

        # Update state variables
        state_vars = [eval(var, globals()) for var in state_variables]

        # Save final state
        if mpm_iteration == pycbg_sim.analysis_params["nsteps"] and save_fstate:
            O.save(rve_directory + "rve{:d}_final_state.{:}yade.bz2".format(rve_id, yade_sha1))

        return (dsigma[0,0], dsigma[1,1], dsigma[2,2], dsigma[0,1], dsigma[1,2], dsigma[0,2], mpm_iteration) + tuple(state_vars)

    # Give the function the name chosen by the user
    function_name = fct

    # Update variables in main script
    loc = locals().copy()
    for loc_var in ["dem_strain_rate", "function_name", "run_on_setup", "fct"]: del loc[loc_var]
    __update_imports(loc)