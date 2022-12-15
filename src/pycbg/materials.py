import gmsh, os, json, pickle, csv, runpy, sys, shutil
import multiprocessing
import numpy as np
import itertools as it  

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

