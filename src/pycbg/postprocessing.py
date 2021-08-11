import os, pickle
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d import Axes3D

class ResultsReader():
    """Load the result of a simulation. Can also load the :func:`~pycbg.preprocessing.Simulation` object used during preprocessing.

    Parameters
    ----------
    directory : str
        Directory in which the input file of the simulation was saved.

    Attributes
    ----------
    ppositions : list of numpy arrays
        Particles' positions for every saved steps. Noting npart the number of particles in the simulations at the ith step, the shape of ``ppositions[i]`` is ``(npart,3)``.
    pvelocities : list of numpy arrays
        Particles' velocities for every saved steps. Noting npart the number of particles in the simulations at the ith step, the shape of ``pvelocities[i]`` is ``(npart,3)``.
    pstresses : list of numpy arrays
        Particles' stresses for every saved steps. Noting npart the number of particles in the simulations at the ith step, the shape of ``pstresses[i]`` is ``(npart,6)``. The columns represent the directions `xx`, `yy`, `zz`, `xy`, `yz` and `xz` respectively.
    pstrains : list of numpy arrays
        Particles' strains for every saved steps. Noting npart the number of particles in the simulations at the ith step, the shape of ``pstrains[i]`` is ``(npart,6)``. The columns represent the directions `xx`, `yy`, `zz`, `xy`, `yz` and `xz` respectively.
    ppressures : list of numpy arrays
        Particles' pressures for every saved steps. Noting npart the number of particles in the simulations at the ith step, the shape of ``ppressures[i]`` is ``(npart,)``.
    pmasses : list of numpy arrays
        Particles' masses for every saved steps. Noting npart the number of particles in the simulations at the ith step, the shape of ``pmasses[i]`` is ``(npart,)``.
    pvolumes : list of numpy arrays
        Particles' volumes for every saved steps. Noting npart the number of particles in the simulations at the ith step, the shape of ``pvolumes[i]`` is ``(npart,)``.
    pmaterials : list of numpy arrays
        Particles' material's id for every saved steps. Noting npart the number of particles in the simulations at the ith step, the shape of ``pmaterials[i]`` is ``(npart,)``.
    raw_data : list of pandas dataframes
        All data saved from the simulation. The data is stored for each time step as a dataframe.
    steps : list of ints
        All saved steps. 
    times : list of floats
        All times corresponding to saved steps.
    """
    def __init__(self, directory):
        if directory[-1] != '/': directory += '/'
        self.main_dir = directory
        self.res_dir = directory + "results/"

        with open(self.main_dir + "input_file.json", 'r') as fil: lines = fil.readlines()
        for l in lines:
            if '"title"' in l: self.title = l.split('"')[-2]
            if '"dt"' in l: dt = float(l.split(' ')[-1][:-2])

        self.data_dir = self.res_dir + self.title + '/'
        
        self.__extract_data()

        self.times = [dt*step for step in self.steps]
    
    def load_simulation(self):
        """Load the simulation object used to write the input files"""
        save_name = self.main_dir + self.title + ".Simulation"
        with open(save_name, 'rb') as fil : sim = pickle.load(fil)
        
        return sim
    
    def __extract_data(self):
        files = os.listdir(self.data_dir)
        files = [f for f in files if f[-3:]=='.h5']

        raw_data, steps = [], []
        for f in files: 
            raw_data.append(pd.read_hdf(self.data_dir + f, 'table'))
            steps.append(int(f[9:-3]))

        sort_ind = np.argsort(steps)
        self.raw_data, self.steps = [], []
        for i in sort_ind: 
            self.raw_data.append(raw_data[i])
            self.steps.append(steps[i])

        self.ppositions, self.pvelocities = [], []
        self.pstresses, self.pstrains, self.ppressures = [], [], []
        self.pmasses, self.pvolumes, self.pmaterials = [], [], []
        for df in self.raw_data : 
            self.ppositions.append(np.array([df[key] for key in ["coord_x", "coord_y", "coord_z"]]).T)
            self.pvelocities.append(np.array([df[key] for key in ["velocity_x", "velocity_y", "velocity_z"]]).T)
            self.pstresses.append(np.array([df[key] for key in ["stress_xx", "stress_yy", "stress_zz", "tau_xy", "tau_yz", "tau_xz"]]).T)
            self.pstrains.append(np.array([df[key] for key in ["strain_xx", "strain_yy", "strain_zz", "gamma_xy", "gamma_yz", "gamma_xz"]]).T)
            self.ppressures.append(df['pressure'].values)
            self.pmaterials.append(df['material_id'].values)
            self.pvolumes.append(df['volume'].values)
            self.pmasses.append(df['mass'].values)
    
__ind_sgmts = [[0,1], [1,2], [2,3], [3,0],
               [4,5], [5,6], [6,7], [7,4],
               [0,4], [1,5], [2,6], [3,7]]

def plot_mesh(mesh, fig=None, ax=None):
    if fig is None and ax is None: 
        fig = plt.figure(1, figsize=(15,10))
        ax = fig.add_subplot(111, projection='3d')

    added_segments = []
    for cell in mesh.cells:
        nodes = [mesh.nodes[i_node] for i_node in cell]
        xs = [[nodes[i][0], nodes[j][0]] for i,j in __ind_sgmts]
        ys = [[nodes[i][1], nodes[j][1]] for i,j in __ind_sgmts]
        zs = [[nodes[i][2], nodes[j][2]] for i,j in __ind_sgmts]
        
        for seg in zip(xs, ys, zs): 
            if seg not in added_segments: 
                ax.plot(*seg, color="black")
                added_segments.append(seg)
    
    ax.set_xlabel(r"x", labelpad=10)
    ax.set_ylabel(r"y", labelpad=10)
    ax.set_zlabel(r"z", labelpad=10)

    return fig 

def load_batch_results(directory):
    if directory[-1] != '/' : directory += '/'

    with open(directory + "parameters_sets.table", 'r') as fil: lines = fil.readlines()

    header = lines[0][:-1].split("\t")
    sims = []
    for line in lines[1:]:
        sl = [float(i) for i in line.split("\t")]
        dic = {}
        for key, val in zip(header, sl) : dic[key] = val
        
        sim_dir = directory + "sim{:d}".format(int(dic["sim_id"]))
        results = ResultsReader(sim_dir)

        sims.append((dic, results))
    return sims


