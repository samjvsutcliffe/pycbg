import os, pickle, io
import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d import Axes3D

import multiprocessing as mp

from PIL import Image

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
        Particles' stresses for every saved steps. Noting npart the number of particles in the simulations at the ith step, the shape of ``pstresses[i]`` is ``(npart,6)``. The columns respectively correspond to `xx`, `yy`, `zz`, `xy`, `yz` and `xz` components.
    pstrains : list of numpy arrays
        Particles' strains for every saved steps. Noting npart the number of particles in the simulations at the ith step, the shape of ``pstrains[i]`` is ``(npart,6)``. The columns respectively correspond to `xx`, `yy`, `zz`, `xy`, `yz` and `xz` components.
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
        files = [f for f in files if f[-4:]=='.csv']

        raw_data, steps = [], []
        for f in files: 
            raw_data.append(pd.read_csv(self.data_dir + f, sep="\t", header=0))
            steps.append(int(f[9:-4]))

        sort_ind = np.argsort(steps)
        self.raw_data, self.steps = [], []
        for i in sort_ind: 
            self.raw_data.append(raw_data[i])
            self.steps.append(steps[i])

        self.ppositions, self.pvelocities = [], []
        self.pstresses, self.pstrains, self.ppressures = [], [], []
        self.pmasses, self.pvolumes, self.pmaterials = [], [], []
        n_mp_init = len(self.raw_data[0]["id"])
        for df in self.raw_data :
            ids = list(df["id"])

            # Positions
            ppos = np.full([n_mp_init, 3], np.nan)
            ppos_df = np.array([df[key] for key in ["coord_x", "coord_y", "coord_z"]]).T
            for i, p_id in enumerate(ids): ppos[p_id,:] = ppos_df[i, :]
            self.ppositions.append(ppos)

            # Velocities
            pvel = np.full([n_mp_init, 3], np.nan)
            pvel_df = np.array([df[key] for key in ["velocity_x", "velocity_y", "velocity_z"]]).T
            for i, p_id in enumerate(ids): pvel[p_id,:] = pvel_df[i, :]
            self.pvelocities.append(pvel)

            # Stresses
            psig = np.full([n_mp_init, 6], np.nan)
            psig_df = np.array([df[key] for key in ["stress_xx", "stress_yy", "stress_zz", "tau_xy", "tau_yz", "tau_xz"]]).T
            for i, p_id in enumerate(ids): psig[p_id,:] = psig_df[i, :]
            self.pstresses.append(psig)
            
            # Strains
            peps = np.full([n_mp_init, 6], np.nan)
            peps_df = np.array([df[key] for key in ["strain_xx", "strain_yy", "strain_zz", "gamma_xy", "gamma_yz", "gamma_xz"]]).T
            for i, p_id in enumerate(ids): peps[p_id,:] = peps_df[i, :]
            self.pstrains.append(peps)

            # Pressures
            ppre = np.full([n_mp_init], np.nan)
            ppre_df = df['pressure'].values
            for i, p_id in enumerate(ids): ppre[p_id] = ppre_df[i]
            self.ppressures.append(ppre)

            # Materials
            pmat = np.full([n_mp_init], np.nan)
            pmat_df = df['material_id'].values
            for i, p_id in enumerate(ids): pmat[p_id] = pmat_df[i]
            self.pmaterials.append(pmat)

            # Volumes
            pvol = np.full([n_mp_init], np.nan)
            pvol_df = df['volume'].values
            for i, p_id in enumerate(ids): pvol[p_id] = pvol_df[i]
            self.pvolumes.append(pvol)

            # Masses
            pmas = np.full([n_mp_init], np.nan)
            pmas_df = df['mass'].values
            for i, p_id in enumerate(ids): pmas[p_id] = pmas_df[i]
            self.pmasses.append(pmas)
    
__ind_sgmts = [[0,1], [1,2], [2,3], [3,0],
               [4,5], [5,6], [6,7], [7,4],
               [0,4], [1,5], [2,6], [3,7]]

__ind_gimp_sgmts = [[8,19], [9, 18], [24,31], [25,30], 
                    [17,21], [16,20], [12,22], [13,23],
                    [14,10], [15,11], [28,26], [29,27],
                    [16,12], [17,13], [20,22], [21,23],
                    [8,24], [9,25], [19,31], [18,30],
                    [14,28], [15,29], [10,26], [11,27],
                    [32,52], [35,53], [34,54], [33,55],
                    [39,50], [36,49], [37,48], [48,51],
                    [47,41], [44,40], [46,42], [45,43]]

def plot_mesh(mesh, fig=None, ax=None):
    if fig is None and ax is None: 
        fig = plt.figure(1, figsize=(15,10))
        ax = fig.add_subplot(111, projection='3d')

    added_segments, added_dotted_segments = [], []
    for cell in mesh.cells:
        nodes = [mesh.nodes[i_node] for i_node in cell]
        # for i, node in enumerate(nodes): ax.text(*node, str(i)) # to print node numbers
        xs = [[nodes[i][0], nodes[j][0]] for i,j in __ind_sgmts]
        ys = [[nodes[i][1], nodes[j][1]] for i,j in __ind_sgmts]
        zs = [[nodes[i][2], nodes[j][2]] for i,j in __ind_sgmts]

        for seg in zip(xs, ys, zs): 
            if seg not in added_segments: 
                ax.plot(*seg, color="black")
                added_segments.append(seg)

        if not mesh.cell_type=="ED3H64G": continue

        xds = [[nodes[i][0], nodes[j][0]] for i,j in __ind_gimp_sgmts]
        yds = [[nodes[i][1], nodes[j][1]] for i,j in __ind_gimp_sgmts]
        zds = [[nodes[i][2], nodes[j][2]] for i,j in __ind_gimp_sgmts]

        for seg in zip(xds, yds, zds): 
            if seg not in added_dotted_segments: 
                ax.plot(*seg, color="gray", linestyle="dotted")
                added_dotted_segments.append(seg)

    
    ax.set_xlabel(r"x", labelpad=10)
    ax.set_ylabel(r"y", labelpad=10)
    ax.set_zlabel(r"z", labelpad=10)

    return fig 

def load_batch_results(directory, fail_flag=True):
    if directory[-1] != '/' : directory += '/'

    with open(directory + "parameters_sets.table", 'r') as fil: lines = fil.readlines()

    header = lines[0][:-1].split("\t")
    sims = []
    for line in lines[1:]:
        sl = [float(i) for i in line.split("\t")]
        dic = {}
        for key, val in zip(header, sl) : dic[key] = val
        
        sim_dir = directory + "sim{:d}".format(int(dic["sim_id"]))
        if fail_flag: results = ResultsReader(sim_dir)
        else:
            try: results = ResultsReader(sim_dir)
            except FileNotFoundError: break

        sims.append((dic, results))
    return sims

def load_results(sim_results):
    """Manage basic plotting of results, mostly material points' positions during all simulation.

    Parameters
    ----------
    results: :class:`~pycbg.postprocessing.ResultsReader` object
            Results to plot graphs of.
    """
    global results, n_mp, n_saved_steps, sim, lcs, mesh_lims, cell_ids, pos_ss
    results = sim_results
    n_mp = results.ppositions[0].shape[0]
    n_saved_steps = len(results.ppositions)
    sim = results.load_simulation()
    lcs = [l/nc for l,nc in zip([sim.mesh.l0, sim.mesh.l1, sim.mesh.l2], [sim.mesh.nc0, sim.mesh.nc1, sim.mesh.nc2])]
    mesh_lims = [(o, o+l) for o,l in zip(sim.mesh.origin, [sim.mesh.l0, sim.mesh.l1, sim.mesh.l2])]
    cell_ids = [[] for i_step in range(len(results.steps))]

    pos_ss = [[[] for i_p in range(n_mp)], [[] for i_p in range(n_mp)], [[] for i_p in range(n_mp)]]
    for i_p in range(n_mp):
        for poss in results.ppositions: 
            pos_ss[0][i_p].append(poss[i_p, 0])
            pos_ss[1][i_p].append(poss[i_p, 1])
            pos_ss[2][i_p].append(poss[i_p, 2])

    for istep, data in enumerate(results.raw_data): 
        for i_p in range(n_mp): cell_ids[istep].append(data["cell_id"].values[i_p])


def plot_positions(step_index, plane, colored_by=None, color_label=None, color_scale="linear", patch_kwargs={"facecolor":"grey", "alpha":.1, "edgecolor":"grey"}, get_cmap_args={"name":"autumn_r", "lut":1e3}, colorbar_args={"orientation":"horizontal"}):
    """Plot material points' projection of their positions on the plane defined by `plane`. Material points can be colored by all variables saved in the csv file. Additional matplotlib commands can be executed through the `tweak_plt` function.

    Parameters
    ----------
    step_index: int
        Index of the step to plot.
    plane: tuple of ints
        Projection axes' indexes. For instance `plane=(0,1)` means that the horizontal axis will be the axis `0` and the vertical axis will be the axis `1`.
    colored_by: str, list of matplotlib colors or None
        Color of all material points. If a string is given, the color will represent the variable in the csv files whose header is this string (e.g. `displacement_x`, `stress_yy`, `svars_5`, ...). If a list of matplotlib colors is given, it should contain as much elements as material points and the index of each color corresponds to the material point's ID. If `None` is given, all material points are black.
    color_label: str or None
        Label for the colorbar. If None is given and `colored_by` is a string, `color_label` is set to `colored_by`. If None is given and `colored_by` is a list, no colorbar is plotted. 
    color_scale: {"linear", "log"}
        Scale to use for colors.
    patch_kwargs: dict
        Keyword arguments for matplotlib.collections.PatchCollection, can be used to change how mesh's cells are displayed on the graph (for instance, if there is many material point per cell one might want to set a low "alpha").
    get_cmap_args: dict
        Arguments for matplotlib's get_cmap function (i.e. `name` and `lut`). 
    colorbar_args: dict
        Arguments for matplotlib.colorbar.ColorbarBase, can be used to set a specific colormap range for several graphs.

    Returns
    -------
    matplotlib.figure.Figure object
        The matplotlib figure object.
    numpy array of matplotlib.axes._subplots.AxesSubplot objects
        The axes of the matplotlib figure. The first element corresponds to the figure itself, the second corresponds to the colorbar if any.
    """
    # Useful variables
    plot_colorbar = (type(colored_by) == str) or (type(color_label) == str) or ("cmap" in colorbar_args)
    height_ratios = [10,1] if plot_colorbar else [1]
    ix, iy = plane
    x_ss, y_ss = pos_ss[ix], pos_ss[iy]
    lcx, lcy = lcs[ix], lcs[iy]
    mesh_cells, mesh_nodes = sim.mesh.cells, sim.mesh.nodes

    # If colored_by is a string
    if type(colored_by)==str: 
        # Check color_label value
        if color_label is None: color_label=colored_by

        # Get colors
        values = results.raw_data[step_index][colored_by].values
        cmap = mpl.cm.get_cmap(**get_cmap_args)
        norm = _get_norm(values, color_scale)
        if not "cmap" in colorbar_args: colorbar_args["cmap"] = cmap
        if not "norm" in colorbar_args: colorbar_args["norm"] = norm
        colored_by = [cmap(norm(values[mp_id])) for mp_id in range(n_mp)]
    
    # Each material point is black by default
    elif colored_by is None: colored_by = ["black" for mp_id in range(n_mp)]

    # Create figure object
    fig, axes = plt.subplots(nrows=2 if plot_colorbar else 1, ncols=1, figsize=(17,10), gridspec_kw={"height_ratios": height_ratios})

    # If there is no colorbar, reformat axes variable
    if not plot_colorbar: axes = np.array([axes])

    # Color cells depending on how much material points are inside
    cells = []
    for cell_id in cell_ids[step_index]:
        if np.isnan(cell_id): continue
        x,y = np.array([(mesh_nodes[i][ix], mesh_nodes[i][iy]) for i in mesh_cells[cell_id]]).min(0)
        cells.append(mpl.patches.Rectangle((x,y), lcx, lcy))
    pc = mpl.collections.PatchCollection(cells, **patch_kwargs)
    axes[0].add_collection(pc)

    # Plot material points
    for mp_id in range(n_mp): axes[0].scatter(x_ss[mp_id][step_index], y_ss[mp_id][step_index], color=colored_by[mp_id])

    # Plot colorbar if necessary
    if plot_colorbar:
        cb = mpl.colorbar.ColorbarBase(axes[1], **colorbar_args)
        cb.set_label(color_label)

    # Improve graph
        # Labels
    all_labels = [r"$x$ (m)", r"$y$ (m)", r"$z$ (m)"]
    x_label, y_label = all_labels[ix], all_labels[iy]
    axes[0].set_xlabel(x_label)
    axes[0].set_ylabel(y_label)

        # Limits
    axes[0].set_xlim(mesh_lims[ix])
    axes[0].set_ylim(mesh_lims[iy])

        # Aspect ratio
    axes[0].set_aspect("equal")

    return fig, axes

def plot_all_positions(plane, colored_by=None, n_cores=1, **plot_positions_kwarg):
    """Plot material points' projection of their positions on the plane defined by `plane` for all time steps. Takes as keyword argument all arguments of the `plot_positions` method.

    Parameters
    ----------
    plane: tuple of ints
        Projection axes' indexes. For instance `plane=(0,1)` means that the horizontal axis will be the axis `0` and the vertical axis will be the axis `1`.
    colored_by: list of list, or str, or list, or None
        If `colored_by` is a list of list, each of its element are passed to `plot_positions` as the `colored_by` argument. If it is a simple list, a string, or `None`, it is directly passed to `plot_positions` as the `colored_by` argument.
    n_cores: int
        Number of cores to use for parallel plotting.
    **plot_positions_kwarg: `plot_positions` keyword arguments
        Arguments to be passed to `plot_positions`.

    Returns
    -------
    list of tuples
        Each element is a couple of figure and axes, as returned by the `plot_positions` method.
    """
    # Manage colored_by input
    if colored_by is None: set_colors = False
    elif type(colored_by)==str or type(colored_by[0])==list: set_colors = True
    else: set_colors = False
    if not set_colors: colored_by_s = [colored_by for i in range(n_saved_steps)]
    else: 
        # Set color_scale if user didn't
        if not "color_scale" in plot_positions_kwarg: plot_positions_kwarg["color_scale"] = "linear"

        # Set get_cmap_args if user didn't
        if not "get_cmap_args" in plot_positions_kwarg: plot_positions_kwarg["get_cmap_args"] = {"name":"autumn_r", "lut":1e3}

        # Get values
        if type(colored_by)==str: values = np.array([results.raw_data[i][colored_by].values for i in range(n_saved_steps)])
        else: values = np.array([colored_by[i] for i in range(n_saved_steps)])
        
        # Get colors
        cmap = mpl.cm.get_cmap(**plot_positions_kwarg["get_cmap_args"])
        norm = _get_norm(values, plot_positions_kwarg["color_scale"])
        colored_by_s = [[cmap(norm(values[i, mp_id])) for mp_id in range(n_mp)] for i in range(n_saved_steps)]

        # Manage colorbar_args argument
        if not "colorbar_args" in plot_positions_kwarg: plot_positions_kwarg["colorbar_args"] = {"cmap":cmap, "norm":norm, "orientation":"horizontal"}
        else:
            if "cmap" not in plot_positions_kwarg["colorbar_args"]: plot_positions_kwarg["colorbar_args"]["cmap"] = cmap
            if "norm" not in plot_positions_kwarg["colorbar_args"]: plot_positions_kwarg["colorbar_args"]["norm"] = norm
        
    args = [((i, plane, colored_by), plot_positions_kwarg) for i, colored_by in enumerate(colored_by_s)]

    # Plot figures
    if n_cores>1:
        p = mp.get_context('fork').Pool(processes=n_cores)
        figs_and_axes = p.map(_plot_positions_expend_args, args)
    else: figs_and_axes = list(map(_plot_positions_expend_args, args))

    return figs_and_axes

def _plot_positions_expend_args(args):
    f_args, f_kwargs = args
    return plot_positions(*f_args, **f_kwargs)

def _get_norm(values, color_scale="linear"):
    # Check color_scale value
    if color_scale=="log": normalizer = mpl.colors.LogNorm
    elif color_scale=="linear": normalizer = mpl.colors.Normalize
    else: raise RuntimeError("color_scale must be set to 'linear' or 'log'")
    
    # Get colors
    vmin, vmax = values.min(), values.max()
    if vmin==vmax: 
        vmax *= 1.01
        vmin *= .99
    norm = normalizer(vmin=vmin, vmax=vmax)

    return norm

def make_gif(figures, filename="video.gif", max_size=(1000, 1000), pil_save_kwargs={"quality":95, "duration":.1, "optimize":True, "loop":0, "save_all":True}, n_cores=1):
    """Make a gif from all figures in `figures` in the specified order.

    Parameters
    ----------
    figures: list of matplotlib.figure.Figure objects
        Figures to be used to make the gif. Their order in the list is preserved in the gif.
    filename: str
        Path to the gif file, relative to the current working directory. 
    max_size: tuple of ints
        Maximum number of pixels in each direction.
    pil_save_kwargs: dict
        Keyword arguments passed to PIL's `Image.save` fuction, see `https://pillow.readthedocs.io/en/stable/handbook/image-file-formats.html#gif` for more details.
    """
    def thumbnail_wrapper(im):
        nonlocal max_size
        im.thumbnail(max_size, Image.ANTIALIAS)

    # Turn figures into PIL images
    if n_cores>1:
        p = mp.get_context('fork').Pool(processes=n_cores)
        all_images = p.map(_convert_mpl_to_pil, figures)
        p.map(thumbnail_wrapper, all_images)
    else: 
        all_images = list(map(_convert_mpl_to_pil, figures))
        map(thumbnail_wrapper, all_images)
    
    all_images = [_convert_mpl_to_pil(fig) for fig in figures]
    for im in all_images: im.thumbnail(max_size, Image.ANTIALIAS)
    img, *imgs = all_images

    img.save(fp=filename, format='GIF', append_images=imgs, **pil_save_kwargs)

def _convert_mpl_to_pil(fig):
    buf = io.BytesIO()
    fig.savefig(buf, bbox_inches="tight")
    buf.seek(0)
    img = Image.open(buf)
    return img
