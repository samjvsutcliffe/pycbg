import argparse, os, subprocess, sys
from . import _version
from _pycbg_definitions import BUILD_DOC_SCRIPT
import pycbg.postprocessing as utl
import matplotlib as mpl

mpl.rcParams['font.size'] = 25 # Setting the size of the fonts in the graphs

def main():
  
    parser = argparse.ArgumentParser(prog='pycbg-gif',
                                     description="Make a gif of material points' positions for the given simulation",
                                     argument_default=argparse.SUPPRESS)
  
    parser.add_argument('sim', metavar='SIMULATION_DIR', type=str, nargs='?',
                        help="Path to the simulation's directory to be plotted")
    
    parser.add_argument('plane', metavar='PROJECTION_PLANE', type=str, nargs='?',
                        help="Plane on which to project material points' positions. For instance, `1,2` will project positions on the (y,z) plane.")
    
    parser.add_argument('-c', '--colored-by', metavar="COLOR_VAR", type=str, nargs='?', dest="colored_by",
                        help="Variable to use for material points' coloring. All variables in the csv results files are available, their name being their column's header. Can be a comma-separated list of several variables. If not provided, all material points are black.")
    
    parser.add_argument('-l', '--color-label', metavar="COLOR_LABEL", type=str, nargs='?', dest="color_label",
                        help="Label for colorbar, can be a latex expression. If a list of variables was specified with the c option, then it should be a comma-separated list of labels. Default to variables' names.")
    
    parser.add_argument('-o', '--output-file', metavar="OUTPUT_FILE", type=str, nargs='?', dest="file_prefix",
                        help="Path of the output file, without the '.gif' extension. If a list of variables was specified with the -c option, all file are suffixed with the name of the coloring variable. Default is 'video'.")

    parser.add_argument('-j', '--parallel-jobs', metavar="N_JOBS", type=str, nargs='?', dest="n_cores",
                        help="Number of cores to use. Default to 1.")

  
    args = parser.parse_args()

    results = utl.ResultsReader(args.sim)
    utl.load_results(results)

    plane = [int(axis) for axis in args.plane.split("(")[-1].split(")")[0].split(",")]
  
    use_default_color = False
    if not hasattr(args, "colored_by"): 
        colored_by_s = [["black" for mp_id in range(utl.n_mp)]]
        use_default_color = True
    else: 
        if not "," in args.colored_by: colored_by_s = [args.colored_by]
        else: colored_by_s = args.colored_by.split(",")
        if hasattr(args, "color_label"):
            if not "," in args.color_label: color_label_s = [fr"{args.color_label}"]
            else: color_label_s = [fr"{lab}" for lab in args.color_label.split(",")]
    
    if not "color_label_s" in locals(): color_label_s = [None]*len(colored_by_s)
    
    if not hasattr(args, "file_prefix"): file_prefix = "video"
    else: file_prefix = args.file_prefix

    if not hasattr(args, "n_cores"): n_cores = 1
    else: n_cores = int(args.n_cores)

    for colored_by, color_label in zip(colored_by_s, color_label_s):
        figs = [res[0] for res in utl.plot_all_positions(plane, colored_by=colored_by, color_label=color_label, n_cores=n_cores)]
        if type(colored_by)==str: suffix = "_" + colored_by
        elif use_default_color: suffix = ""
        else: raise RuntimeError("Can't make color label")
        
        filename = file_prefix + suffix + ".gif"
        utl.make_gif(figs, n_cores=n_cores, filename=filename)

