import argparse, os, subprocess
from . import _version
from _pycbg_definitions import BUILD_DOC_SCRIPT

def main():
  
    parser = argparse.ArgumentParser(prog ='pycbg',
                                     description ='Manage CG-Geo MPM simulations using PyCBG Python module.',
                                     argument_default=argparse.SUPPRESS)
  
    parser.add_argument('script', metavar='PYCBG_SCRIPT', type=str, nargs='?',
                        help='%(prog)s script to be run. By default, the following import lines are added at the top of the file: `from pycbg.preprocessing import *`, `from pycbg.postprocessing import *`. To deactivate this behaviour, use the -n (or --no-import) option.')
    
    parser.add_argument('-n', '--no-import', action='store_true', default=False, dest="import_pycbg",
                        help="deactivate automatic import of %(prog)s when running PYCBG_SCRIPT")
    
    parser.add_argument('-d', '--build-doc', type=str, nargs='?', dest="build_doc",
                        help="build %(prog)s's documentation in the directory whose path is specified, the path being relative to the current working directory. If the directory already exists, it is removed without prompt before building the doc. If no path is specified, the directory will be named `pycbg_doc` and located in the current directory. If -d and PYCBG_SCRIPT are specified, the documentation is build before running the script.")

    parser.add_argument('-i', '--interactive', action='store_true', default=False, dest="interactive",
                        help="run the script in an interactive IPython session")

    parser.add_argument('-v', '--version', action='version', version=_version.get_versions()['version'],
                        help="print %(prog)s version")
  
    args = parser.parse_args()
  
    if hasattr(args, "build_doc"):
        directory = "pycbg_doc" if args.build_doc is None else args.build_doc
        subprocess.check_call([BUILD_DOC_SCRIPT, directory])

    if hasattr(args, "script"):
        with open(args.script, 'r') as fil: lines = fil.readlines()
        if args.import_pycbg: str_script = ""
        else: str_script = "from pycbg.preprocessing import *\nfrom pycbg.postprocessing import *\n"
        
        for line in lines: str_script += line

        exec(str_script, globals())
    globals().update(locals())

    if args.interactive: 
        from IPython.terminal.embed import InteractiveShellEmbed
        ipshell = InteractiveShellEmbed()

        if not hasattr(args, "script"): exec("from pycbg.preprocessing import *\nfrom pycbg.postprocessing import *", globals())
        globals().update(locals())

        ipshell()


        
