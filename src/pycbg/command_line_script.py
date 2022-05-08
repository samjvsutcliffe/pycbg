import argparse, os, subprocess, sys
from . import _version
from _pycbg_definitions import BUILD_DOC_SCRIPT

def main():
  
    parser = argparse.ArgumentParser(prog='pycbg',
                                     description='Manage CB-Geo MPM simulations using PyCBG Python module',
                                     argument_default=argparse.SUPPRESS)
  
    parser.add_argument('-v', '--version', action='version', version=_version.get_versions()['version'],
                        help="print %(prog)s version")
    
    parser.add_argument('-p', '--pip-show', action='store_true', dest="pip_show",
                        help="alias for `pip show pycbg`")

    parser.add_argument('script', metavar='PYCBG_SCRIPT', type=str, nargs='?',
                        help='%(prog)s script to be run. By default, the following import lines are added at the top of the file: `from pycbg.preprocessing import *`, `from pycbg.postprocessing import *` and `from pycbg.MPMxDEM import *`. To deactivate this behaviour, use the -n (or --no-import) option')
    
    parser.add_argument('-i', '--interactive', action='store_true', default=False, dest="interactive",
                        help="run in an interactive IPython session. Using both the -i and -n options simply creates a IPython interactive session")

    parser.add_argument('-n', '--no-import', action='store_true', default=False, dest="import_pycbg",
                        help="deactivates automatic import of %(prog)s")

    parser.add_argument('-d', '--build-doc', metavar="BUILD_DIR", type=str, nargs='?', dest="build_doc",
                        help="build %(prog)s's documentation in BUILD_DIR, its path being relative to the current working directory. If the directory already exists, it is removed without prompt before building the doc. If BUILD_DIR isn't specified, it will be set to `${PWD}/pycbg_doc`. If -d and PYCBG_SCRIPT are specified, the documentation is build before running the script")

  
    args = parser.parse_args()
  
    if hasattr(args, "pip_show"):
        pip_show_output = os.popen("python3 -m pip show pycbg")
        print(pip_show_output.read())

    if hasattr(args, "build_doc"):
        directory = "pycbg_doc" if args.build_doc is None else args.build_doc
        subprocess.check_call([BUILD_DOC_SCRIPT, directory])

    if hasattr(args, "script"):
        with open(args.script, 'r') as fil: lines = fil.readlines()
        if args.import_pycbg: str_script = ""
        else: str_script = "from pycbg.preprocessing import *\nfrom pycbg.postprocessing import *\nfrom pycbg.MPMxDEM import *\n"
        
        for line in lines: str_script += line

        exec(str_script, globals())
    globals().update(locals())

    if args.interactive or len(sys.argv) <= 1 or (len(sys.argv)==2 and sys.argv[1]=="-n"):
        print("Welcome to PyCBG \033[0;32m{:}\033[0m !".format(_version.get_versions()['version'])) 
        print("PyCBG's documentation can be accessed either locally by building it with the \033[0;35mpycbg -d\033[0m command, or at \033[0;36mhttps://pycbg.readthedocs.io/en/latest/\033[0m\n")
        from IPython.terminal.embed import InteractiveShellEmbed
        ipshell = InteractiveShellEmbed()

        if not args.import_pycbg and not hasattr(args, "script"): 
            exec("from pycbg.preprocessing import *\nfrom pycbg.postprocessing import *\nfrom pycbg.MPMxDEM import *", globals())
            print("Everything has been imported from \033[0;33mpycbg.preprocessing\033[0m, \033[0;33mpycbg.postprocessing\033[0m and \033[0;33mpycbg.MPMxDEM\033[0m\n")
        globals().update(locals())

        ipshell()