"""
Type checks all .py files inside src/phylojunction/.
Must be called from PhyloJunction/

Usage: python3 typecheck.py
"""

import pkgutil
import ast
import os
import subprocess
import phylojunction

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"

module_pkgs_or_mod_name_list = list()
module_mod_name_list_mypy_args = dict()

for _, mod_name, _ in pkgutil.walk_packages(path=phylojunction.__path__, prefix=phylojunction.__name__ + '.', onerror=lambda x: None):
    parsed_pkg_or_mod_name = "src/" + mod_name.replace(".", "/")
    module_pkgs_or_mod_name_list.append(parsed_pkg_or_mod_name)
    potential_mod_name = parsed_pkg_or_mod_name + ".py"
    
    if os.path.isfile(potential_mod_name):
        module_mod_name_list_mypy_args[potential_mod_name] = ""

imported_modules = set()
  
def visit_Import(node):
    for name in node.names:
        # imported_modules.add(name.name.split(".")[0])
        imported_modules.update(name.name.split("."))

def visit_ImportFrom(node):
    # if node.module is missing it's a "from . import ..." statement
    # if level > 0 it's a "from .submodule import ..." statement
    if node.module is not None and node.level == 0:
        imported_modules.add(node.module.split(".")[0])

node_iter = ast.NodeVisitor()
node_iter.visit_Import = visit_Import
node_iter.visit_ImportFrom = visit_ImportFrom

# if calling script from different paths, sometimes it is necessary to use these arguments to guide package discovery
# for mod_name in module_mod_name_list_mypy_args:
#     imported_modules = set()

#     with open(mod_name, "r") as f:
#         node_iter.visit(ast.parse(f.read()))

#     deepest_pkg_name = os.path.dirname(mod_name).split("/")[-1]

    # print("\nReading " + mod_name + ":")
    # print("    imported modules: " + " ".join(imported_modules))
    # print("    deepest_pkg_name: " + deepest_pkg_name)

    # if deepest_pkg_name in imported_modules:
    #     module_mod_name_list_mypy_args[mod_name] = ["--namespace-packages", "--explicit-package-bases"]
    

def call_mypy(module_args_dict):
    for module_name, arg_list in module_args_dict.items():
        # module_name_within_pj = module_name.replace("src/phylojunction/", "")
        module_name_within_pj = module_name
        cmd_list = ["mypy"]
        cmd_list.append(module_name_within_pj)

        if arg_list: cmd_list.extend(arg_list)

        # see command
        print(" ".join(cmd_list))

        print("Type checking " + module_name_within_pj)
        p = subprocess.Popen(cmd_list, cwd="./")
        p.wait()

call_mypy(module_mod_name_list_mypy_args)
