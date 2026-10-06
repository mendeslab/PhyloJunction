# Notes for developers

## Making packages and modules discoverable on VSCode

The first thing you must do before you can start coding within the `PhyloJunction` framework
is to set it up within an integrated development environment (IDE).
The IDEs of choice here are VSCode or PyCharm, which arevery popular, clean and fast IDEs.

### VS Code

The first step is to go under "File" > "Open Folder..." and choose `PhyloJunction/`, the git repository root folder.
If you click around for a bit, you will notice that the module files (`.py` files) are inside
the `src/` directory, nested within different packages (e.g., `phylojunction/`, `phylojunction/distribution/`).
We must point our python interpreter to the topmost package, which is called `phylojunction/` so that all modules can be found.
Here, we must create a `launch.json` file, which VSCode can help us with.

Click the button on the left that has a bug with a triangle ("Run and Debug"), and then "Create launch.json file", choosing module.
This should create a that looks like this:

```
{
    // Use IntelliSense to learn about possible attributes.
    // Hover to view descriptions of existing attributes.
    // For more information, visit: https://go.microsoft.com/fwlink/?linkid=830387
    "version": "0.2.0",
    "configurations": [
        {
            "name": "Python: Current File",
            "type": "python",
            "request": "launch",
            "program": "${file}",
            "console": "integratedTerminal",
            "justMyCode": true
        }
    ]
}
```

In order to tell the interpreter to look inside `src/phylojunction/` for modules, add the following
line below `"justMyCode"`:

```
"env": {"PYTHONPATH": "${workspaceRoot}/src/phylojunction/"}
```

After doing this, you might notice if you choose one of the `.py` files that the import statements at the top of the file will be marked as problematic.
This is because VSCode itself cannot discover the imported packages and modules.
We must fix that, which we do by creating a `settings.json` file inside hidden folder `.vscode`.
Enter shortcut "CMD + Shift + P", type `settings.json` and choose "Open Workspace Settings (JSON)".

Now inside `settings.json`, add this:

```
{
    "python.analysis.extraPaths": [ "${workspaceFolder}/src/phylojunction/" ]
}
```

**Configuring unit tests (using unittest) on VSCode**

If you have an Apple machine, the shortcut on VSCode is "CMD + Shift + P", "Python: Configure Tests".
This will add to your `settings.json` file inside hidden folder `.vscode`.
Choose `.` as the directory containing tests, and `test_*` as the prefix for tests.

Your `settings.json` will look like this:

```
{
    "python.analysis.extraPaths": [ "${workspaceFolder}/src/phylojunction/" ],
    "python.testing.unittestArgs": [
        "-v",
        "-s",
        ".",
        "-p",
        "test_*.py"
    ],
    "python.testing.pytestEnabled": false,
    "python.testing.unittestEnabled": true
}
```

If VSCode refuses to find the tests, one thing that might be set up wrong is the Python interpreter VSCode is using.
Unit tests usually have many dependencies -- if those dependencies are not found because the wrong interpreter is used (the native one that comes with the machine, say, instead of an interpreter installed with Homebrew), then the tests will not be found.

### PyCharm

Once you Open Project, you must tell PyCharm where the source files are: PyCharm > Settings > Project: PhyloJunction > Project Structure, the click on `src/` and click on the folder "Sources" icon.

## Type checking

`mypy`'s documentation suggests calling mypy on entire packages, e.g.,

```
cd src/phylojunction/
mypy data/
mypy calculation/
```

But it's also possible to call it on individual files:

```
mypy data/tree.py
mypy utility/helper_functions.py --namespace-packages --explicit-package-bases
mypy distribution/dn_discrete_sse.py --namespace-packages --explicit-package-bases
```

`mypy` will verify type hinting in source files against mypy stubs (files ending with `.pyi`).
Stubs must be kept up-to-date by developers so `mypy` can find incongruencies.

Sometimes I have noticed that calling `mypy` on individual modules that depend on other modules within the same package fails (e.g., `helper_functions.py` is in the same package as, and depends on `exception_classes.py`).
Fixing it can be done by guiding `mypy`'s import discovery; 
This is what the flags `--namespace-packages` together with `--explicit-package-bases` do above.
They are necessary even when an `__init__.py` is present in the same folder as the modules being type checked.

PJ has a script that calls `mypy` automatically on everything:

```
cd PhyloJunction/
python3 typecheck.py
```

This script should be called and be successful every time new code is written, before committing and pushing.

## Style convention: PEP8

In order to check for PEP8 style convention rules, you need the `pycodestyle` (previously called `pep8`) python package installed.
Checking one file can be done with:

```
pycodestyle [path to .py]
```

It is recommended to do that from within VS Code's terminal, for example, that way it is easy to jump immediately to the violating line.

A couple of rules are always broken because of type hinting, and because of an automatic newline character probably added by VS Code.
These are rules E701 and W391.
Sometimes, excessive indentation (which should be avoided anyway) causes lines to be > 79 characters.
This violates rule E501.
Another pair of problematic rules (W503 and W504) we can ignore is is "line break after binary operator" (it pops up for no obvious reason sometimes).
We can also ignore indentation in comments (E116) and no newline at the end of file (W292)
We can ask pycodestyle to ignore them:

```
pycodestyle --ignore=E116,E501,W503,W504,E701,W391 [path to .py]
```


## Editable pip install for developers

Building and installing PhyloJunction is done using `setuptools`, which is specified in `pyproject.toml`.
Instructions for `setuptools` are provided in `setup.cfg`.

Start by cloning the PhyloJunction repository, and going to its root (i.e., the `PhyloJunction` folder).

If you are on an Apple machine that obstrusively curtails file system access, run the following command:

```
cd PhyloJunction/
python3.11 -m pip install --user -e .
```

The `--user` flag is the way around said restrictions.
This will tell `pip` to place the entry point executables (`pjcli` and `pjgui`) in  `/Users/user_name/Library/Python/3.9/bin/` if you are using Python 3.9, say.
Remember to add this path to the PATH system variable if you want to call the executables from anywhere on your file system.
Another option is to run the following command:

```
python3.11 -m pip install --prefix ~/.local -e .
```

Which creates (or writes inside) directories `bin/` (placing the executables therein) and `lib/python3.X/site-packages/` (placing the egg-link therein) inside `~/.local`.
This is an attractive option if you normally already have `~/.local/bin` as part of your PATH variable.

Add a shell variable for PYTHONPATH to your bash profile, for example
```
PYTHONPATH="/Users/[username]/projects/PhyloJunction/src"
```
(Note: the file path associated with this `PYTHONPATH` variable ends with the `src` directory, not `src/phylojunction` -- i.e., irrespective of "env" variable in your `launch.json` file in VS Code).

You can test things worked by trying, from any directory:

```
pjgui
```

Important notes:
    (1) On Apple machines, sometimes when Homebrew updates itself and starts updating python3.X.Y, say, and depending on your computer architecture (e.g., M1 or M2 chips). Make sure you update your interpreter (command + shift + P, Python: Select Interpreter) to the newest version (e.g., `/opt/homebrew/Cellar/python@3.9/3.9.16/Frameworks/Python.framework/Versions/3.9/bin/python3`)

    (2) On Linux machines, depending on the python version (e.g., 3.9.7), the PySimpleGUI GUI will have its Helvetica fonts replaced by a Courier-like font and be all messed up. It might have to do with running the GUI within a conda environment though. The way around is to reinstall PJ (with `pip` as described above) with another Python version on the Terminal, by quitting the conda environment (`conda deactivate`) if necessary

### Known issues

* Installing `scipy` for python3.10 under M1 architecture is complicated (it seems issues were fixed, because a pip install for 3.11 worked for scipy)

## Releasing PJ to TestPyPI

First, create an account on https://test.pypi.org/.
We will use TestPyPI to test a package before releasing it for real with PyPI.

Then start by logging in to TestPyPI, and go to Account Settings.
Then set up 2FA with an app like Authy, for example.
With the 2FA set up, get an API token -- this will save time every time the package is pushed to TestPyPI.
Then place the token inside $HOME/.pypirc, like so:

```
[testpypi]
  username = __token__
  password = <API token here>
```

Now from Phylojunction/ (make sure you have the `build` Python package installed):

```
python3.11 -m build --sdist
python3.11 -m build --wheel
```

These two commands will create a dist/ directory inside PhyloJunction/, and place a .tar.gz (containing all source .py/.pyi files, license, readme, setup files) and a wheel file inside.
You can check the .tar.gz to see if all modules are in there.
If one of them is missing, make sure an `__init__.py` file is inside the module's directory.

Then to release to TestPyPI, do:

```
twine upload --repository testpypi dist/*
```

And to test the package, do:

```
python3.11 -m pip install --index-url https://test.pypi.org/simple/ phylojunction
```

### Next release

When a new release must be made, in this order:
1. Delete the files inside `dist/`,
2. Change the version of the package in `setup.cfg`,
3. Dun the build commands above again,
4. Push to TestPyPI.

## Releasing PJ to PyPI

To release to PyPI, we just change the twine command in the previous section to

```
twine upload dist/*
```

To see where pip installs things:

```
python3 -m pip list -v
``` 

The PyPI link to PJ is: https://pypi.org/project/PhyloJunction/

## Icon for launching PJ 

On Linux, there are two ways to have icons for PJ:

(1) You add the icon to the dash (command key) and/or the taskbar (by making the icon a favorite)

You do so by creating a file `/home/foo/.local/share/applications/.desktop`, and putting the following inside:

```
[Desktop Entry]
Version=0.0.2
Name=PhyloJunction
Comment=An environment for prototyping, teaching and learning with evolutionary models
Exec=nohup pjgui &
Icon=/path/to/PhyloJunction/some_icon.svg
Terminal=false
Type=Application
Categories=Application;
```

When you click on the dash or taskbar icon, even with `Terminal=false`, a Terminal window will still be spawned.
The `nohup pjgui &` prevents that.

(2) You add an icon file to the desktop

You do so by creating a file `/home/foo/Desktop/PhyloJunction.desktop', and again, put the following inside:

```
[Desktop Entry]
Type=Application
Name=PhyloJunction
Comment=An environment for prototyping, teaching and learning with evolutionary models
Exec=pjgui
Icon=/path/to/PhyloJunction/some_icon.svg
Terminal=false
Categories=Application;
```

Then, you have to right-click the icon that will appear on your Desktop, and make sure you click "Allow launching".

Important notes:
    (1) It can happen that you mess up something like a module import, but the GUI still functions on the Terminal (`pjgui`) or inside VS Code (`python3 PhyloJunction/src/phylojunction/interface/pj_gui.py`).
    But then when you double-click its icon, nothing happens. 
    In order to see what the error is, you need to set `Terminal=true` inside the `.desktop` file, and then make sure that, in your Terminal app, when you go to "Preferences" > "Profiles", that option "Hold the terminal open" is chosen for "When command exits".
    After the GUI bombs, the Terminal will stay open with the error messages and you can then debug it.


On Mac OS X, PySimpleGUI provides the `use_custom_titlebar=True` (which automatically means `no_titlebar=True`). 
These allow `titlebar_icon="some_icon.png"`, which places a .png file on the titlebar, but simultaneously makes the window non-mini/maximizable on Mac OS X.
So PJ stays cross-platform with windows resizability, I won't be adding any titlebar icons.

## GUI development

PhyloJunction has two functional GUIs:

(1) with PySimpleGUI (`src/interface/pj_gui.py`), and
(2) with PySide6 (`src/pysidegui/pj_gui_pyside`).

The GUI built with PySimpleGUI was coded entirely by hand, and is out-of-date at this point.

The GUI built with PySide6 had its main structure coded by hand (`pysidegui/content_main_window.py`), but then individual pages built with Qt Creator / Qt Designer (an app for automating PySide6 GUI development).
These pages are saved in `pjguipages/pjgui_pages.ui`. From the repository root, regenerate
its Python module using the installed PySide6 tool (no PyQt import conversion is needed):

```sh
.venv/bin/pyside6-uic src/phylojunction/interface/pysidegui/pjguipages/pjgui_pages.ui -o src/phylojunction/interface/pysidegui/pjguipages/gui_pages.py
```

After generation, retain the resource registration import alongside the custom widget imports:

```python
from phylojunction.interface.pysidegui.images.icons import resources
```

Check that the draw and clear icons use `:/draw.svg` and `:/icon_clear.svg`, respectively.
The UI file supplies the custom widget module paths. Review generated changes separately from
layout changes when switching generator versions.

### Qt Creator / Qt Designer

This can be a finnicky program, so here are a few notes.

(1) Adding and deleting new pages

In order to add new pages to a QStackedWidget (the class of the main container), you need to right-click the QStackedWidget on the right menu, and do "Insert page".
You can then "Change page order" after the new page is created.

If you want to delete a page (e.g., one of the two pages that come automatically with QStackedWidget), you need to navigate to the page you want to delete.
To do that, you must right-click the QStackedWidget, and click "Next Page" until you reach the page to be deleted.
Then above the "Insert page" menu (from right-clicking QStackedWidget), there will be a menu saying "Page 2 out of 3", say. 
There you will find a "Delete" button.

(2) Adding a custom widget (that makes use of some .py code the developer wrote) to Qt Creator's .ui pages

If you have a custom widget (e.g., QPushButton) that you want to place on a page being designed with Qt Creator, you can start by placing a QWidget on the page.
Then you right-click the QWidget > "Promote to...". 
Here, you can add the name ("Promoted class name") of the custom Python class implementing the custom widget.
For "Header file", you want to put the entire Python path to the custom class, as it would appear in the imports of a Python script.
E.g., `phylojunction.interface.pysidegui.pjguiwidgets.file_defining_class`

pyuic5 or pyuic6 creates a .py from Qt Creator's ui, it will place import statements at the bottom of the file, so that your module can see the classes used by the custom widgets.

### Resources 

Certain types of data files like icon images for the GUI must be listed in a `.qrc` resource file, like so:

```
<!DOCTYPE RCC>
<RCC version="1.0">
    <qresource>
        <file>icon_clear.svg</file>
        <file>icon_clear_over.svg</file>
        <file>icon_clear_pressed.svg</file>
        <file>draw.svg</file>
        <file>draw_pressed.svg</file>
        <file>draw_over.svg</file>
    </qresource>
</RCC>
```

The above content is file `resource.qrc` inside `pysidegui/images/icons`. 
This resource file must then be converted into a `.py` file, like so:

```
pyside6-rcc resources.qrc -o resources.py
```

To make these images visible to the GUI, we need to import the resource file in two locations.
First, in `pj_buttons.py`, we need to set the paths to images like so:

```
clear_dag_normal = QIcon(QPixmap(":/icon_clear.svg"))
clear_dag_over = QIcon(QPixmap(":/icon_clear_over.svg"))
clear_dag_pressed = QIcon(QPixmap(":/icon_clear_pressed.svg"))
```

The above code specifies the path to the icon for clearing the DAG.
The `:/` syntax makes use of the resource file.

The Designer file uses the resource URLs `:/draw.svg` and `:/icon_clear.svg` for these icons.
After regenerating the Python module, register the resources with the import described above.
This keeps the icons independent of the working directory used to launch the GUI.

## Documentation

PhyloJunction is documented automatically with `sphinx`, so make sure you have this program installed:

```
python3.11 -m pip install sphinx
python3.11 -m pip install sphinx-new-tab-link
```

`sphinx` is a Python program that converts docstrings into .rst files, and then from those and other .rst files it creates it can produce .html files for online documentation.

The first time documentation was done, a `docs/` directory was manually created in the git repo's root, `PhyloJunction/`.
Then:

```
cd docs/
sphinx-quickstart
```

The answer was "yes" for separate source and build folders and "PhyloJunction" for project name.

Then we needed to activate the Napoleon extension of sphinx, so that we can make use of the Google-style documentation strings PhyloJunction uses.
Inside `source/conf.py` we added:

```
extensions = [ 'sphinx.ext.napoleon' ]
```

With the following options (mostly the defaults):

```
napoleon_google_docstring = True
napoleon_numpy_docstring = False
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False # not using admonitions                                                                                                                          
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = False
napoleon_type_aliases = None
napoleon_attr_annotations = True
```

Because PhyloJunction's code is not in `docs/`, but in fact in `src/phylojunction/`, we must tell sphinx via `conf.py`.
At the top of `conf.py`, add:

```
import os
import sys
sys.path.insert(0, os.path.abspath("../src/phylojunction/"))
```

### Theme and other configurations

We will use the "Read the Docs" (RTD) theme to document PJ, as it is straightforward and familiar.
First, install it:

```
python3.11 -m pip install sphinx-rtd-theme
```

And the following must be added to `conf.py`:

```
html_theme_path = sphinx_bootstrap_theme.get_html_theme_path()
html_theme = "sphinx_rtd_theme"
```

We will also add a few .rst directives to `conf.py`.
These work like "macros" that will be available for all .html files rendered by sphinx, and help us replace strings and add hyperlinks to them automatically.

```
rst_prolog = """
.. |fkm| replace:: Fabio K. Mendes
.. |pj| replace:: PhyloJunction"""

# related to installation, requirements and executing 
rst_prolog += """\
.. |Python| replace:: Python
.. _Python: http://www.python.org
.. |Python3| replace:: Python 3
.. |dendropy| replace:: DendroPy
.. _DendroPy: https://dendropy.org/
.. |msprime| replace:: msprime
.. _msprime: https://tskit.dev/msprime
.. |numpy| replace:: Numpy
.. _NumPy: https://numpy.org/
.. |scipy| replace:: SciPy
.. _SciPy: https://scipy.org
.. |matplotlib| replace:: Matplotlib
.. _Matplotlib: https://matplotlib.org
.. |pandas| replace:: pandas
.. _pandas: https://pandas.pydata.org
.. |pip| replace:: pip
.. _pip: https://pypi.python.org/pypi/pip
.. |pysimplegui| replace:: PySimpleGUI
.. _PySImpleGUI: https://pysimplegui.org/en/latest/"""

# other
rst_prolog += """\
.. _R: https://www.r-project.org
.. |java| replace:: Java
.. _Java: https://www.java.com/en
.. |beast| replace:: BEAST
.. _BEAST: https://beast.community
.. |beast2| replace:: BEAST 2
.. _BEAST 2: https://www.beast2.org
.. |rb| replace:: RevBayes
.. _RevBayes: https://revbayes.github.io
.. |cpp| replace:: C++
.. _C++: https://isocpp.org/"""
```

Now that the configuration files are ready, we can call the program that automatically generate .rst files from all docstrings in PJ:

```
cd docs/
sphinx-apidoc -o ./source ../src/phylojunction
```

Do not forget to run this command again if docstrings are updated.
Finally, we call the command that produces the .html files:

```
make clean && make html
```

This command will need to be called every time we want to update the .html files once the .rst files are modified.

If for some reason sphinx/autodoc starts to complain certain modules cannot be found, it may be that the .rst files sphinx produces automatically need a refresh (module may have been renamed).
The refresh may not be happening because make things they do not need updating.
Just delete all .rst files inside docs/source/, and call the doc generation commands again.

### Documentation website structure

The .rst file that will be converted into `index.html`, the welcome page, is `source/index.rst` (the master .rst file).
This file contains the information for citing and installing PJ.

The remaining .rst files are organized within manually created folders:

* user/index.rst
* developer/index.rst

These files are included in the website via the `.. toctree::` .rst directive in the master .rst file.
