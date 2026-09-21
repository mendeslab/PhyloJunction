# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys
sys.path.insert(0, os.path.abspath('../../src')) # PJ's src .py files inside

project = 'PhyloJunction'
copyright = '2023, Fabio K. Mendes'
author = 'Fabio K. Mendes'
release = '0.0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["sphinx.ext.autodoc", "sphinx.ext.napoleon", "sphinx_new_tab_link"]
# "sphinx_new_tab_link" makes new links be opened in separate tabs; it needs module sphinx-new-tab-link

# necessary to load this custom css, which 
# (i) places equation numbers to their right, as opposed to the top
#     of the equations (the default of sphinx)
def setup(app):
    app.add_css_file('css/pj_custom.css')

templates_path = ['_templates']
exclude_patterns = []

#####################
# Napoleon settings #
#####################
napoleon_google_docstring = True
napoleon_numpy_docstring = False
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False # for _membername (e.g., __private_member(self))
napoleon_include_special_with_doc = True # for __membername__ (e.g., __str__(self))
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False # no admonitions in Google style
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = False
napoleon_type_aliases = None
napoleon_attr_annotations = True

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ['_static']

# these two are for placing the logo on the top-left corner
# on the RTD website menu
html_logo = 'images/pj_logo_doc.svg'
html_theme_options = {
    'logo_only': True,
    'display_version': False,
}


rst_prolog = """
.. |fkm| replace:: Fabio K. Mendes
.. |pj| replace:: PhyloJunction"""

# related to installation, requirements and executing 
rst_prolog += """
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
.. _pip: https://pypi.python.org/pypi/pip
.. |pysimplegui| replace:: PySimpleGUI
.. _PySImpleGUI: https://pysimplegui.org/en/latest/
.. |vs| replace:: Visual Studio Code
.. |pych| replace:: PyCharm
.. _PyCharm: https://www.jetbrains.com/pycharm/"""

# other
rst_prolog += """
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

