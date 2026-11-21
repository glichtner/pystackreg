import os
import time
from pathlib import Path

# The version is read from the source file so that the configuration does not
# depend on which pystackreg is importable; autodoc itself needs an installed
# (compiled) pystackreg.
__version__ = ""
exec((Path(__file__).parent.parent / "pystackreg" / "version.py").read_text())

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.coverage",
    "sphinx.ext.doctest",
    "sphinx.ext.extlinks",
    "sphinx.ext.ifconfig",
    "sphinx.ext.napoleon",
    "sphinx.ext.todo",
    "sphinx.ext.viewcode",
    "sphinx.ext.mathjax",
]

exclude_patterns = ["_build"]
autodoc_typehints = "none"


if os.getenv("SPELLCHECK"):
    extensions += ("sphinxcontrib.spelling",)
    spelling_show_suggestions = True
    spelling_lang = "en_US"

source_suffix = {".rst": "restructuredtext"}
master_doc = "index"
project = "pystackreg"
author = "Gregor Lichtner"
year = f'2018-{time.strftime("%Y")}'
copyright = "{0}, {1}".format(year, author)
version = release = __version__

pygments_style = "trac"
templates_path = ["."]
extlinks = {
    "issue": ("https://github.com/glichtner/pystackreg/issues/%s", "#%s"),
    "pr": ("https://github.com/glichtner/pystackreg/pull/%s", "PR #%s"),
}
html_theme = "sphinx_rtd_theme"
html_last_updated_fmt = "%b %d, %Y"
html_split_index = False
html_sidebars = {
    "**": ["searchbox.html", "globaltoc.html", "sourcelink.html"],
}
html_short_title = "%s-%s" % (project, version)

napoleon_use_ivar = True
napoleon_use_rtype = True
napoleon_use_param = True
