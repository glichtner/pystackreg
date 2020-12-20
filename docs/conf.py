import os
import time

# cannot directly import because __init__.py imports pystackreg which imports the
# compiled plugin, which is not available before setup.py is run
__version__ = ""  # placeholder for linters
exec(open("../pystackreg/version.py").read())

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.coverage',
    'sphinx.ext.doctest',
    'sphinx.ext.extlinks',
    'sphinx.ext.ifconfig',
    'sphinx.ext.napoleon',
    'sphinx_autodoc_typehints',
    'sphinx.ext.todo',
    'sphinx.ext.viewcode',
]


if os.getenv('SPELLCHECK'):
    extensions += 'sphinxcontrib.spelling',
    spelling_show_suggestions = True
    spelling_lang = 'en_US'

source_suffix = '.rst'
master_doc = 'index'
project = 'pystackreg'
author = "Gregor Lichtner, Philippe Thevenaz"
year = f'2018-{time.strftime("%Y")}'
copyright = '{0}, {1}'.format(year, author)
version = release = __version__

pygments_style = 'trac'
templates_path = ['.']
extlinks = {
    'issue': ('https://github.com/glichtner/pystackreg/%s', '#'),
    'pr': ('https://github.com/glichtner/pystackreg/pull/%s', 'PR #'),
}
# on_rtd is whether we are on readthedocs.org
on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

if not on_rtd:  # only set the theme if we're building docs locally
    html_theme = 'sphinx_rtd_theme'

html_use_smartypants = True
html_last_updated_fmt = '%b %d, %Y'
html_split_index = False
html_sidebars = {
   '**': ['searchbox.html', 'globaltoc.html', 'sourcelink.html'],
}
html_short_title = '%s-%s' % (project, version)

napoleon_use_ivar = True
napoleon_use_rtype = True
napoleon_use_param = True
