# Configuration of the pages of the website of the prototype (Sphinx, MyST
# Markdown). The documentation of TFEL, MGIS and mfem-mgis is built
# separately and copied next to these pages by scripts/build-site.sh.

project = "MMM-GPU-Proto"
copyright = "2026"
language = "en"
extensions = ["myst_parser"]
source_suffix = {".md": "markdown"}
exclude_patterns = ["_build", "requirements.txt"]
myst_heading_anchors = 3
html_theme = "sphinx_rtd_theme"
html_title = "MMM-GPU-Proto"
html_show_sourcelink = False
html_copy_source = False
html_static_path = ["_static"]
html_css_files = ["custom.css"]
