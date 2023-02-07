from os.path import dirname, basename, isfile, join
import glob

# allows you to just do "from chemquest_website.routes import *" in __init__.py
modules = glob.glob(join(dirname(__file__), "*.py"))
__all__ = [ basename(f)[:-3] for f in modules if isfile(f) and not f.endswith('__init__.py')]