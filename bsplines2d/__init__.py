# ###############
# __version__
# ###############


from . import _version
__version__ = _version.version
__version_tuple__ = _version.version_tuple


# ######################
# sub-packages
# ######################


from ._class03_Bins import Bins as BSplines2D
from ._saveload import load
from . import tests
