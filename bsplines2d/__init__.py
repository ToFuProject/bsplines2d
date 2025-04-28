# ###############
# __version__
# ###############


from . import _version
__version__ = _version.version
__version_tuple__ = _version.version_tuple


# ######################
# sub-packages
# ######################


from ._class04_Bins import Bins as Collection
from ._saveload import load
from . import tests
