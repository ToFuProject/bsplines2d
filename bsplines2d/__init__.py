

from .version import __version__

# #####################
# Debugging unit tests - START
# #####################

# import sys
# import os
# _PATH_HERE = os.path.dirname(__file__)
# _PATH_PROJ = os.path.dirname(os.path.dirname(_PATH_HERE))
# _PATH_DS = os.path.join(_PATH_PROJ, 'datastock')
# sys.path.insert(0, _PATH_DS)
# import datastock as ds
# sys.path.pop(0)

# ######################
# Debugging unit tests - END
# ######################


from ._class03_Bins import Bins as BSplines2D
from . import tests
