
import logging

LOG = logging.getLogger(__name__)
_HANDLER = logging.StreamHandler()
_FORMATTER = logging.Formatter(fmt="%(asctime)s : %(module)s : %(levelname)s : %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
_HANDLER.setFormatter(_FORMATTER)
LOG.addHandler(_HANDLER)
del(_HANDLER)
LOG.setLevel(logging.DEBUG)  # Default to DEBUG level.

from .scheme import *
from .shock_capturing import *
from .averaging import *
from .weno import *
from .teno import *
