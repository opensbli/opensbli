#    OpenSBLI: An automatic code generator for solving differential equations.
#    Copyright (C) 2016 Satya P. Jammy, Christian T. Jacobs

#    This file is part of OpenSBLI.

#    OpenSBLI is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    OpenSBLI is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with OpenSBLI.  If not, see <http://www.gnu.org/licenses/>.

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
from .hybrid import *
