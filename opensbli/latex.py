#!/usr/bin/env python

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

from sympy import *
from sympy.printing.latex import *
from sympy.parsing.sympy_parser import *
init_printing(use_latex=True)
import re
import textwrap

import logging
LOG = logging.getLogger(__name__)


class LatexWriter(LatexPrinter):

    """ Handles writing of equations and arbitrary strings in LaTeX format. """

    def open(self, path):
        """ Open the LaTeX file for writing.

        :arg str path: the path to the LaTeX file.
        :returns: None
        """
        self.f = open(path, "w")
        return

    def close(self):
        """ Close the LaTeX file.

        :returns None
        """
        self.f.close()
        return

    def write_string(self, s):
        """ Write an arbitrary user-provided string.

        :arg str s: a user-provided string.
        :returns: None
        """
        self.f.write(s)
        return

    def write_header(self, metadata):
        """ Write the header of the LaTeX article file.

        :arg dict metadata: a dictionary containing the title, author and institution.
        :returns: None
        """

        header = []
        header.append('\\documentclass{article}')
        header.append('\\title{%s}\n\\author{%s\\\\ %s}' % (metadata["title"], metadata["author"], metadata["institution"]))
        header.append('\\date{\\today}')
        header.append('\\usepackage{color}\n\\usepackage{breqn}')
        header.append('\\begin{document}')
        header.append('\\maketitle')
        header.append('')
        header = '\n'.join(header)

        self.f.write(header)
        return

    def write_footer(self):
        """ Write the footer of the LaTeX article file.

        :returns: None
        """

        footer = '\\end{document}'
        self.f.write(footer)
        return

    def latexify_expression(self, expression):
        """ Format a SymPy expression as LaTeX.

        :arg expression: The SymPy expression to format.
        :returns: The LaTeX string representation of the expression.
        :rtype: str
        """

        self._settings['mode'] = "dmath"
        self._settings['long_frac_ratio'] = 3
        self._settings['mul_symbol'] = "cdot"

        tex = Printer.doprint(self, expression)

        if self._settings['mode'] == 'plain':
            output = tex
        elif self._settings['mode'] == 'inline':
            output = r"$%s$" % tex
        elif self._settings['itex']:
            output = r"$$%s$$" % tex
        else:
            env_str = self._settings['mode']
            output = r"\begin{%s}%s\end{%s}" % (env_str, tex, env_str)
        return output

    def _print_Indexed(self, expr):
        tex = '{%s}' % self._print(expr.base)+'_{%s}' % ','.join( map(self._print, expr.indices))
        return tex

    def write_expression(self, expression, substitutions = {}):
        """ Convert a single expression or list of expressions to LaTeX format and then write it/them to file.

        :arg expression: a single SymPy expression of list of SymPy expressions to format and write.
        :returns: None
        """

        output = []

        if isinstance(expression, list):
            for e in expression:
                output.append(self.latexify_expression(e))
        else:
            output.append(self.latexify_expression(expression))

        # Perform any user-defined LaTeX substitutions
        for key, value in substitutions.iteritems():
            for i in range(len(output)):
                output[i] = output[i].replace(key, value)

        output = ' \n'.join(output)
        output = ' \n'.join(textwrap.wrap(output, width=70, break_long_words=False))

        self.f.write(output)
        return

def write_descritization():
    
    return