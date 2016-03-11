#!/usr/bin/env python

#    AutoFD: An automatic code generator for solving differential equations.
#    Copyright (C) 2016 Satya P. Jammy, Christian T. Jacobs

#    This file is part of AutoFD.

#    AutoFD is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    AutoFD is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with AutoFD.  If not, see <http://www.gnu.org/licenses/>.

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
        # header.append('\\date{\\today}')
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
        
        self._settings['mode'] = "equation"
        self._settings['long_frac_ratio'] = 3
        self._settings['mul_symbol'] = "cdot"
    
        tex = Printer.doprint(self, expression)
        
#        # Replace local tokens with corresponding LaTeX expressions.
#        for o in range(len(out)):
#            out[o] = out[o].replace('gama', '\gamma').replace('rhou', '\\rho u').replace('rhoo', '\\rho o')
#            out[o] = out[o].replace('rhoE', '\\rho E').replace('equation', 'dmath')
#            out[o] = out[o].replace('_{OPS_ACC%d', '{').replace('x0', 'i').replace('x1', 'j').replace('x2', 'k')
#            out[o] = out[o].replace('{i,j,k}', '{}')
#            # out[o] = out[o].replace('_{OPS_ACC%d','{').replace('x0,','').replace('x1,','').replace('x2','')
#            out[o] = re.sub(r'\bd\b', '\partial', out[o])
#            out[o] = out[o].replace('u_{0}', 'u').replace('u_{1}', 'v').replace('u_{2}', 'w')
#            out[o] = out[o].replace('x_{0}', 'x').replace('x_{1}', 'y').replace('x_{2}', 'z')
#            out[o] = out[o].replace('\operatorname', '')

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
        
    def write_expression(self, expression):
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

        output = ' \n'.join(output)
        output = ' \n'.join(textwrap.wrap(output, width=70, break_long_words=False))

        self.f.write(output)
        return
        
#    def write_equations(self, equations, varform=None):
#        """ Write out the equations.

#        :arg equations: the list of equations.
#        :returns: None
#        """

#        # Transform equations
#        out = []

#        if isinstance(equations, dict):
#            out = out + ['The input is a Dictionary\n\n']
#            for key, value in equations.iteritems():
#                if isinstance(value, list):
#                    out = out + ['The values are a list\n\n']
#                    for l in value:
#                        out = out + [latex(l, mode='equation', mul_symbol="ldot")]
#                    out = out + ['the list of values ends here\n\n']
#                else:
#                    # out = out + ['The values are not a list\n\n']
#                    if (value.is_Equality):
#                        temp = latex(value, mode='equation', mul_symbol="ldot")
#                    else:
#                        temp = '\\begin{dmath}' + latex(key, mul_symbol="ldot") + '=' + latex(value, mul_symbol="ldot") + '\\end{dmath}'
#                    out = out + [temp]
#                    # out = out + ['The values are not a list and end here\n\n']

#        elif isinstance(equations, list):
#            out = out + ['The values are a list\n\n']
#            for e in equations:
#                out = out + [latex(e, mode='equation', long_frac_ratio=3, mul_symbol="ldot")]
#            out = out + ['The list ends here\n\n']

#        # Replace local tokens with corresponding LaTeX expressions.
#        for o in range(len(out)):
#            out[o] = out[o].replace('gama', '\gamma').replace('rhou', '\\rho u').replace('rhoo', '\\rho o')
#            out[o] = out[o].replace('rhoE', '\\rho E').replace('equation', 'dmath')
#            out[o] = out[o].replace('_{OPS_ACC%d', '{').replace('x0', 'i').replace('x1', 'j').replace('x2', 'k')
#            out[o] = out[o].replace('{i,j,k}', '{}')
#            # out[o] = out[o].replace('_{OPS_ACC%d','{').replace('x0,','').replace('x1,','').replace('x2','')
#            out[o] = re.sub(r'\bd\b', '\partial', out[o])
#            out[o] = out[o].replace('u_{0}', 'u').replace('u_{1}', 'v').replace('u_{2}', 'w')
#            out[o] = out[o].replace('x_{0}', 'x').replace('x_{1}', 'y').replace('x_{2}', 'z')
#            out[o] = out[o].replace('\operatorname', '')
#        out = out + ['']
#        out = ' \n'.join(out)
#        out = ' \n'.join(textwrap.wrap(out, width=70, break_long_words=False))
#        self.f.write(out)

#        return
        
l = LatexWriter()
l.open("test.tex")
l.write_expression(["test"])
