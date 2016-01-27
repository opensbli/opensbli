""" This contains utility routines not related to code generationm, such as writing out expanded equations in LaTeX format. """

from sympy import *
from sympy.parsing.sympy_parser import standard_transformations, implicit_application
transformations = standard_transformations + (implicit_application,)
import re
import textwrap

import logging
LOG = logging.getLogger(__name__)


class LatexWriter(object):

    """ Handles writing of equations and arbitrary strings in LaTeX format. """

    def __init__(self):
        return

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
        footer_article = '\\end{document}'
        return

    def write_footer(self):
        """ Write the footer of the LaTeX article file.

        :returns: None
        """

        footer = '\\end{document}'
        self.f.write(footer)
        return

    def write_equations(self, equations, varform=None):
        """ Write out the equations.

        :arg equations: the list of equations.
        :returns: None
        """

        # Transform equations
        out = []

        if isinstance(equations, dict):
            out = out + ['The input is a Dictionary\n\n']
            for key, value in equations.iteritems():
                if isinstance(value, list):
                    out = out + ['The values are a list\n\n']
                    for l in value:
                        out = out + [latex(l, mode='equation', mul_symbol="ldot")]
                    out = out + ['the list of values ends here\n\n']
                else:
                    # out = out + ['The values are not a list\n\n']
                    if (value.is_Equality):
                        temp = latex(value, mode='equation', mul_symbol="ldot")
                    else:
                        temp = '\\begin{dmath}' + latex(key, mul_symbol="ldot") + '=' + latex(value, mul_symbol="ldot") + '\\end{dmath}'
                    out = out + [temp]
                    # out = out + ['The values are not a list and end here\n\n']

        elif isinstance(equations, list):
            out = out + ['The values are a list\n\n']
            for e in equations:
                out = out + [latex(e, mode='equation', long_frac_ratio=3, mul_symbol="ldot")]
            out = out + ['The list ends here\n\n']

        # Replace local tokens with corresponding LaTeX expressions.
        for o in range(len(out)):
            out[o] = out[o].replace('gama', '\gamma').replace('rhou', '\\rho u').replace('rhoo', '\\rho o')
            out[o] = out[o].replace('rhoE', '\\rho E').replace('equation', 'dmath')
            out[o] = out[o].replace('_{OPS_ACC%d', '{').replace('x0', 'i').replace('x1', 'j').replace('x2', 'k')
            out[o] = out[o].replace('{i,j,k}', '{}')
            # out[o] = out[o].replace('_{OPS_ACC%d','{').replace('x0,','').replace('x1,','').replace('x2','')
            out[o] = re.sub(r'\bd\b', '\partial', out[o])
            out[o] = out[o].replace('u_{0}', 'u').replace('u_{1}', 'v').replace('u_{2}', 'w')
            out[o] = out[o].replace('x_{0}', 'x').replace('x_{1}', 'y').replace('x_{2}', 'z')
            out[o] = out[o].replace('\operatorname', '')
        out = out + ['']
        out = ' \n'.join(out)
        out = ' \n'.join(textwrap.wrap(out, width=70, break_long_words=False))
        self.f.write(out)

        return

