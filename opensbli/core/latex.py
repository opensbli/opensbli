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
        self._settings['mul_symbol_latex'] = r" \,.\, "
        #self._settings['mul_symbol_latex'] = \
            #mul_symbol_table[self._settings['mul_symbol']]

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
        ind = list(set(expr.indices))
        if len(ind) ==1 and ind[0] == 0:
            tex = '{%s}' % self._print(expr.base)
        else:
            tex = '{%s}' % self._print(expr.base)+'[{%s}]' % ','.join(map(self._print, expr.indices))
        return tex
    def _print_Pow(self, expr):
        base, exponent = expr.as_base_exp()
        tex = '{%s}^{%s}'%(self._print(base),self._print(exponent))
        return tex
    def _print_EinsteinTerm(self, expr):
        return str(expr)

    def _print_CentralDerivative(self, expr):
        return r'\left. %s \right|_{{%s }}' % (self._print(Derivative(*expr.args)), "Central")
    def _print_MetricDerivative(self, expr):
        return r'\left. %s \right|_{{%s }}' % (self._print(Derivative(*expr.args)), "Metric")
    def _print_TemporalDerivative(self, expr):
        return r'\left. %s \right|_{{%s }}' % (self._print(Derivative(*expr.args)), "Temporal")
    def _print_WenoDerivative(self, expr):
        return r'\left. %s \right|_{{%s }}' % (self._print(Derivative(*expr.args)), "Weno")

    def write_expression(self, expression, substitutions={}):
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
        output =  output + '\n\n'
        self.f.write(output)
        return
#from .kernel import Kernel
class WriteAlgorithm():
    def __init__(self, code):
        self.substitutions = {'gama':'\gamma ', 'deltat':'\delta t ','rhoE':'\\rho E','deltai0':'\Delta_{i0}',\
            'rhou0':'\\rho u0', 'rhou1':'\\rho u1', 'rhou2':'\\rho u2', 'deltai1':'\Delta_{i1}','deltai2':'\Delta_{i2}'}
        self.latex = LatexWriter()
        self.latex.open(code.CODE_DIR + "/algorithm_%s_%dd.tex"%(code.simulation_parameters["name"], code.ndim))
        self.write_header(code)
        self.write_computation_kernels(code)
        self.latex.write_footer()
        self.latex.close()
        return
    def write_header(self, code):
        """Writes the algorithm from the code object, this code object would be the one returned by
        OPSC or other programming languages. File name is automatically generated from the simulation
        name given in simulation parameters"""
        metadata = {"title": "Algorithm", "author": "Satya P Jammy", "institution": "University of Southampton"}
        self.latex.write_header(metadata)

        return
    def write_footer(self, code):

        return
    def write_computation_kernels(self, code):

        equations_code = []
        for block in range(code.nblocks):
            # Get all the computations to be performed. Add computations as needed.
            block_computations = []
            self.latex.write_string('The following are the computations performed on block %s\n\n'%str(block))
            if code.spatial_discretisation[block].computations:
                block_computations += code.spatial_discretisation[block].computations
            if code.temporal_discretisation[block].computations:
                block_computations += code.temporal_discretisation[block].computations
            if code.temporal_discretisation[block].start_computations:
                block_computations += code.temporal_discretisation[block].start_computations
            if code.temporal_discretisation[block].end_computations:
                block_computations += code.temporal_discretisation[block].end_computations
            if code.initial_conditions[block].computations:
                block_computations += code.initial_conditions[block].computations
            if code.diagnostics:
                for inst in code.diagnostics[block]:
                    block_computations += inst.computations
            if code.boundary_condition[block].computations:
                block_computations += [t for t in code.boundary_condition[block].computations if isinstance(t, Kernel)]

            for computation in block_computations:
                self.write_kernel(computation)
        return
    def write_kernel(self, computation):
        print(computation.computation_type)
        self.latex.write_string("The computation performed is %s \n\n"%computation.computation_type)
        self.latex.write_string("The name of the kernel given in the code is %s\n\n"%computation.name)
        self.latex.write_string("The grid range on which the kernel is executed is %s\n\n"%computation.ranges)
        if isinstance(computation.equations, list):
            for eq in computation.equations:
                self.latex.write_expression(eq, self.substitutions)
        else:
            self.latex.write_expression(computation.equations, self.substitutions)
        self.latex.write_string("\n \n")
        return


#def create_latex_kernel(kernels):
    #latex = LatexWriter()
    #latex.open('./kernels.tex')
    #metadata = {"title": "Characteristic boundary conditions", "author": "Jammy", "institution": ""}
    #latex.write_header(metadata)
    #for ker in kernels:
        ##pprint(ker.__dict__)
        #latex.write_string('The kernel is %s'%ker.computation_name)
        #for index, eq in enumerate(ker.equations):
            #if isinstance(eq, Matrix):
                #for n, e in enumerate(eq):
                    #v = Eq(Symbol('dW0%d' %n), e)
                    #latex.write_expression(v)
            #else:
                #latex.write_expression(eq)
    #latex.write_footer()
    #latex.close()
    #return