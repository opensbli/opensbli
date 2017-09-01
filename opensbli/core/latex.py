from sympy import Derivative
from sympy.printing.latex import LatexPrinter, Printer
import textwrap
import logging
import os
LOG = logging.getLogger(__name__)


class LatexWriter(LatexPrinter):

    """ Handles writing of equations and arbitrary strings in LaTeX format. """

    def open(self, path, title):
        """ Open the LaTeX file for writing.

        :arg str path: the path to the LaTeX file.
        :returns: None
        """
        dir = "./" + "latex_output"
        if not os.path.exists(dir):
            os.makedirs(dir)
        file_path = dir + '/' + path
        if os.path.exists(file_path):
            print ("The latex file %s already exists overwriting it" % file_path)
        self.f = open(file_path, "w")
        self.title = title
        self.write_header
        return

    def close(self):
        """ Close the LaTeX file.

        :returns: None
        """
        self.write_footer
        self.f.close()
        return

    def write_string(self, s):
        """ Write an arbitrary user-provided string.

        :arg str s: a user-provided string.
        :returns: None
        """
        self.f.write("\\noindent %s\\\\" % s)
        return

    @property
    def write_header(self):
        """ Write the header of the LaTeX article file.

        :arg dict metadata: a dictionary containing the title, author and institution.
        :returns: None
        """

        header = []
        header.append('\\documentclass{article}')
        header.append('\\title{%s}\n\\author{Satya P Jammy \& OpenSBLI developers team}' % (self.title))
        header.append('\\date{\\today}')
        header.append('\\usepackage{color}\n\\usepackage{amsmath}\n\\usepackage{breqn}')
        header.append('\\begin{document}')
        header.append('\\maketitle')
        header.append('')
        header = '\n'.join(header)
        self.f.write(header)
        return

    @property
    def write_footer(self):
        """ Write the footer of the LaTeX article file.

        :returns: None
        """

        footer = '\\end{document}'
        self.f.write(footer)
        return

    def latexify_expression(self, expression, mode):
        """ Format a SymPy/OpenSBLI expression as LaTeX.

        :arg expression: The expression to convert.
        :returns: The LaTeX string representation of the expression.
        :rtype: str
        """

        self._settings['mode'] = mode
        self._settings['long_frac_ratio'] = 3
        self._settings['mul_symbol_latex'] = r" \,.\, "
        # self._settings['mul_symbol_latex'] = \
        # mul_symbol_table[self._settings['mul_symbol']]

        tex = Printer.doprint(self, expression)

        if self._settings['mode'] == 'plain':
            output = r"\noindent$%s$\\" % tex
        elif self._settings['mode'] == 'inline':
            output = r"$%s$" % tex
        elif self._settings['itex']:
            output = r"$$%s$$" % tex
        else:
            env_str = self._settings['mode']
            output = r"\begin{%s}%s\end{%s}" % (env_str, tex, env_str)
        return output

    def _print_DataSetBase(self, expr):
        tex = "%s{_{B%s}}" % (self._print(expr.label), self._print(expr.blocknumber))
        return tex

    def _print_DataSet(self, expr):
        ind = list(expr.indices)[:-1]
        tex = '{%s}' % self._print(expr.base)+'[{%s}]' % ','.join(map(self._print, ind))
        return tex

    def _print_Indexed(self, expr):
        tex = '{%s}' % self._print(expr.base)+'[{%s}]' % ','.join(map(self._print, expr.indices))
        return tex

    def _print_Pow(self, expr):
        base, exponent = expr.as_base_exp()
        tex = '\left(%s \\right)^{%s}' % (self._print(base), self._print(exponent))
        return tex

    def _print_EinsteinTerm(self, expr):
        return str(expr)

    def _print_CentralDerivative(self, expr):
        # return r'\left. %s \right|_{{%s }}' % (self._print(Derivative(*expr.args)), "Central")
        return r'%s' % (self._print(Derivative(*expr.args)))

    def _print_MetricDerivative(self, expr):
        return r'\left. %s \right|_{{%s }}' % (self._print(Derivative(*expr.args)), "Metric")

    def _print_TemporalDerivative(self, expr):
        return r'\left. %s \right|_{{%s }}' % (self._print(Derivative(*expr.args)), "Temporal")

    def _print_WenoDerivative(self, expr):
        return r'\left. %s \right|_{{%s }}' % (self._print(Derivative(*expr.args)), "Weno")

    def _print_TenoDerivative(self, expr):
        return r'\left. %s \right|_{{%s }}' % (self._print(Derivative(*expr.args)), "Teno")

    def _print_KD(self, expr):
        # print expr
        # print expr._latex_no_arg(self)
        # for arg in expr.args:
            # print self._print(arg)
        # exit()
        return expr._latex_no_arg(self) + '_{%s}' % (' '.join([self._print(a) for a in expr.args]).replace('_', ''))

    def write_expression(self, expression, mode=None, substitutions={}):
        """ Convert a single expression or list of expressions to LaTeX format and then write it/them to file.

        :arg expression: a single SymPy expression of list of SymPy expressions to format and write.
        :arg mode: weather the latex is an inline expression or equation type or any that user wants. defaults to `dmath` mode of LaTeX
        :arg substitutions: additional substitutions the user wants in the output
        :returns: None
        """

        output = []
        if not mode:
            mode = "dmath"

        if isinstance(expression, list):
            for e in expression:
                output.append(self.latexify_expression(e, mode))
        else:
            output.append(self.latexify_expression(expression, mode))

        # Perform any user-defined LaTeX substitutions
        for key, value in substitutions.iteritems():
            for i in range(len(output)):
                output[i] = output[i].replace(key, value)

        output = ' \n'.join(output)
        output = ' \n'.join(textwrap.wrap(output, width=70, break_long_words=False))
        output = output + '\n\n'
        self.f.write(output)
        return
