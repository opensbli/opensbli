''' This contains the routines used for algorithm generation
'''
from sympy import *
from sympy.parsing.sympy_parser import (parse_expr, standard_transformations,implicit_application)
transformations = standard_transformations + (implicit_application,)
import re
def latex_article_header(inp):
  header_article = []
  header_article.append('\\documentclass{article}')
  header_article.append('\\title{%s}\n\\author{%s\\\\ %s}'%(inp[0],inp[1],inp[2]))
  #header_article.append('\\date{\\today}')
  header_article.append('\\usepackage{color}\n\\usepackage{breqn}')
  header_article.append('\\begin{document}')
  header_article.append('\\maketitle')
  header_article.append('')
  header_article = '\n'.join(header_article)

  end_article = '\\end{document}'
  return header_article,end_article

def write_latex(fname,inp,varform=None):
  out = []

  if isinstance(inp, dict):
    out = out + ['The input is a Dictionary\n\n']
    for key, value in inp.iteritems():
      if isinstance(value, list):
        out = out + ['The values are a list\n\n']
        for l in value:
          out = out +[latex(l, mode='equation',mul_symbol="ldot")]
        out = out + ['the list of values ends here\n\n']
      else:
        #out = out + ['The values are not a list\n\n']
        if (value.is_Equality):
          temp =  latex(value, mode='equation',mul_symbol="ldot")
        else:
          temp = '\\begin{dmath}' + latex(key,mul_symbol="ldot") + '=' + latex(value,mul_symbol="ldot") +'\\end{dmath}'
        out = out + [temp]
        #out = out + ['The values are not a list and end here\n\n']

  elif isinstance(inp, list):
    out = out + ['The values are a list\n\n']
    for l in inp:
      out = out +[latex(l, mode='equation',  long_frac_ratio=2,mul_symbol="ldot")]
    out = out + ['The list ends here\n\n']

  for o in range(len(out)):
    out[o] = out[o].replace('gama','\gamma').replace('rhou','\\rho u').replace('rhoo','\\rho o')
    out[o] = out[o].replace('rhoE','\\rho E').replace('equation','dmath')
    out[o] = out[o].replace('_{OPS_ACC%d','{').replace('x0','i').replace('x1','j').replace('x2','k')
    out[o] = out[o].replace('{i,j,k}','{}')
    #out[o] = out[o].replace('_{OPS_ACC%d','{').replace('x0,','').replace('x1,','').replace('x2','')
    out[o] = re.sub(r'\bd\b', '\partial', out[o])
    out[o] = out[o].replace('u_{0}','u').replace('u_{1}','v').replace('u_{2}','w')
    out[o] = out[o].replace('x_{0}','x').replace('x_{1}','y').replace('x_{2}','z')
    out[o] = out[o].replace('\operatorname','')
  out = out + ['']
  write = '\n'.join(out)
  fname.write(write)
  return