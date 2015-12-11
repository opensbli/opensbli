''' This contains the routines used for codegeneration
'''
from sympy import *
from sympy.parsing.sympy_parser import (parse_expr, standard_transformations,implicit_application)
transformations = standard_transformations + (implicit_application,)
import re
import textwrap
#from algorithm import eqto_dict as eqtd
line_comment = {}
line_comment['OPSC'] = '//'
line_comment['F90'] = '!'
lend = {}
lend['OPSC'] = ';'
lend['F90'] = ''
def eqto_dict(inp):
  lh = list(eq.lhs for eq in inp)
  rh = list(eq.rhs for eq in inp)
  dict_list = zip(lh,rh)
  diction = dict(dict_list)
  #print(diction)
  return diction
def OPSC_write_kernel(eqs, inp):
  def get_kernel(evals):
    lh = flatten(list(list(eq.lhs.atoms(Indexed)) for eq in evals))
    rh = flatten(list(list(eq.rhs.atoms(Indexed)) for eq in evals))
    tot_indexed = list(set(lh+rh))
    libs = set(list(i.base.label for i in lh))
    ribs = set(list(i.base.label for i in rh))
    inouts = libs.intersection(ribs); ins = ribs.difference(inouts)
    outs = libs.difference(inouts)
    inouts= list(inouts); ins =list(ins); outs = list(outs)
    tot_base = ins+outs+inouts
    symbs = flatten(list(list(eq.atoms(Symbol)) for eq in evals))
    Idxs = list(set(list(i.indices for i in (lh+rh))))
    Idxs = flatten(list(i) for i in Idxs)
    labes = list(set(list(i.label for i in Idxs)))
    for i in Idxs:
      labes = labes + [i.lower, i.upper]
    symbs = list(set(symbs).difference(set(labes)).difference(set(inp.const)).difference(set(tot_base)))
    symbs = list(set(symbs))
    out = []
    symdec = []
    evdict = eqto_dict(evals)
    for sym in symbs:
      symdec = symdec + ['double %s;'%sym]
      if evdict.get(sym):
        pass
      else:
        raise ValueError("I dont know the formula for %s"%sym)
    # grid range
    lower = []
    upper = []
    for dim in range(inp.ndim):
      lower = lower + [inp.blockdims[dim].lower - inp.halos[dim]]
      upper = upper + [inp.blockdims[dim].upper + inp.halos[dim]+1]
    #print(lower,upper)
    for ev in evals:
      code = ccode(ev)
      code = code.replace('==','=') + lend['OPSC']
      out =out + [code]
    kercall = []; kerheader = []; kernel = []
    kername = inp.kername%inp.kernel_ind
    inp.kernel_ind = inp.kernel_ind + 1
    kerca = 'ops_par_loop(%s, \"%s\", %s[%s], %d, %%s'%(kername,kername,inp.blkname, inp.block, inp.ndim)
    head = 'void %s('%kername
    if tot_base:
      for ind,v in enumerate(tot_base):
        if v in ins:
          opstype = 'OPS_READ'
          headty = 'const double *%s'
        elif v in outs:
          opstype = 'OPS_WRITE'
          headty = 'double *%s'
        elif v in inouts:
          opstype = 'OPS_RW'
          headty = 'double *%s'
        else:
          raise ValueError('Dont know what the base is %s'%v)
        varib = flatten(list(v1 for v1 in tot_indexed if v1.base.label == v))
        varib = list(set(varib))
        variabind = flatten(list(v1.indices) for v1 in varib)
        variabind = list(set(variabind))

        if all(va.upper == inp.block.upper for va in variabind ):
          indexes = list(va for va in variabind)
          for dim in range(inp.ndim):
            indexes = list(str(te).replace('x%d'%dim,'0') for te in indexes)
          indexes = list(parse_expr(v) for v in indexes)
          for inde in range(len(indexes)):
            for ou in range(len(out)):
              if isinstance(indexes[inde], tuple):
                new = '%s[OPS_ACC%d%s]'%(v,ind,indexes[inde])
              else:
                new = '%s[OPS_ACC%d(%s)]'%(v,ind,indexes[inde])
              old = ('%s\[%s\]'%(v,variabind[inde])).replace('+','\+')
              out[ou] = re.sub(r"\b(%s)"%old, new, out[ou])
              #out[ou] = out[ou].replace(old, new)
          # get the stencil name to be written
          indexes = indexes + [parse_expr(', '.join(list(str(0) for dim in range(inp.ndim))))]
          indexes =  list(set(indexes))
          if inp.ndim > 1 :
            for dim in range(inp.ndim):
              indexes = sorted(indexes, key=lambda indexes: indexes[dim])
            temp = flatten(list(list(t) for t in indexes))
          else:
            indexes = [sorted(indexes)]
            temp = flatten(list(t) for t in indexes)

          sten = ','.join(list(str(t) for t in temp))
          if inp.stencils.get(sten):
            sten_na = inp.stencils.get(sten)
          else:
            sten_na = inp.sten_name%inp.sten_ind
            inp.stencils[sten] = sten_na
            inp.sten_ind = inp.sten_ind +1

          # update range on which the loop to be iterated
          if len(indexes) == 1:
            for dim in range(inp.ndim):
              lower[dim] = lower[dim]- indexes[0][dim]
              upper[dim] = upper[dim]- indexes[0][dim]
          else:
            for dim in range(inp.ndim):
              lower[dim] = lower[dim]- indexes[0][dim]
              upper[dim] = upper[dim]- indexes[-1][dim]
          datatype = 'double'
          arg_call = '%%s(%%s[%s], 1, %%s, \"%%s\", %%s)'%inp.block
          call = arg_call%('ops_arg_dat',v,sten_na,datatype,opstype)
          kercall=  kercall + [call]
          kerheader = kerheader + [headty%v]
        else:
          indexes = list(va for va in variabind)
          indexes = list(str(te).replace(str(te),'0') for te in indexes)
          indexes = list(parse_expr(v) for v in indexes)

          for inde in range(len(indexes)):
            for ou in range(len(out)):
              temp = [indexes[inde]]
              temp = list(str(te) for te in temp)
              new = '%s[%s]'%(v,','.join(temp))
              old = str(varib[inde])
              out[ou] = out[ou].replace(old, new)
          datatype = 'double'
          arg_call = '%%s(&%%s[%s], 1, \"%%s\", %%s)'%variabind[0]
          call = arg_call%('ops_arg_gbl',v,datatype,opstype)
          kercall=  kercall + [call]
          kerheader = kerheader + [headty%v]
      iter_range = []
      for dim in range(inp.ndim):
        iter_range = iter_range + [str(lower[dim])] + [str(upper[dim])]
      iter_range = ','.join(iter_range)
      kercall.insert(0, kerca%'iter_range%d'%inp.iterrange);
      for indno in range(len(kercall)-1):
        kercall[indno] = kercall[indno] + ','
      kercall[-1] = kercall[-1] + ');'
      #kercall = ',\n'.join(kercall)
      kercall =  ['int iter_range%d[] = {%s};\n'%(inp.iterrange,iter_range)] +kercall
      inp.iterrange = inp.iterrange + 1
      kerheader = head + ', '.join(kerheader) + '){'
      kernel = [kerheader] + symdec + out + ['}']; #kernel = '\n'.join(kernel)
    else:
      print(tot_base)
      pass
    return kercall, kernel
  allcalls = [];  allkernels = []
  if isinstance(eqs, dict):
    for key, value in eqs.iteritems():
      if isinstance(value, list):
        call, comp = get_kernel(value)
        allcalls = allcalls + [call]; allkernels = allkernels + [comp]
      else:
        call, comp = get_kernel([value])
        allcalls = allcalls + [call]; allkernels = allkernels + [comp]
  elif isinstance(eqs, list):
    call, comp = get_kernel(eqs)
    allcalls = allcalls + [call]; allkernels = allkernels + [comp]
  else:
    call, comp = get_kernel([eqs])
    allcalls = allcalls + [call]; allkernels = allkernels + [comp]
  #pprint('\n kernel is')
  #print('\n\n'.join(allkernels))
  #pprint('\n Call is')
  #print('\n\n'.join(allcalls))
  return allcalls, allkernels


def header_code(inp,alg):
  out = []
  lang = alg.lang
  if alg.lang == 'OPSC':
    out = out + ['#include <stdlib.h>']
    out = out + ['#include <string.h>']
    out = out + ['#include <math.h>']

    out.append('%s Global Constants in the equations are'%line_comment[lang])
    dobs = []
    ints = []
    for con in inp.const:
      if isinstance(con, Symbol):
        if con.is_integer:
          ints = ints + ['%s'%con]
        else:
          dobs = dobs + ['%s'%con]
      elif isinstance(con, Indexed):
        tot = 0
        for ind in con.shape:
          tot = tot + ind
        if con.is_integer:
          ints = ints + ['%s[%d]'%(con.base,tot)]
        else:
          dobs = dobs + ['%s[%d]'%(con.base,tot)]
    if ints:
      out = out + ['int %s %s'%(', '.join(ints),lend[lang])]
    if dobs:
      out = out + ['double %s %s'%(', '.join(dobs),lend[lang])]

    out= out +['// OPS header file']
    out = out + ['#define OPS_%sD'%inp.ndim]
    out = out + ['#include "ops_seq.h"']
    out = out + ['#include "auto_kernel.h"']
    out = out + ['%s main program start'%line_comment[lang]]
    out = out + ['int main (int argc, char **argv) {']
  elif lang == 'F90':
    dobs = []
    ints = []
    for con in inp.const:
      if isinstance(con, Symbol):
        if con.is_integer:
          ints = ints + ['%s'%con]
        else:
          dobs = dobs + ['%s'%con]
      elif isinstance(con, Indexed):
        tot = 0
        print(fcode(con))
        for ind in con.shape:
          tot = tot + ind
        if con.is_integer:
          ints = ints + ['%s(%d)'%(con.base,tot)]
        else:
          dobs = dobs + ['%s(%d)'%(con.base,tot)]
    if ints:
      out = out + ['integer :: %s %s'%(', '.join(ints),lend[lang])]
    if dobs:
      out = out + ['real(8) :: %s %s'%(', '.join(dobs),lend[lang])]
    inp.module.append('\n'.join(out))
    out = []
    #TODO spj change module name later as of now using the same module name parammod
    out = out + ['program codegen']
    out = out + ['use param_mod']
    out = out + ['IMPLICIT NONE']

  else:
    raise ValueError('Implement %s in the header declaration'%lang)

  return out
def loop(indices,alg):
  ''' this reuqires an index as input'''
  lstart = []
  lend = []
  if alg.lang == 'OPSC':
    comm = '%s loop start %s'%(line_comment.get(alg.lang),indices[0])
    lstart.append(comm)
    for dim in range(0,len(indices)):
      temp = indices[dim]
      lstart.append('for(int %s=%s; %s<%s; %s++){'%(temp,temp.lower,temp,temp.upper,temp))
      lend.append('}')
    comm = '%s loop end %s'%(line_comment.get(alg.lang),indices[0])
    lend.append(comm)
  elif alg.lang == 'F90':
    comm = '%s loop start %s'%(line_comment.get(alg.lang),indices[0])
    lstart.append(comm)
    for dim in reversed(range(0,len(indices))):
      temp = indices[dim]
      lstart.append('do %s=%s,%s'%(temp,temp.lower,temp.upper))
      lend.append('enddo')
    comm = '%s loop end %s'%(line_comment.get(alg.lang),indices[0])
    lend.append(comm)

  else:
    raise ValueError('Implement %s in the loop declaration'%alg.lang)
  return lstart,lend
def defdec_lang(inp,alg):
  out = []
  lang = alg.lang
  if alg.lang == 'OPSC':
    # Define inputs to the code
    blkname = inp.blkname
    totblock = inp.block.upper - inp.block.lower
    ind = inp.blockdims
    inputs = []
    inputs = inputs + ['int %s = %d %s'%(totblock,inp.nblocks,lend[lang])]
    ## change this
    inputs = inputs + ['int %s %s'%(','.join(list('%s[%s]'%('nx%dp'%dim,totblock) for dim in range(inp.ndim))),lend[lang])]


    if not inp.MB:
      inputs = inputs + ['int %s = %d%s'%(inp.block, inp.nblocks-1,lend[lang])]
      inputs = inputs + ['\n'.join(list('%s = ;'%inp.grid[dim+1] for dim in range(inp.ndim)))]
    else:
      inputs = inputs + ['%s Write the block dimensions here'%(line_comment[lang])]
      inputs = inputs + ['\n\n']
      inputs = inputs + ['%s Writing the block dimensions  ends here'%(line_comment[lang])]
    #inputs

    # Declare Constants in OPS format
    out = out + ['ops_init(argc,argv,1)%s'%lend[lang]]
    for con in inp.const:
      if isinstance(con, Symbol):
        inputs = inputs + ['%s = %s'%(con, lend[lang])]
      elif isinstance(con, Indexed):
        tot = 0
        for inde in con.shape:
          tot = tot + inde
        for no in range(tot-1):
          inputs = inputs + ['%s[%d] = %s'%(con.base,no, lend[lang])]
    for con in inp.const:
      if con.is_Symbol:
        if con.is_integer:
          dtype = 'int'
        else:
          dtype = 'double'
        out = out + ['ops_decl_const(\"%s\" , 1, \"%s\", &%s)%s'%(con,dtype,con,lend[lang])]
    # Declare block
    out.append('ops_block *%s = (ops_block *)malloc(%s*sizeof(ops_block*));'%(blkname,totblock))
    out.append('char buf[100];')
    if inp.MB:
      stloop,enloop = loop([inp.block],alg)
    else:
      stloop = ['\n']; enloop = ['\n']
    out = out + stloop
    out.append('sprintf(buf,\"%s[%%d]\",%s);'%(blkname,inp.block))
    out.append('%s[%s] = ops_decl_block(%d,buf);'%(blkname,inp.block,inp.ndim))

    out = out + enloop


    # Declare data files
    for da in inp.dats:
      out.append('ops_dat *%s = (ops_dat *)malloc(%s*sizeof(ops_dat*));'%(da.base,totblock))
    out =  out + ['int d_p[%s]   = {%s};'%(inp.ndim,','.join(list('%d'%(inp.halos[dim]) for dim in range(inp.ndim))))]
    out =  out + ['int d_m[%d]   = {%s}; '%(inp.ndim,','.join(list('-%d'%(inp.halos[dim]) for dim in range(inp.ndim))))]
    out = out + ['int base[%d]  = {%s};'%(inp.ndim,','.join(list('%d'%(0) for dim in range(inp.ndim))))]
    out = out + ['double* temp = NULL;']
    out = out + stloop

    size = ','.join(list('%s'%(ind[dim].upper -ind[dim].lower + 1) for dim in range(inp.ndim)))
    out = out + ['int size[%d] = {%s};'%(inp.ndim,size)]
    for da in inp.dats:
      out.append('sprintf(buf,\"%s[%%d]\",%s);'%(da.base,inp.block))
      out = out + ['%s[%s] = ops_decl_dat(%s[%s], 1, size, base, d_m, d_p, temp, "double", buf);'%(da.base,inp.block,blkname,inp.block)]
    out = out + enloop

    for key,value in inp.stencils.iteritems():
      npts = divmod(len(key.split(',')),inp.ndim)
      if npts[1] == 0:
        pass
      else:
        raise ValueError('number of points are not a multiple of dimensions')
      tname = 'sten_%s'%value
      out = out + ['int %s[] = {%s}%s'%(tname,key,lend[lang])]
      sname = 'ops_stencil %s = ops_decl_stencil(%d, %d, %s, \"%s\")%s'
      out = out + [sname%(value,inp.ndim, npts[0], tname, tname, lend[lang])]
    out = out + ['\n\n']+inp.bcdecl + ['\n\n'] + ['ops_partition("");']
    out = inputs + out

  elif lang =='F90':
    inputs = []
    # these are the grid dimensions etc..
    inputs = inputs + ['integer :: %s %s'%(','.join(list('%s'%('nx%dp'%dim) for dim in range(inp.ndim))),lend[lang])]
    inputs = inputs + ['\n'.join(list('%s = '%inp.grid[dim+1] for dim in range(inp.ndim)))]
    for con in inp.const:
      if isinstance(con, Symbol):
        inputs = inputs + ['%s = %s'%(con, lend[lang])]
      elif isinstance(con, Indexed):
        tot = 0
        for inde in con.shape:
          tot = tot + inde
        for no in range(tot):
          inputs = inputs + ['%s(%d) = %s'%(con.base,no, lend[lang])]
    out = []
    dimen = ','.join(list(':' for dim in range(inp.ndim)))

    for da in inp.dats:
      out.append('%s'%(da.base))
    out = 'real(8), allocatable, dimension (%s) :: '%(dimen) + ', '.join(out)
    out = ' &\n'.join(textwrap.wrap(out,width=70,break_long_words=False))
    inp.module.append(out)
    out = []
    # allocate stuff
    ind = inp.blockdims
    sz = []
    for dim in range(inp.ndim):
      l = ind[dim].lower -inp.halos[dim]
      u = ind[dim].upper + inp.halos[dim]
      s = str(l) +':' + str(u)
      sz.append(s)
    sz = ', '.join(sz)
    for da in inp.dats:
      out = out + ['allocate (%s(%s))'%(da.base,sz)]
    out = inputs + out
    out = list(' '*6+ou for ou in out)
    #out = '\n'.join(out)
    #print(out)
    # TODO write this to take care of long allocations
  else:
    raise ValueError('Implement %s in the Definitions and declaration'%alg.lang)


  return out

def footer_code(inp,alg):
  out = []
  lang =alg.lang
  if alg.lang == 'OPSC':
    out = out + ['ops_printf(\" finished running the code\\n\");', 'ops_exit();', '}']
  elif lang == 'F90':
    out = out + ['end program']
  else:
    raise ValueError('Implement %s in the footer code'%alg.lang)
  return out

def indent_code(self, code):
  """Accepts a string of code or a list of code lines"""

  if isinstance(code, string_types):
    code_lines = self.indent_code(code.splitlines(True))
    return ''.join(code_lines)

  tab = "   "
  inc_token = ('{', '(', '{\n', '(\n')
  dec_token = ('}', ')')

  code = [ line.lstrip(' \t') for line in code ]

  increase = [ int(any(map(line.endswith, inc_token))) for line in code ]
  decrease = [ int(any(map(line.startswith, dec_token)))
         for line in code ]

  pretty = []
  level = 0
  for n, line in enumerate(code):
    if line == '' or line == '\n':
      pretty.append(line)
      continue
      level -= decrease[n]
      pretty.append("%s%s" % (tab*level, line))
      level += increase[n]
  return pretty
def bc_periodic_OPSC():


  return
def bcs(inp, alg):
  inp.bcdecl = []
  inp.bccall = []
  inp.bcexind = 0
  finbc = {}
  lang = alg.lang
  if len(alg.bcs) == inp.ndim:
    pass
  elif (len(alg.bcs) > inp.ndim):
    raise ValueError('There are more boundary conditions than the number of dimensions')
  elif (len(alg.bcs) < inp.ndim):
    raise ValueError('There are less boundary conditions than the number of dimensions')
  inp.bc_appl = []
  for dim in range(inp.ndim):
    bc = alg.bcs[dim]
    if bc[0] == bc[1] and bc[0] =='periodic' and alg.lang == 'OPSC' and not inp.MB:
      out = []
      iter_size = list( te for te in inp.gridhalo)

      print('periodic bc in x%d- direction'%dim)
      halo = inp.halos[dim] # get the number of halo points
      iter_size[dim] = halo # the no of halos to be transfered
      from_base = list(-halo for te in range(inp.ndim))
      to_base = list(-halo for te in range(inp.ndim))
      # now copy from data at the end of domain to the first
      l = inp.blockdims[dim].lower - halo
      u = inp.blockdims[dim].upper - halo +1
      from_base[dim] = u
      to_base[dim] = l
      iter_size = ','.join(list(str(te) for te in iter_size))
      fro = ','.join(list(str(te) for te in from_base))
      to = ','.join(list(str(te) for te in to_base))
      #dire = ','.join(['1','2'])
      dire = ','.join(str(i) for i in range(1,inp.ndim+1))
      halo_count = 0
      stloop,enloop = loop([inp.block],alg)
      out = out + stloop
      out = out + ['int off = 0;']
      out = out + ['int halo_iter[] = {%s}%s'%(iter_size,lend[alg.lang])]
      out = out + ['int from_base[] = {%s}%s'%(fro,lend[alg.lang])]
      out = out + ['int to_base[] = {%s}%s'%(to,lend[alg.lang])]
      out = out + ['int dir[] = {%s}%s'%(dire,lend[alg.lang])]
      haloname = 'halos_x%d_%s'%(dim,bc[0])
      #halos[off++] = ops_decl_halo(u[i-1+ngrid_x*j], u[i+ngrid_x*j], halo_iter, base_from, base_to, dir, dir);
      halo_decl = '%s[off++] = ops_decl_halo(%%s[%%s],%%s[%%s],halo_iter, from_base, to_base, dir, dir);'%haloname
      for con in inp.conser:
        out = out + [halo_decl%(con.base.label,inp.block,con.base.label,inp.block)]
        halo_count = halo_count +1
      # copy the first n terms to the Rightend halos
      l = inp.blockdims[dim].lower
      u = inp.blockdims[dim].upper +1
      from_base[dim] = l; to_base[dim] = u
      fro = ','.join(list(str(te) for te in from_base))
      to = ','.join(list(str(te) for te in to_base))
      out = out + ['from_base[%d] = %s%s'%(dim,l,lend[alg.lang])]
      out = out + ['to_base[%d] = %s%s'%(dim,u,lend[alg.lang])]
      for con in inp.conser:
        out = out + [halo_decl%(con.base.label,inp.block,con.base.label,inp.block)]
        halo_count = halo_count +1
      out = ['ops_halo *%s = (ops_halo *)malloc(%d*sizeof(ops_halo *));'%(haloname,halo_count)] + out
      out = out + enloop
      out = out + ['ops_halo_group grp_%s = ops_decl_halo_group(%s,%s);'%(haloname,halo_count,haloname)]
      inp.bcdecl.append([out])
      inp.bccall.append('ops_halo_transfer(grp_%s);'%haloname)
    else:
      for loc_bound,bound in enumerate(bc):
        out = []
        halo_count = 0
        boundary = 'Bc_x%d_%%s_%s'%(dim,bound)
        iter_size = list( te for te in inp.gridhalo)
        halo = inp.halos[dim]
        from_base = list(-halo for te in range(inp.ndim))
        to_base = list(-halo for te in range(inp.ndim))
        if loc_bound == 0:
          bcloc = boundary%'min'
          val = inp.blockdims[dim].lower ; dire =-1
        elif loc_bound == 1:
          bcloc = boundary%'max'
          val = inp.blockdims[dim].upper; dire =1
        else:
          raise ValueError('undefined bc')
        stloop,enloop = loop([inp.block],alg)
        out = out + ['%s\n%s Boundary condition %s\n%s'%(line_comment[lang],line_comment[lang],bcloc,line_comment[lang])]
        out = out + stloop
        te1 =  ','.join(list(str(te) for te in iter_size))
        out = out + ['int off = 0;']
        out = out + ['int halo_iter[] = {%s}%s'%(te1,lend[lang])]
        te1 =  ','.join(list(str(te) for te in from_base))
        out = out + ['int from_base[] = {%s}%s'%(te1,lend[lang])]
        te1 =  ','.join(list(str(te) for te in to_base))
        out = out +['int to_base[] = {%s}%s'%(te1,lend[lang])]
        out = out +['int dir[] = {%s}%s'%(','.join(str(i) for i in range(1,inp.ndim+1)),lend[lang])]
        halo_decl = '%s[off++] = ops_decl_halo(%%s[%%s],%%s[%%s],halo_iter, from_base, to_base, dir, dir);'%bcloc
        print(bcloc)
        if bound == 'symmetry':
          raise ValueError('Implementation %s bc'%bound)
        elif bound == 'zero_grad':
          # change the values depending on the direction
          iter_size[dim] = 1
          structu = '%s[%d] = %s%s'
          out = out + [structu%('halo_iter',dim,iter_size[dim],lend[lang])]
          print(bound,halo,val)
          print('from is', val - dire*1)
          for d in range(0, halo+1):
            from_base[dim] = val - dire*1
            to_base[dim] = val + dire*d
            out = out + [structu%('from_base',dim,from_base[dim],lend[lang])]
            out = out + [structu%('to_base',dim,to_base[dim],lend[lang])]
            for con in inp.conser:
              out = out + [halo_decl%(con.base.label,inp.block,con.base.label,inp.block)]
              halo_count = halo_count +1

        else:
          raise ValueError('dont know implementation of the bc')

        out = ['ops_halo *%s = (ops_halo *)malloc(%d*sizeof(ops_halo *));'%(bcloc,halo_count)] + out
        out = out + enloop
        out = out + ['ops_halo_group grp_%s = ops_decl_halo_group(%s,%s);'%(bcloc,halo_count,bcloc)]
        inp.bcdecl.append([out])
        inp.bccall.append('ops_halo_transfer(grp_%s);'%bcloc)

  return


def after_time(inp, alg):
  out = []
  lang = alg.lang
  if alg.lang =='OPSC':
    writing = 'ops_print_dat_to_txtfile(%s[%s], "%s_final.dat");'
    for con in inp.conser:
      out = out + [writing%(con.base.label,inp.block,con.base.label)]
  elif lang =='F90':
    out = out + ['%swrite data file output here'%line_comment[lang]]
  else:
    raise ValueError('Implement output writing for lanuage %s'%alg.lang)
  return out
import itertools

def write_final_code(template, codes, main, routines, lang):
  import collections
  od = collections.OrderedDict(sorted(template.items()))
  main_code = []
  for key,val in od.iteritems():
    #good go ahead with writing code
    if codes.get(val):
      va = codes.get(val)
      va=  flatten(va)
      main_code = main_code + va + ['\n']
    else:
      main_code = main_code + ['%s There is no code provided for %s part in the algorithm\n'%(line_comment[lang], val)] + ['\n']

  main_code = flatten(main_code)
  main.write('\n'.join(main_code))

  if lang == 'OPSC':
    temp = ['#ifndef kernels_KERNEL_H \n#define kernels_KERNEL_H\n'] #+ codes['kernels'] + ['\n#endif']
    routines.write('\n'.join(temp))
    temp = flatten(codes['kernels'])
    temp = flatten(temp)
    routines.write('\n'.join(temp))
    routines.write('\n#endif')
  elif lang == 'F90':
    temp = flatten(codes['kernels'])
    #temp = flatten(temp)
    routines.write('\n'.join(temp))
  else:
    raise ValueError('Implement %s in write_final_code'%lang)

  return
def write_module(mods, modfile):
  temp = []
  pprint(mods)
  temp = temp + ['module param_mod'] + mods + ['end module param_mod']
  modfile.write('\n'.join(temp))

  return