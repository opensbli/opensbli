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
from sympy.printing.ccode import CCodePrinter
from sympy.parsing.sympy_parser import parse_expr
import re
import os
# AutoFD functions
from .codegen_utils import END_OF_STATEMENT_DELIMITER

import logging
LOG = logging.getLogger(__name__)
BUILD_DIR = os.getcwd()

COMMENT_DELIMITER = {"OPSC": "//", "F90": "!"}
END_OF_STATEMENT_DELIMITER = {"OPSC": ";", "F90": ""}

class GenerateCode():
    '''
    This System generate the code for solving the equations
    '''
    def AlgorithmMismatch(Exception):
        return
    def stencil_to_string(self, stencil):
        string = ','.join([str(s) for s in stencil])
        return string

    def __init__(self, grid, spatial_solution, temporal_soln, boundary, Ics, IO, simulation_parameters,code = None):
        '''
        The default code generation is OPSC
        '''
        self.check_consistency(grid, spatial_solution, temporal_soln, boundary, Ics,IO)
        self.gridsizes = []
        self.constants = []
        self.arrays = []
        self.get_arrays_constants()
        self.stencils = {}
        assumptions = {'precission':'double', 'iter_range':0, 'stencil_name':'stencil%d'}
        assumptions['stencil_index']= 0
        assumptions['exchange_self']= 0
        simulation_params = {'niter': 1000}
        language = OPSC(grid.shape)
        # inner time loop calls and kernels 
        time_calls, time_kernels, assumptions = self.get_inner_timeloop_code(language, assumptions)
        # get the Bc's code
        boundary_calls, boundary_kernels = self.get_bc_code(language,assumptions)
        # get the Initialization code and kernels
        init_calls, init_kernels = self.get_initial_condition(language, assumptions)
        if temporal_soln.nstages >1:
            time_start_call, time_start_kernel = self.time_start(language,assumptions)
            time_end_call, time_end_kernel = self.time_end(language,assumptions)
        return
    def time_end(self, language, assumptions):
        tend_call = []
        tend_kernel = []
        for block_number, temp in enumerate(self.temporal_soln):
            blk_code = []
            blk_calls = []
            for comp in temp.start_computations:
                blk_code += language.kernel_computation(comp,['temp_start',block_number],**assumptions)
                stencils = language.get_stencils(comp)
                for sten,value in stencils.iteritems():
                    key = self.stencil_to_string(value)
                    if key not in self.stencils.keys():
                        self.stencils[key] = stencil_name%(stencil_index);
                        stencil_index = stencil_index+1
                    # update the stencil name in stencils
                    stencils[sten] = self.stencils[key]
                call, assumptions = language.kernel_call(comp,stencils, **assumptions)
                blk_calls += call
            tend_kernel += blk_code
            tend_call += blk_calls
        return tend_call, tend_kernel
    def time_start(self, language, assumptions):
        tstart_call = []
        tstart_kernel = []
        for block_number, temp in enumerate(self.temporal_soln):
            blk_code = []
            blk_calls = []
            for comp in temp.start_computations:
                blk_code += language.kernel_computation(comp,['temp_start',block_number],**assumptions)
                stencils = language.get_stencils(comp)
                for sten,value in stencils.iteritems():
                    key = self.stencil_to_string(value)
                    if key not in self.stencils.keys():
                        self.stencils[key] = stencil_name%(stencil_index);
                        stencil_index = stencil_index+1
                    # update the stencil name in stencils
                    stencils[sten] = self.stencils[key]
                call, assumptions = language.kernel_call(comp,stencils, **assumptions)
                blk_calls += call
            tstart_kernel += blk_code
            tstart_call += blk_calls
        return tstart_call, tstart_kernel
    def get_initial_condition(self, language, assumptions):
        '''
        All works fine for single block, need to do something for multiblock sense
        '''
        init_kernel = []
        init_call = []
        for block_number,ic in enumerate(self.Ics):
            blk_code = []
            blk_calls = []
            for comp_number,comp in enumerate(ic.computations):
                blk_code += language.kernel_computation(comp,['initial',block_number],**assumptions)
                stencils = language.get_stencils(comp)
                for sten,value in stencils.iteritems():
                    key = self.stencil_to_string(value)
                    if key not in self.stencils.keys():
                        self.stencils[key] = stencil_name%(stencil_index);
                        stencil_index = stencil_index+1
                    # update the stencil name in stencils
                    stencils[sten] = self.stencils[key]
                call, assumptions = language.kernel_call(comp,stencils, **assumptions)
                blk_calls += call
            init_kernel += blk_code
            init_call += blk_calls

        return init_call, init_kernel
    def get_bc_code(self, language, assumptions):
        ''' loop over boundary conditions instances'''
        bc_calls = []
        bc_kernels = []
        for block_number in range(len(self.boundary)):
            # loop over boundaries of the instance
            call_block = []
            code_block = []
            ninst = len(self.boundary[block_number].type_of_boundary)
            boundary_inst = self.boundary[block_number]
            for inst in range(ninst):
                if boundary_inst.type_of_boundary[inst] == 'exchange_self':
                    call, code = language.bc_exchange(boundary_inst.transfers[inst], assumptions)
                    call_block += call
                    code_block += code
            bc_calls += call_block
            bc_kernels += code_block
        return bc_calls, bc_kernels
    def get_inner_timeloop_code(self, language, assumptions):
        '''
        This will evaluate all the inner time loop if exists or this is the time loop
        Returns a list of lists depending on the number of blocks
        '''
        stencil_name = assumptions.get('stencil_name')
        stencil_index = assumptions.get('stencil_index')
        inner_tloop_calls = []
        inner_tloop_kernels = []
        for block_number in range(len(self.spatial_solution)):
            allinst = [self.spatial_solution[block_number]] + [self.temporal_soln[block_number]]
            calls = []
            kernels = []
            # get the spatial and temporal solutions
            for soln in allinst:
                for comp_number,comp in enumerate(soln.computations):
                    code = language.kernel_computation(comp,[block_number,comp_number],**assumptions)
                    kernels += code
                    stencils = language.get_stencils(comp)
                    for sten,value in stencils.iteritems():
                        key = self.stencil_to_string(value)
                        if key not in self.stencils.keys():
                            self.stencils[key] = stencil_name%(stencil_index);
                            stencil_index = stencil_index+1
                        # update the stencil name in stencils
                        stencils[sten] = self.stencils[key]
                    call, assumptions = language.kernel_call(comp,stencils, **assumptions)
                    calls += call
            inner_tloop_kernels += kernels
            inner_tloop_calls += calls
        assumptions['stencil_index'] = stencil_index

        return inner_tloop_calls, inner_tloop_kernels, assumptions
    def get_arrays_constants(self):
        allinst = self.spatial_solution+self.temporal_soln+self.boundary+self.Ics
        arrays = set()
        const = set()
        stencils = set()
        for inst in allinst:
            for comp in inst.computations:
                if comp != None:
                    arrays = arrays.union(set(comp.inputs.keys())).union(set(comp.outputs.keys())).union(set(comp.inputoutput.keys()))
                    #pprint(["COMPUTATION IS ",comp.equations, comp.computation_type])
                    const = const.union(set(comp.constants))
        # The arrays are the ones that are to be declared on the grid
        self.arrays = [arr for arr in arrays if not isinstance(arr,str) if arr.is_grid]
        self.constants = [con for con in const] + [arr for arr in arrays if not isinstance(arr,str) if not arr.is_grid]
        return

    def check_consistency(self,  grid, spatial_solution, temporal_soln, boundary, Ics, IO):
        self.grid = self.listvar(grid)
        length = len(self.grid)

        self.spatial_solution = self.listvar(spatial_solution);
        if len(self.spatial_solution) != length:
            raise AlgorithmMismatch("The length of spatial solution doesnot match the grid")

        self.temporal_soln = self.listvar(temporal_soln);
        if len(self.temporal_soln) != length:
            raise AlgorithmMismatch("The length of temporal solution doesnot match the grid")

        self.boundary = self.listvar(boundary);
        if len(self.boundary) != length:
            raise AlgorithmMismatch("The length of boundary doesnot match the grid")

        self.Ics = self.listvar(Ics);
        if len(self.Ics) != length:
            raise AlgorithmMismatch("The length of Ics doesnot match the grid")

        self.IO = self.listvar(IO);
        if len(self.IO) != length:
            raise AlgorithmMismatch("The length of IO doesnot match the grid")

        return
    def listvar(self, var):
        if isinstance(var, list):
            return var
        else:
            return [var]

class OPSCCodePrinter(CCodePrinter):
    def __init__(self, Indexed_accs):
        settings = {}
        CCodePrinter.__init__(self, settings)
        # Indexed access numbers are required in dictionary
        self.Indexed_accs = Indexed_accs
    def _print_Rational(self, expr):
        p, q = int(expr.p), int(expr.q)
        return '%d.0/%d.0' %(p,q)
    def _print_Indexed(self, expr):
        # find the symbols in the indices of the expression
        syms = flatten([list(at.atoms(Symbol)) for at in expr.indices])
        # Replace the symbols in the indices with `zero'
        for x in syms:
            expr = expr.subs({x: 0})
        if self.Indexed_accs[expr.base] != None:
            out = "%s[%s(%s)]"%(self._print(expr.base.label) \
            , self.Indexed_accs[expr.base], ','.join([self._print(ind)  for ind in expr.indices]))
        else:
            out = "%s[%s]"%(self._print(expr.base.label) \
            , ','.join([self._print(ind)  for ind in expr.indices]))
        return out

def ccode(expr, Indexed_accs):
    if isinstance(expr, Eq):
        return OPSCCodePrinter(Indexed_accs).doprint(expr.lhs) \
        + ' = ' + OPSCCodePrinter(Indexed_accs).doprint(expr.rhs)
    return OPSCCodePrinter(Indexed_accs).doprint(expr)

class OPSC(object):
    '''
    '''
    # OPS Access types, used for kernel call
    ops_access = {'inputs':'OPS_READ', 'outputs':'OPS_WRITE', 'inputoutput':'OPS_RW'}
    # OPS kenel header
    ops_header = {'inputs':'const %s *%s', 'outputs':'%s *%s', 'inputoutput':'%s *%s', 'Idx':'const int *%s'}
    # Line comments
    line_comment = "//"
    # end of line delimiter
    end_of_statement = ";"
    # Commonly used brackets
    open_brace = "{"; close_brace = "}"
    open_parentheses = "("; close_parentheses = ")"
    block_name = 'auto_block_OPSC'
    def __init__(self,shape):
        self.ndim = len(shape)
        return
    def kernel_computation(self, computation, number, **assumptions):
        precission = assumptions['precission']
        header = []
        if computation.name == None:
            computation.name = self.get_kernel_name(number)
        name = computation.name + OPSC.open_parentheses
        # process inputs
        gridbased = ([OPSC.ops_header['inputs']%(precission,inp) for inp in computation.inputs.keys() if inp.is_grid] + \
            [OPSC.ops_header['outputs']%(precission,inp) for inp in computation.outputs.keys() if inp.is_grid ] + \
                [OPSC.ops_header['inputoutput']%(precission,inp) for inp in computation.inputoutput.keys() if inp.is_grid])
        # nongrid based inputs are
        nongrid = ([OPSC.ops_header['inputs']%(precission,inp) for inp in computation.inputs.keys() if not inp.is_grid] + \
            [OPSC.ops_header['outputs']%(precission,inp) for inp in computation.outputs.keys() if not inp.is_grid ] + \
                [OPSC.ops_header['inputoutput']%(precission,inp) for inp in computation.inputoutput.keys() if not inp.is_grid])
        header = gridbased + nongrid
        if computation.Idx:
            header += [OPSC.ops_header['Idx']%('idx') ]
        header = [name + ' , '.join(header) + OPSC.close_parentheses ]
        header += [OPSC.open_brace]
        code =  header
        ops_accs = self.get_OPS_ACCESS_number(computation)
        for eq in computation.equations:
            code += [ccode(eq,ops_accs)+ OPSC.end_of_statement]
        code += [OPSC.close_brace] + ['\n']
        return code
    def get_OPS_ACCESS_number(self, computation):
        '''
        Returns a dictionary of OPS_ACC's
        '''
        ops_accs = {}
        allidbs = list(computation.inputs.keys()) +  list(computation.outputs.keys()) + list(computation.inputoutput.keys())
        gridbased = [al for al in allidbs if al.is_grid]
        # all grid based ops_accs
        for no,inp in enumerate(gridbased):
            ops_accs[inp] = 'OPS_ACC%d'%no
        # non grid based stuff are
        nongrid = set(allidbs).difference(set(gridbased))
        for no,inp in enumerate(nongrid):
            ops_accs[inp] = None


        return ops_accs
    def get_kernel_name(self, number):
        return 'void kernel_block_%s_computation_%s' %(str(number[0]), str(number[1]))
    def get_stencils(self, computation):
        allidb = list(computation.inputs.keys()) +  list(computation.outputs.keys()) + list(computation.inputoutput.keys())
        stencils = {}
        dicts = [computation.inputs, computation.outputs,computation.inputoutput]
        for d in dicts:
            for key, value in d.iteritems():
                if key.is_grid:
                    stencils[key] = self.relative_stencil(value)
        # REturn stencil as a string
        return stencils
    def relative_stencil(self, value):
        '''
        This returns the relative stencil wrt the grid location
        i.e. grid indices eg(i0,i1,i2) are replaced with (0,0,0)
        TODO Need to check if OPS also requires the grid location
        '''
        if isinstance(value,list):
            pass
        else:
            value = [value]
        retun_val = []
        for va in value:
            out = []
            for number, v in enumerate(va):
                outv = v
                for a in v.atoms(Symbol):
                    outv = outv.subs(a,0)
                out.append(outv)
            retun_val.append(out)
        retun_val = self.sort_stencil_indices(retun_val)

        return retun_val
    def sort_stencil_indices(self, indexes):
        '''
        helper function for relative_stencil, sorts the relative stencil in
        the directions
        '''
        if len(indexes[0]) > 1:
            for dim in range(len(indexes[0])):
                indexes = sorted(indexes, key=lambda indexes: indexes[dim])
            temp = flatten(list(list(t) for t in indexes))
        else:
            indexes = [sorted(indexes)]
            temp = flatten(list(t) for t in indexes)
        return temp
    def kernel_call(self, computation, stencils, **assumptions):
        precission = assumptions['precission']
        iterrange_ind = assumptions['iter_range']
        iterrange = 'iter_range%d'%iterrange_ind
        kercall = []
        # range of iterations
        range_main = 'int %s[] = {%s}'%(iterrange, ', '.join([str(r) for ran in computation.ranges for r in ran]))
        range_main = range_main +  OPSC.end_of_statement
        kercall += ['ops_par_loop(%s, \"%s\", %s, %s, %s' % (computation.name,\
            computation.computation_type,OPSC.block_name, self.ndim, iterrange)]
        assumptions['iter_range'] = iterrange_ind+1
        # do the inputs first gridbased
        gridbased = [self.ops_argument_call(inp, stencils[inp],precission, OPSC.ops_access['inputs'])\
            for inp in computation.inputs.keys() if inp.is_grid ] + \
                [self.ops_argument_call(inp, stencils[inp],precission, OPSC.ops_access['outputs'])\
                    for inp in computation.outputs.keys() if inp.is_grid ] + \
                        [self.ops_argument_call(inp, stencils[inp],precission, OPSC.ops_access['inputoutput'])\
                            for inp in computation.inputoutput.keys() if inp.is_grid ]
        # Globals
        nongrid = [self.ops_global_call(inp, value, precission, OPSC.ops_access['inputs'])\
            for inp,value in computation.inputs.iteritems() if not inp.is_grid ] + \
                [self.ops_global_call(inp, value, precission, OPSC.ops_access['outputs'])\
                    for inp, value in computation.outputs.iteritems() if not inp.is_grid ] + \
                        [self.ops_global_call(inp, value, precission, OPSC.ops_access['inputoutput'])\
                            for inp, value in computation.inputoutput.iteritems() if not inp.is_grid ]
        if computation.Idx:
            nongrid += [self.grid_index_call()]
        kercall = kercall + gridbased + nongrid
        call = [k+',' for k in kercall[:-1]]
        call = [range_main] + call + [kercall[-1] + OPSC.close_parentheses + OPSC.end_of_statement] + ['\n']
        return call,assumptions
    def grid_index_call(self):
        return 'ops_arg_idx()'
    def ops_global_call(self, array, indices, precission, access_type):
        arr = array[tuple(indices[0])]
        template = 'ops_arg_gbl(&%s, %d, \"%s\", %s)'
        return template%(arr,1, precission, access_type)
    def ops_argument_call(self, array, stencil, precission, access_type):
        template = 'ops_arg_dat(%s, %d, %s, \"%s\", %s)'
        return template%(array,1,stencil, precission, access_type)
    def bc_exchange(self, instance, assumptions):
        off = 0; halo = 'halo'
        #name of the halo exchange
        name = 'halo_exchange_self%d'%(assumptions['exchange_self'])
        # update exchange self index
        assumptions['exchange_self']= assumptions['exchange_self']+1
        code = []
        code += ['ops_halo_group %s %s'%(name, OPSC.end_of_statement)]
        code += [OPSC.open_brace]
        code += ['int halo_iter[] = {%s}%s'%(', '.join([str(s) for s in instance.transfer_size]), OPSC.end_of_statement)]
        code += ['int from_base[] = {%s}%s'%(', '.join([str(s) for s in instance.transfer_from]), OPSC.end_of_statement)]
        code += ['int to_base[] = {%s}%s'%(', '.join([str(s) for s in instance.transfer_to]), OPSC.end_of_statement)]
        # dir in OPSC not sure what it is but 1to ndim works
        code += ['int dir[] = {%s}%s'%(', '.join([str(ind+1) for ind in range(len(instance.transfer_to))]), OPSC.end_of_statement)]
        # now process the arrays
        for arr in instance.transfer_arrays:
            code += ['ops_halo %s%d = ops_decl_halo(%s, %s, halo_iter, from_base, to_base, dir, dir)%s'\
                %(halo, off, arr, arr, OPSC.end_of_statement)]
            off = off+1
        code += ['ops_halo grp[] = {%s}%s'%(','.join([str('%s%s'%(halo, of)) for of in range(off)]),OPSC.end_of_statement )]
        code += ['%s = ops_decl_halo_group(%d,grp)'%(name, off)]
        code += [OPSC.close_brace]
        # finished OPS halo exchange, now get the call
        call = ['ops_halo_transfer(%s)%s'%(name,OPSC.end_of_statement)]
        return call, code
    def loop_open(self, var, range_of_loop):
        return 'for (int %s=%d; %s<=%d, %s++)%s'%(var, range_of_loop[0], var, range_of_loop[1], var,\
            OPSC.open_brace)
    def loop_close(self):
        return OPSC.close_brace
    def header(self):
        return