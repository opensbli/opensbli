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
from string import Template

import logging
LOG = logging.getLogger(__name__)
BUILD_DIR = os.getcwd()


class GenerateCode():
    '''
    This System generate the code for solving the equations
    '''
    def AlgorithmMismatch(Exception):
        return
    def stencil_to_string(self, stencil):
        string = ','.join([str(s) for s in stencil])
        return string

    def __init__(self, grid, spatial_solution, temporal_soln, boundary, Ics, IO, simulation_parameters, Diagnostics = None):
        '''
        The default code generation is OPSC
        '''
        self.check_consistency(grid, spatial_solution, temporal_soln, boundary, Ics,IO)
        self.stencils = {}
        assumptions = {'precission':'double', 'iter_range':0, 'stencil_name':'stencil%d'}
        assumptions['stencil_index']= 0
        assumptions['exchange_self']= 0
        simulation_params = {'niter': 1000}
        language = OPSC(grid.shape)
        code_template,code_dictionary  = self.template(language,assumptions)
        # Process the main time loop
        time_calls, time_kernels, assumptions = self.get_inner_timeloop_code(language, assumptions)
        code_dictionary['time_calls'] = '\n'.join(time_calls)

        time_start_call, time_start_kernel = self.time_start(language,assumptions)
        code_dictionary['time_start_call'] = '\n'.join(time_start_call)

        time_end_call, time_end_kernel = self.time_end(language,assumptions)
        code_dictionary['time_end_call'] = '\n'.join(time_end_call)

        # Process the Boundary conditions
        boundary_calls, boundary_kernels = self.get_bc_code(language,assumptions)
        code_dictionary['bccall'] = '\n'.join(boundary_calls)

        # Process the initialization
        init_calls, init_kernels = self.get_initial_condition(language, assumptions)
        code_dictionary['init_call'] = '\n'.join(init_calls)

        # Process IO calls
        io_calls, io_kernels = self.get_IO_code(language, assumptions)
        code_dictionary['io_calls'] = '\n'.join(io_calls)

        # FIXME presently copying the stencils from this to OPSC, stencils should be implicit to OPSC
        language.stencils = self.stencils

        # Get the header code
        code = (language.header(**assumptions))
        code_dictionary['header'] = '\n'.join(code)

        # get the footer code
        code = language.footer()
        code_dictionary['footer'] = '\n'.join(code)
        # get the language declarations, allocations and soon
        code = language.defdec(grid,assumptions)
        code_dictionary['defdec'] = '\n'.join(code)

        code_template = code_template.safe_substitute(code_dictionary)
        main_file = open(BUILD_DIR+'/%s.cpp' % simulation_parameters["name"], 'w')
        main_file.write(code_template)
        main_file.close()
        # write the final code
        return
    def get_IO_code(self, language, assumptions):
        IO_call = []
        IO_kernel = []
        # As of now only arrays at the end of the simulation
        for block_number, temp in enumerate(self.IO):
            IO_call += language.HDF5_array_fileIO(temp)

        return IO_call, IO_kernel
    def time_end(self, language, assumptions):
        tend_call = []
        tend_kernel = []
        for block_number, temp in enumerate(self.temporal_soln):
            blk_code = []
            blk_calls = []
            if temp.end_computations:
                for comp in temp.end_computations:
                    blk_code += language.kernel_computation(comp,['time_end',block_number],**assumptions)
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
                blk_code += language.kernel_computation(comp,['time_start',block_number],**assumptions)
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
    def template(self,language, assumptions):
        '''
        # header contains all the header information for that language
        # defdec contains all the declarations and initializations of
        #
        # This is the template of the code for a n stage Runge-Kutta method
        '''
        code_template = '''$header \n$defdec \n$init_call \n$bccall \n$timeloop \n$time_start_call
        \n$innerloop \n$time_calls \n $bccall \n$end_inner_loop \n$time_end_call \n$time_io \n$end_time_loop
        \n$io_calls \n$footer'''
        templated_code = Template(code_template)
        if self.temporal_soln[0].nstages >1:
            diction = {}
            ns = self.temporal_soln[0].nstages
            var = self.temporal_soln[0].coeff.stage
            diction['innerloop'] = language.loop_open(var,(0,ns))+ '\n'
            diction['end_inner_loop'] = language.loop_close()

        else:
            diction = {'innerloop':'', 'end_inner_loop':''}
        # Process the time loop
        # FIXME to read iteration ranges from the simulation parameters
        var = 'iteration' # name for iteration
        niter = 100 # Number of iterations, this need to be read in from simulation parameters
        diction['timeloop'] = language.loop_open(var,(0,ns)) + '\n'
        diction['end_time_loop'] = language.loop_close()


        return templated_code, diction


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
    def __init__(self,shape, nblocks=None):
        if nblocks == None or nblocks == 1:
            self.MultiBlock = False
            self.nblocks = 1
        else:
            self.MultiBlock = True
            self.nblocks = nblocks
        self.ndim = len(shape)
        # grid based arrays used for declaration and definition in OPSC format
        self.grid_based_arrays = set()
        # The global constants that are to be declared in c
        self.constants = set()
        # OPS constants, These are the constants in the above list to be defined into OPS format
        self.ops_constant = set()
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
        if computation.has_Idx:
            header += [OPSC.ops_header['Idx']%('idx') ]
        header = [name + ' , '.join(header) + OPSC.close_parentheses ]
        header += [OPSC.open_brace]
        code =  header
        ops_accs = self.get_OPS_ACCESS_number(computation)
        for eq in computation.equations:
            code += [ccode(eq,ops_accs)+ OPSC.end_of_statement]
        code += [OPSC.close_brace] + ['\n']
        self.update_definitions(computation)
        return code
    def update_definitions(self, computation):
        arrays = set([inp for inp in computation.inputs.keys() if inp.is_grid] + \
            [inp for inp in computation.outputs.keys() if inp.is_grid ] + \
                [inp for inp in computation.inputoutput.keys() if inp.is_grid])
        constant_arrays = set([inp for inp in computation.inputs.keys() if not inp.is_grid] + \
            [inp for inp in computation.outputs.keys() if not inp.is_grid ] + \
                [inp for inp in computation.inputoutput.keys() if not inp.is_grid])
        constants = set(computation.constants)
        #constants need to think??
        self.grid_based_arrays = self.grid_based_arrays.union(arrays)
        self.constants = self.constants.union(constant_arrays).union(constants)
        return
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
        if computation.has_Idx:
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
        code = ['%s Boundary condition exchange code'%OPSC.line_comment]
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
        call = ['%s Boundary condition exchange calls'%OPSC.line_comment,'ops_halo_transfer(%s)%s'%(name,OPSC.end_of_statement)]
        return call, code
    def loop_open(self, var, range_of_loop):
        return 'for (int %s=%d; %s<%d, %s++)%s'%(var, range_of_loop[0], var, range_of_loop[1], var,\
            OPSC.open_brace)
    def loop_close(self):
        return OPSC.close_brace
    def header(self, **assumptions):
        precision = assumptions['precission']
        code = []
        code += ['#include <stdlib.h>']
        code += ['#include <string.h>']
        code += ['#include <math.h>']
        code += ['%s Global Constants in the equations are' % OPSC.line_comment]
        for con in self.constants:
            if isinstance(con, IndexedBase):
                code += ['%s %s[%d]%s'%(precision, con, con.ranges, OPSC.end_of_statement)]
            else:
                code += ['%s %s%s'%(precision, con, OPSC.end_of_statement)]
        # Include constant declaration
        code += ['// OPS header file']
        code += ['#define OPS_%sD' % self.ndim]
        code += ['#include "ops_seq.h"']
        # Include the kernel file names
        code += ['#include "auto_kernel.h"']
        code += ['%s main program start' % OPSC.line_comment]
        code += ['int main (int argc, char **argv) {']
        return code
    def ops_init(self, diagnostics_level=None):
        '''
        the default diagnostics level in 1 which is the best performance
        refer to ops user manual
        '''
        out = ['%s Initializing OPS '%OPSC.line_comment]
        if diagnostics_level:
            self.ops_diagnostics = True
            return out + ['ops_init(argc,argv,%d)%s'%(diagnostics_level, OPSC.end_of_statement)]
        else:
            self.ops_diagnostics = False
            return out + ['ops_init(argc,argv,%d)%s'%(1, OPSC.end_of_statement)]
    def ops_diagnostics(self):
        '''
        untested OPS diagnostics output need to check if it gives the result or not
        '''
        if self.ops_diagnostics:
            return ['ops diagnostic output()']
        else:
            return []
        return
    def ops_exit(self):
        return ['%s Exit OPS '%OPSC.line_comment,'ops_exit()']
    def ops_partition(self):
        return ['%s Init OPS partition'%OPSC.line_comment,'ops_partition(\" \")']
    def ops_timers(self):
        st = ["cpu_start", "elapsed_start"]
        en = ["cpu_end", "elapsed_end"]
        timer_start = ["double %s, %s%s"%(st[0],st[1],OPSC.end_of_statement)]\
            + ["ops_timers(&%s, &%s)%s"%(st[0], st[1], OPSC.end_of_statement)]
        timer_end = ["double %s, %s%s"%(en[0],en[1],OPSC.end_of_statement)]\
            + ["ops_timers(&%s, &%s)%s"%(en[0], en[1], OPSC.end_of_statement)]
        timing_eval = self.ops_print_timings(st, en)
        return timer_start, timer_end, timing_eval
    def ops_print_timings(self, st, en):
        code = []
        code += ["ops_printf(\"\\nTimings are:\\n\")%s"%OPSC.end_of_statement]
        code += ["ops_printf(\"-----------------------------------------\\n\")%s"%OPSC.end_of_statement]
        code += ["ops_printf(\"Total Wall time %%lf\\n\",%s-%s)%s"%(en[1], st[1],OPSC.end_of_statement)]
        return
    def define_block(self):
        code = ['%s Defining block in OPS Format'%(OPSC.line_comment)]
        if not self.MultiBlock:
            # no dynamic memory allocation required
            code += ['ops_block  %s%s'\
            % (OPSC.block_name, OPSC.end_of_statement)]
        else:
            code += ['ops_block *%s = (ops_block *)malloc(%s*sizeof(ops_block*))%s'\
                % (OPSC.block_name, self.nblocks, OPSC.end_of_statement )]
        #print('\n'.join(code))
        return code
    def initialize_block(self):
        code = ['%s Initializing block in OPS Format'%(OPSC.line_comment)]
        if not self.MultiBlock:
            code += ['%s = ops_decl_block(%d, \"%s\")%s'\
                %(OPSC.block_name,self.ndim, OPSC.block_name,OPSC.end_of_statement)]
        else:
            raise NotImplementedError("Multi block is not implemented")
        #print('\n'.join(code))
        return code
    def define_dat(self):
        code = ['%s Define data files'%(OPSC.line_comment)]
        if not self.MultiBlock:
            def_format = 'ops_dat %%s%s'% OPSC.end_of_statement
            code += [def_format%arr for arr in self.grid_based_arrays]
        else:
            raise NotImplementedError("Multi block is not implemented")
        #print('\n'.join(code))
        return code
    def initialize_dat(self, grid, assumptions):
        code = ['%s Initialize/ Allocate data files'%(OPSC.line_comment)]
        precision = assumptions['precission']
        dtype_int = 'int'
        if not self.MultiBlock:
            code += [self.helper_array(dtype_int, 'halo_p', [halo[1] for halo in grid.halos])]
            code += [self.helper_array(dtype_int, 'halo_m', [halo[0] for halo in grid.halos])]
            code += [self.helper_array(dtype_int, 'size', grid.shape)]
            code += [self.helper_array(dtype_int, 'base', [0 for g in grid.shape])]
            code += ['%s* val= Null;'%(precision)]
            init_format = '%%s = ops_decl_dat(%s, 1, size, base, halo_m, halo_p, val, \"%%s\", \"%%s\")%s'\
            % (OPSC.block_name, OPSC.end_of_statement)
            inits = [init_format%(arr, precision, arr) for arr in self.grid_based_arrays]
            code = code + inits
        else:
            raise NotImplementedError("Multi block is not implemented")
        #print('\n'.join(code))
        return code
    def helper_array(self, dtype, name, values):
        ''' Helper function to declare inline arrays in OPSC
        dtype: data type
        name: name of the array
        size: size of the array
        vals: list of values
        '''
        return '%s %s[] = {%s}%s'%(dtype, name, ', '.join([str(s) for s in values]),\
             OPSC.end_of_statement)
    def declare_stencils(self):
        '''
        This declares all the stencils used in the code.
        We donot differentiate between the stencils for each block.
        returns the code
        '''
        code = ['%s Declare all the stencils used '%(OPSC.line_comment)]
        dtype_int = 'int'
        sten_format = 'ops_stencil %%s = ops_decl_stencil(%%d,%%d,%%s,\"%%s\")%s'%(OPSC.end_of_statement)
        for key, value in self.stencils.iteritems():
            count = len(key.split(','))/ self.ndim
            # value is the name in the stencils format
            code += [self.helper_array(dtype_int, value, [key])]
            code += [sten_format%(value, self.ndim, count, value, key)]
        return code
    def HDF5_array_fileIO(self,instance):
        code = []
        # to do generate file name automatically
        block_to_hdf5 = ["ops_fetch_block_hdf5_file(%s, \"state.h5\")%s" % (OPSC.block_name, OPSC.end_of_statement)]
        code += block_to_hdf5
        # Then write out each field.
        for c in instance.save_arrays:
            convervative_variable_to_hdf5 = ["ops_fetch_dat_hdf5_file(%s, \"state.h5\")%s" \
                % (c, OPSC.end_of_statement)]
            code += convervative_variable_to_hdf5
        return code
    def defdec(self, grid, assumptions):
        defdec = []
        defdec += self.ops_init() + ['\n'] + self.define_block()+ ['\n'] + self.initialize_block()+ ['\n'] + self.define_dat()+ ['\n']
        defdec += self.initialize_dat(grid, assumptions) + ['\n']
        defdec += self.declare_stencils() + ['\n']
        return defdec
    def footer(self):
        '''
        This writes out the footer code in OPSC this is a call to OPS_exit and closing the main loop
        '''
        code = self.ops_exit()
        code += [OPSC.close_brace]
        return code
