from opensbli.core.grid import Grid
from sympy import pprint, Equality
from opensbli.core.bcs import BoundaryConditionTypes
from opensbli.core.opensbliobjects import ConstantObject, DataObject, DataSetBase, GroupedPiecewise
from opensbli.equation_types.opensbliequations import OpenSBLIEq, ConstituentRelations
from opensbli.equation_types.metric import MetricsEquation
from sympy import flatten, eye
_known_equation_types = (GroupedPiecewise, OpenSBLIEq)
from opensbli.schemes.spatial.scheme import CentralHalos_defdec


class KernelCounter():
    """A Counter for the kernels, this stores the current kernel number for a block,
    and is used to name the kernels
    """

    def __init__(self):
        self.kernel_counter = 0
        self.stored_counter = 0

    @property
    def reset_kernel_counter(self):
        """Resets the kernel counter to zero
        """
        self.kernel_counter = 0
        return

    @property
    def increase_kernel_counter(self):
        """Increases the kernel counter by 1
        """
        self.kernel_counter = self.kernel_counter + 1
        return

    @property
    def store_kernel_counter(self):
        """Stores the current values to a variables
        """
        self.stored_counter = self.kernel_counter
        return

    @property
    def reset_kernel_to_stored(self):
        """Updates the kernel counter to the previously stored value
        """
        self.kernel_counter = self.stored_counter
        return


class SimulationBlock(Grid, KernelCounter, BoundaryConditionTypes):
    """This represents the grid on which the equations, boundary conditions etc are
    set
    """

    def __init__(self, ndim, block_number=None):
        if block_number:
            self.blocknumber = block_number
        else:
            self.blocknumber = 0
        self.ndim = ndim
        # Instantiate the kernel counter
        KernelCounter.__init__(self)
        # Instantiate grid class
        Grid.__init__(self)
        # Empty sets for the boundary conditions, these halo type for the discretisation
        # scheme will be added to these sets depending on the derivatives in the equations
        self.boundary_halos = [[set(), set()] for d in range(self.ndim)]
        # place holders to store various parameters
        self.block_datasets = {}
        self.constants = {}
        self.Rational_constants = {}
        self.block_stencils = {}
        self.InputOutput = []
        self.list_of_equation_classes = []
        return

    @property
    def blockname(self):
        """Name of the block
        """
        return 'opensbliblock%02d' % self.blocknumber

    @property
    def blocknumber(self):
        """Stores the block number as a attribute"""
        return self.__blocknumber

    @blocknumber.setter
    def blocknumber(self, number):
        self.__blocknumber = number
        return

    def set_block_boundaries(self, bclist):
        """Sets the boundary conditions for the block
        """
        self.set_boundary_types(bclist, self)
        return

    def set_block_boundary_halos(self, direction, side, types):
        """
        """
        self.boundary_halos[direction][side].add(types)
        return

    def dataobjects_to_datasets_on_block(self, eqs):
        """
        """
        store_equations = flatten(eqs)[:]
        consts = set()

        for no, eq in enumerate(store_equations):
            if isinstance(eq, _known_equation_types):
                store_equations[no] = eq.convert_to_datasets(self)
                consts = consts.union(eq.atoms(ConstantObject))
            elif isinstance(eq, Equality):
                new_eq = OpenSBLIEq(eq.lhs, eq.rhs)
                consts = consts.union(new_eq.atoms(ConstantObject))
                store_equations[no] = new_eq.convert_to_datasets(self)
            elif isinstance(eq, DataObject):
                store_equations[no] = self.location_dataset(str(eq))
            else:  # Integers and Floats from Eigensystem entering here
                pass
        # Convert all equations into the format of input equations WARNING crude way
        out = []
        out_loc = 0
        for no, e in enumerate(eqs):
            if isinstance(e, list):
                out += [store_equations[out_loc:out_loc+len(e)]]
                out_loc += len(e)
            else:
                out += [store_equations[out_loc]]
                out_loc += 1
        # Add the constants to ConstantsToDeclare
        from .kernel import ConstantsToDeclare as CTD
        for c in consts:
            CTD.add_constant(c)
        return out

    def copy_block_attributes(self, otherclass):
        """Set the attributes blocknumber, name and number of dimensions of the
        block on the other block.
        """
        # TODO V2, why the name for blocknumber is different, we should be consistent
        otherclass.block_number = self.blocknumber
        otherclass.ndim = self.ndim
        otherclass.block_name = self.blockname
        return

    @property
    def known_datasets(self):
        """Collects all the datasets evaulated by each of the equation classes,
        used for sorting the relations for a equation class
        """
        known_dsets = set()
        for eq in self.list_of_equation_classes:
            known_dsets = known_dsets.union(eq.evaluated_datasets)
        for io in self.InputOutput:
            known_dsets = known_dsets.union(io.evaluated_datasets)
        return known_dsets

    def discretise(self):
        """Discretises the equations by recursively calling the ``spatial_discretisation`` function
        of the schemes provided if the scheme is of type spatial. If the scheme is of type temporal
        then it recursively calls the temporal discretisation for each equation class.
        See, schemes folder for further information
        """
        # perform the spatial discretisation of the equations using schemes
        for eq in self.list_of_equation_classes:
            eq.spatial_discretisation(self)
            eq.apply_boundary_conditions(self)

        # get the temporal schemes set on the block and discretise
        temporal = self.get_temporal_schemes
        for t in temporal:
            for eq in self.list_of_equation_classes:
                self.discretisation_schemes[t.name].discretise(eq, self)

        # Update the data sets that are to be read from HDF5
        for io in self.InputOutput:
            io.set_read_from_hdf5_arrays(self)
        return

    def create_datasetbase(self, name):
        # A wrapper to create a dataset base on a block, used in discretisation and
        # while applying boundary conditions
        return DataSetBase(str(name), self.shape, self.blocknumber)

    def location_dataset(self, name):
        """Creates a dataset at the location (base location) of the dataset base, OpenSBLI
        uses relative indexing, see opensbliobjects for further details

        :param name: name of the Dataset, can be type str or opensbli:DataSetBase
        :returns: Current location dataset
        :rtype: DataSet
        """
        if isinstance(name, DataSetBase):
            return name[name.location]
        else:
            base = self.create_datasetbase(name)
            return base[base.location]

    def apply_boundary_conditions(self, arrays):
        """Function that applies the boundary conditions on the arrays provided, called from
        boundary conditions for equation classes
        """
        kernels = []
        for no, b in enumerate(self.boundary_types):
            for side in [0, 1]:
                k = self.apply_bc_direction(no, side, arrays)
                if isinstance(k, list):
                    kernels += k
                else:
                    kernels += [k]
        return kernels

    def apply_bc_direction(self, direction, side, arrays):
        """Apply boundary condition for a direction and side
        """
        kernel = self.boundary_types[direction][side].apply(arrays, self)
        return kernel

    def set_equations(self, list_of_equations):
        """Sets the equation classes to be used on a particular block
        """
        for eq in list_of_equations:
            block_eq = self.dataobjects_to_datasets_on_block(eq.equations)
            eq.equations = block_eq
        self.list_of_equation_classes += list_of_equations
        self.sort_equation_classes
        return

    @property
    def sort_equation_classes(self):
        """Sorts the equation classes so that the constituent relations are first in the list,
        the order of discretisation of the equation classes is set here. Any other dependencies
        should be added here
        """

        for no, eq in enumerate(self.list_of_equation_classes):
            if isinstance(eq, ConstituentRelations):
                cr = self.list_of_equation_classes.pop(no)
                self.list_of_equation_classes.insert(0, cr)

        return

    def set_discretisation_schemes(self, schemes):
        """Sets the discretiation schemes to be used for the equations on the block
        """
        self.discretisation_schemes = schemes
        return

    @property
    def get_constituent_equation_class(self):
        """Returns the classes for constituent relations, there can be many constituent
        relations classes on a particular block
        """
        CR_classes = []
        for sc in self.list_of_equation_classes:
            if isinstance(sc, ConstituentRelations):
                CR_classes += [sc]
        return CR_classes

    @property
    def get_temporal_schemes(self):
        """Returns the temporal schemes set on the block
        """
        temporal = []
        for sc in self.discretisation_schemes:
            if self.discretisation_schemes[sc].schemetype == "Temporal":
                temporal += [self.discretisation_schemes[sc]]
        return temporal

    @property
    def collect_all_spatial_kernels(self):
        """ TODO check if we are using this
        """
        all_kernels = []
        for scheme in self.get_temporal_schemes:
            for key, value in scheme.solution.iteritems():  # These are equation classes
                if key.order >= 0 and key.order < 100:  # Checks if the equation classes are part of the time loop
                    all_kernels += key.all_spatial_kernels
                else:
                    print('NOPE')
        return all_kernels

    def setio(self, list_of_ios):
        """Sets the input and output for the block
        """
        self.add_io(list_of_ios)
        return

    def add_io(self, list_of_ios):
        """adds the io classes
        """
        if isinstance(list_of_ios, list):
            self.InputOutput += list_of_ios
        else:
            self.InputOutput += [list_of_ios]
        for io in self.InputOutput:
            io.arrays = self.dataobjects_to_datasets_on_block(io.arrays)
        return

    def add_metric(self, metric_params):
        """ TODO V2 check if we are using
        """
        self.metric_transformations = metric_params
        return

    def get_all_scheme_halos(self):
        """This will be used later to change the halo sizes accordingly, Currently
        hard coded. If you want to manually change the halo points for arrays, change
        the value to the desired value in CentralHalos_defdec class of schemes file
        """
        spatialschemes = []
        for sc in self.discretisation_schemes:
            if self.discretisation_schemes[sc].schemetype == "Spatial":
                spatialschemes += [self.discretisation_schemes[sc]]
        halos = set([CentralHalos_defdec()])
        return halos

    @property
    def get_metric_class(self):
        """Returns the metric classes set on the block
        """
        metric_class = []
        for EqClass in self.list_of_equation_classes:
            if isinstance(EqClass, MetricsEquation):
                metric_class += [EqClass]
        if len(metric_class) == 0:
            return None
        if len(metric_class) > 1:
            raise ValueError("more than one metric class found in the equations")
        else:
            return metric_class[0]

    @property
    def fd_metrics(self):
        """Returns the DatSets where the metrics are evaluated to
        """
        def _convert(x):
            return self.dataobjects_to_datasets_on_block([x])[0]
        metric = self.get_metric_class
        if metric:
            return metric.FD_metrics.applyfunc(_convert)
        else:
            return eye(self.ndim)


def sort_constants(constants_dictionary):
    """ TODO V2 we are not using this
    """
    known_constants, unknown_constants = [], []
    pprint(constants_dictionary)
    for const in constants_dictionary.values():
        if const.is_input:
            known_constants.append(const)
        else:
            unknown_constants.append(const)
    while len(unknown_constants) != 0:
        set_of_known = set(known_constants)
        for const in unknown_constants:
            requires = const.value.atoms(ConstantObject)
            if requires.issubset(set_of_known):
                print("const: ", const, " has formula: ", const.value, " requires: ", requires)
                known_constants.append(const)
                unknown_constants = [x for x in unknown_constants if not const]
            else:
                print(const, "is missing", " it requires", requires)
    return known_constants
