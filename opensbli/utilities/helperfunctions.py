from sympy import *
from opensbli.core.opensblifunctions import *
from opensbli.core import *


def increment_dataset(expression, direction, value):
    """ Increments an expression containing datasets by the given increment and direction.
    arg: object: expression: A SymPy expression containing datasets to have their indices updated.
    arg: int: direction: The integer direction to apply the increment. (which DataSet axis to apply to)
    arg: int: value: The positive or negative change to apply to the DataSet's index.
    returns: object: expression: The original expression updated to the new DataSet location.
    """
    for dset in expression.atoms(DataSet):
        loc = list(dset.indices)
        loc[direction] = loc[direction] + value
        new_dset = dset.base[loc]
        expression = expression.replace(dset, new_dset)
    return expression


def dot(v1, v2):
    out = 0
    if isinstance(v1, list):
        if len(v1) == len(v2):
            for i in range(len(v1)):
                out += v1[i]*v2[i]
            return out
        else:
            raise ValueError("")
    else:
        return v1*v2


def decreasing_order(s1, s2):
    return cmp(len(s2.atoms(CoordinateObject)), len(s1.atoms(CoordinateObject)))


def increasing_order(s1, s2):
    return cmp(len(s1.atoms(CoordinateObject)), len(s2.atoms(CoordinateObject)))


def sort_funcitons(fns, increasing_order=True):
    """Sorts the functions based on the number of arguments in
    increasing order
    """
    if increasing_order:
        return (sorted(fns, cmp=increasing_order))
    else:
        return (sorted(fns, cmp=decreasing_order))


def get_inverse_deltas(delta):
    from opensbli.core.codegeneration.opsc import rc
    if delta in rc.existing:
        return rc.existing[delta]
    else:
        name = rc.name
        b, exp = delta.as_base_exp()
        rc.name = "inv_%d"
        inv_delta_name = rc.get_next_rational_constant(delta)
        rc.name = name
        return inv_delta_name
