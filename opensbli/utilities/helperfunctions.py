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
        new_dset =  dset.base[loc]
        expression = expression.replace(dset, new_dset)
    return expression