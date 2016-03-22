from sympy import *

class FileIO(object):

    """ Saves the arrays provided after every n iterations. These will eventually be dumped into an HDF5 file. """
    
    def __init__(self, arrays, niter = None):
        """ Setup the 'save' arrays to dump to an HDF5 file.
        
        :arg arrays: The arrays to save to a file.
        :arg int niter: The number of iterations that should pass before the arrays are saved to a file. If niter is None, the arrays are saved at the end of the simulation.
        :returns: None
        """
        
        self.save_after = []
        self.save_arrays = []
        self.save_after.append(Symbol('niter', integer=True))
        if isinstance(arrays, list):
            self.save_arrays += arrays
        else:
            self.save_arrays.append(arrays)
        return

