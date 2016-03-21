class Scheme(object):

    """ A numerical discretisation scheme. """

    def __init__(self, name, order):
        """ Initialise the scheme.
        
        :arg str name: The name of the scheme.
        :arg int order: The order of the scheme.
        """
        
        self.name = name
        self.order = order
        return
