import matplotlib
matplotlib.use('Agg')
import pylab
params = {'text.fontsize': 18,
                 'axes.labelsize' : 18,
                          'legend.fontsize': 16,
                                   'xtick.labelsize': 18,
                                            'ytick.labelsize': 18,
                                                     'text.usetex': True,
                                                              'lines.markersize': 8.0}
pylab.rcParams.update(params)
matplotlib.rc('font',**{'family':'serif','serif':['Times']})
matplotlib.rc('lines',**{'linewidth': 1.5})
matplotlib.rc('text', usetex=True)

