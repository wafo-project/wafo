"""
    Modify this file if another plotbackend is wanted.
"""
import numpy as np
if False:
    try:
        from scitools import easyviz as plotbackend
        np.disp('wafo.wafodata: plotbackend is set to scitools.easyviz')
    except:
        np.disp('wafo: Unable to load scitools.easyviz as plotbackend')
        plotbackend = None
else:
    try:
        from matplotlib import pyplot as plotbackend
        np.disp('wafo.wafodata: plotbackend is set to matplotlib.pyplot')
    except:
        np.disp('wafo: Unable to load matplotlib.pyplot as plotbackend')
        plotbackend = None