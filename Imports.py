'''Ce module nous permet de regrouper tous les imports
de bibliothèques Python dans un seul module et d'importer
le module dans chaque fonctionnalité.'''

import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
from pytest import *
from scipy.integrate import solve_ivp
from functools import partial
import PIL
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
import math as mp
from mayavi import mlab
from tvtk.api import tvtk
from scipy.optimize import minimize
from functools import partial