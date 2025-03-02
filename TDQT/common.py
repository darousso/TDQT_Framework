import numpy as np
import matplotlib.pyplot as plt
import sys
import scipy.linalg
import time
from scipy.special import expit
import os
from scipy.sparse.linalg import gmres

hartree_to_eV = 27.2114
hbar=4.135667696e-15 #eV.s
kB = 8.617e-5  # Boltzmann constant in eV/K
e=1.60218e-19
invfs=1e15
fs=1e-15