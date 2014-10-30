'''
Created on 30.07.2014

@author: naxobIdeaPad
'''

import matplotlib.pyplot as plt
from copy import deepcopy
from numpy import mean,std
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
import math

if __name__ == '__main__':
    fig = plt.figure(1)
    plt.plot([3,2,None,4,3])
    plt.xlim(0,10)
    plt.ylim(0,10)
    plt.show()
