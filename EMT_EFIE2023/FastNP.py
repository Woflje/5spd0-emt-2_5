# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 17:05:18 2023

@author: 20192866
"""

import numpy as np
def fastCross(v1,v2):
    #reshape vectors into shape accepted by np.einsum
    v1=np.reshape(v1,(1,3))
    v2=np.reshape(v2,(1,3))
    #use einstein sum notation to perform cross product
    eijk=eijk = np.zeros((3, 3, 3))
    eijk[0, 1, 2] = eijk[1, 2, 0] = eijk[2, 0, 1] = 1
    eijk[0, 2, 1] = eijk[2, 1, 0] = eijk[1, 0, 2] = -1
    s=np.einsum('iuk,vk->uvi', np.einsum('ijk,uj->iuk', eijk, v1), v2)
    return np.reshape(s,(3)) #reshape result so the existing fucntions accept it

def fastNorm(v):
    return np.sqrt(np.sum(v**2, axis=-1))