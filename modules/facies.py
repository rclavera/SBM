# -*- coding: utf-8 -*-
import numpy as np
import modules.var as v

try:
    ff = open('facies.txt', 'r')
  
    for i in range(v.nNode):
        for jj in range(v.jti):
            sedi[:] = sediPC[i, jj, :]
            for fa in range():
                if np.greater_equal(sedi, minimfacies[fa]) and np.less_equal(sedi, maximfacies):
                    facies[i][jj] == fa
    ff.close()
except:
    print("The facies file does not exist")
    quit()

