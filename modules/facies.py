

import numpy as np


try:
  ff = open('facies.txt','r')
  
  for i in range(v.nNode):
    for jj in range(v.jti):
      sedi[:]=sediPC[i,jj,:]
      for fa in range()
        if np.greater_equal(sedi,minimfacies[fa]) and np.less_equal(sedi,maximfacies): facies[i][jj]==fa
           
except:
  print "El fichero no existe"

