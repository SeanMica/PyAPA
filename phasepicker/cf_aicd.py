import numpy as np
from .util import *
from .aicd import aicd
#import SeanAICD-1.0-py3.5.egg-info as sean

class AicDeriv():
  
  def __init__(self, trace):
    self.tr = trace  
    
  def _statistics(self):
    npts = self.tr.stats.npts
    data = self.tr.data;
    delta = 1.0/self.tr.stats.sampling_rate
    
    AIC, AIC_deriv = aicd(data,npts)
    #AIC = np.array(AIC)
    
    #AIC_deriv = []
    #for i in range(npts-1):
    #  b = np.abs(AIC[i+1]-AIC[i])
    #  AIC_deriv.append(b)
      
    #AIC_deriv.insert(0,AIC_deriv[0])
    #AIC_deriv = np.array(AIC_deriv)

    return AIC, AIC_deriv