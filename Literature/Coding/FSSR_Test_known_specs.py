# equations for a FSSR 1D spectrometer

from mmpxrt import mmpxrt
import numpy as np

R = 100
n = 2
twod = 19.915
Ecentral = 1400
theta = np.arcsin(n*(12398/Ecentral)/twod)
psi = 90 - theta
b = R*np.sin(theta)
a = np.absolute(R*np.cos(theta)/np.cos(2*(theta)))

p=mmpxrt.init()
p['simulation']['numraysE']=3
p['simulation']['num_processes']=1
p['simulation']['name']="FSSRFirstTestLowBandwidth"
p['simulation']['numerical_intersect_halving']=1

p['source']['EcentralRay']=Ecentral
p['source']['EmaxBandwidth']= 40
#p['source']['divergenceFWHM']= 100/180*np.pi

p['crystal']['d2']=twod
p['crystal']['mosaicity']=0
p['crystal']['width']=8
p['crystal']['length']=26
p['crystal']['maxThickness']=0.1
p['crystal']['thickness']=0
p['crystal']['radius_w']=R
p['crystal']['radius_l']=R

p['geometry']['detRot']=0
p['geometry']['CrystalSource']=a
p['geometry']['CrystalDetector']=b
p['geometry']['ThBragg0']=theta
rrrs = mmpxrt.spectrometer(p)
mmpxrt.spectrometer_evaluate(p,rrrs)