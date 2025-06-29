from mmpxrt import mmpxrt
from numpy import sin, pi

# simulation parameters
SimName = "ADP final"
SimLength = 6 # amount of rays, exponent of exponential
a = 174.317 # crystal to source
l = 40 # length of crystal
n = 1 # first order
twod = 10.64/n # lattice spacing of ADP in angstrom

thetaC = 47.5815 * pi/180 # give central angle from inventor file
phi = pi/2-thetaC
b =174.317 


phimin = (90-46)*pi/180
phimax = (90-49.14)*pi/180

Emax = 12398/(sin(pi/2-phimin)*twod)
Emin = 12398/(sin(pi/2-phimax)*twod)
print(Emin, Emax)
Bandwidth = Emax-Emin
Emath = Bandwidth/2+Emin
print(Emath)
print(Bandwidth)
intreflect = 0.00000232 # taken from Ferrari paper in ADP folder
rcadp = 0.0008 # approximated from paper from Rajesh in ADP Data folder. 
# in actuality is 0.000165, but that seems to be too small for the simulation to deal
# with. This setting is already sigficantly better than Mica. is adjusted until it 
# becomes the minimum that still lets simulation function. 

p=mmpxrt.init()
p['simulation']['numraysE']=SimLength
p['simulation']['name']=SimName
# p['simulation']['numerical_intersect_halving']=1 
p['simulation']['num_processes']=1

p['source']['EcentralRay']=Emath
p['source']['EmaxBandwidth']= Bandwidth # limit to exactly the range in real spectrometer
p['source']['size']= 100e-3 # set source size to cube with length of 100 micrometer

p['crystal']['d2']=twod
p['crystal']['width']=30
p['crystal']['length']=l
p['crystal']['maxThickness']=2
p['crystal']['thickness']=-1 # use exp func
p['crystal']['mosaicity']=0 # 0 bc monocrystal
# p['crystal']['rockingCurveFWHM']=rcadp # include bc monocrystal
# p['crystal']['integrated_reflectivity']=intreflect # taken from same paper as rc
p['crystal']['rockingCurveFWHM']=rcadp 
p['crystal']['integrated_reflectivity']=intreflect

p['geometry']['detRot']=-47.58
p['geometry']['CrystalSource']=a
p['geometry']['CrystalDetector']=b
p['geometry']['detectorLength']=27.6
p['geometry']['detectorWidth']=6.9
p['geometry']['detectorPxSize']=13.5e-3

rrrs = mmpxrt.spectrometer(p)
mmpxrt.spectrometer_evaluate(p,rrrs)