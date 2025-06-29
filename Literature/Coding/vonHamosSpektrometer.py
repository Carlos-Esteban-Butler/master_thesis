from mmpxrt import mmpxrt

p=mmpxrt.init()
p['simulation']['numraysE']=3
p['simulation']['name']="DistancesCalculated_FullBandwidth_Test"

p['source']['EcentralRay']=13500
p['source']['EmaxBandwidth']= 13000

p['crystal']['d2']=6.708
p['crystal']['mosaicity']=0.4
p['crystal']['width']=10
p['crystal']['length']=150
p['crystal']['maxThickness']=0.1
p['crystal']['thickness']=-1
p['crystal']['radius_w']=22

p['geometry']['detRot']=-1
#p['geometry']['CrystalSource']=100
#p['geometry']['CrystalDetector']=100
p['geometry']['detectorLength']=300
p['geometry']['detectorWidth']=16
#p['geometry']['detectorPxSize']=200

rrrs = mmpxrt.spectrometer(p)
mmpxrt.spectrometer_evaluate(p,rrrs)
