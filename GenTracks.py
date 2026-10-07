import numpy as np
import matplotlib.pyplot as plt
import ROOT
import math as m

B = 1.0

um = 1e-6
xsize = 300 * um
ysize = 50 * um
zsize = 100 * um

###### Classi per i dati ######

class hit:
    def __init__(self, x, y, z, trkind, layind):
        self.x = x
        self.y = y
        self.z = z
        self.trkind = trkind
        self.layind = layind
    def __call__(self):
        return [self.x, self.y, self.z, self.trkind, self.layind]
    def tostr(self):
        return f"{self.x:.6f} {self.y:.6f} {self.z:.6f} {self.trkind} {self.layind}"

class track:
    def __init__(self, pt, phi, theta, charge):
        self.pt = pt
        self.phi = phi
        self.theta = theta
        self.charge = charge
        
###### Funzione per il calcolo dell'intersezione elica-piano ######

def intersec(pos, trk, charge, xl):
    #Intersection at the center of pixel
    xl_center = xl+xsize/2.0

    #Track parameters
    h = charge
    R = 0.3*B*trk.Pt()
    phi0 = trk.Phi()+charge*ROOT.TMath.Pi()/2.0

    xc = pos.X()-R*m.cos(phi0)
    yc = pos.Y()-R*m.sin(phi0)
    
    ang  = m.acos((xl-pos.X())/R+m.cos(phi0))
    if phi0<0:
        ang  = -ang

    s    = R*(ang-phi0)/(-h*m.sin(trk.Theta()))

    x    = pos.X()+R*(m.cos(phi0-h*s*m.sin(trk.Theta())/R)-m.cos(phi0))
    y    = pos.Y()+R*(m.sin(phi0-h*s*m.sin(trk.Theta())/R)-m.sin(phi0))
    z    = pos.Z()+s*m.cos(trk.Theta())

    pos.SetX(x)
    pos.SetY(y)
    pos.SetZ(z)

    phitn = m.atan2(yc-y,xc-x)+h*ROOT.TMath.Pi()/2

    trk.SetPhi(phitn)

    #normalizzo x,y,z
    y = m.floor(y/ysize)*ysize+ysize/2
    z = m.floor(z/zsize)*zsize+zsize/2

    return x,y,z

###### Generazione e disegno degli "hit" ######

#Geometria
Xlayer = [0.05, 0.1, 0.15, 0.25, 0.35, 0.45]
lim = 0.7

mPi = 0.1396  # Massa del pione carico (GeV)
rnd = ROOT.TRandom3()
rnd.SetSeed(1234567)

f = open('dati.txt', 'w')
#MC per eventi
for jevt in range(0,10): #n eventi
    f.write(str(jevt) + '\n')
    ntrack  = int(rnd.Gaus(100,20))
    for j in range(0,ntrack): #n tracce

        trk   = ROOT.TLorentzVector()
        costh = rnd.Rndm()*0.4-0.2
        zref  = 0 #rnd.Rndm()*0.1-0.05
        
        pt    = 10
        phi   = rnd.Rndm()*0.2-0.1
        theta = m.acos(costh)
        
        v = ROOT.TVector3()
        v.SetMagThetaPhi(pt/m.sin(theta),theta,phi)
        trk.SetVectM(v,mPi)
        
        charge = 1
        if rnd.Rndm()<0.5:
            charge = -1
            
        pos = ROOT.TVector3(0,0,zref)
        for i in range(0,len(Xlayer)):
            hits = intersec(pos, trk, charge, Xlayer[i])
            the_hit = hit(*hits, j, i)
            f.write(the_hit.tostr() + '\n')
            
f.close()
