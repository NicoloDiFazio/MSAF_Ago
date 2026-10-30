import numpy as np
import matplotlib.pyplot as plt
import ROOT
import math as m

#Geometria e parametri

Xlayer = [0.1, 0.2, 0.4, 0.6, 0.8]
pts = [0.5, 1, 1.5, 2, 5, 10] #inizio da 1 se no hit non hittano

mPi = 0.1396  # Massa del pione carico (GeV)
B = 1.0

um = 1e-6
Xsizes = [300 * um] * 5
Ysizes = [50 * um, 50 * um, 200 * um, 200 * um, 200 * um]
Zsizes = [100 * um] * 5

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

def intersec(pos, trk, charge, l_i):
    xl = Xlayer[l_i]
    xsize = Xsizes[l_i]
    ysize = Ysizes[l_i] 
    zsize = Zsizes[l_i]
    #Intersection at the center of pixel
    xl_center = xl+xsize/2.0

    #Track parameters
    h = charge
    R = 10/3 * trk.Pt()/(h*B)
    phi0 = trk.Phi()+charge*ROOT.TMath.Pi()/2.0

    xc = pos.X()-R*m.cos(phi0)
    yc = pos.Y()-R*m.sin(phi0)
    
    ang  = m.acos((xl_center-pos.X())/R+m.cos(phi0))
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
rnd = ROOT.TRandom3()
rnd.SetSeed(1234567)

ft = open('test.txt', 'w')
#for jevt, pt in enumerate(pts): #n eventi
#ft.write(str(jevt) + '\n')
ntracce = 1000
for j in range(0, ntracce):
    traccia = ROOT.TLorentzVector()
    
    phi   = (rnd.Rndm()*2-1)*0.1
    theta = 0.5
    pt = rnd.Rndm()*(max(pts)-min(pts)) + min(pts)

    vec = ROOT.TVector3()
    vec.SetMagThetaPhi(pt/m.sin(theta), theta, phi)
    traccia.SetVectM(vec, mPi)

    carica = 1

    posizione = ROOT.TVector3(0, 0, 0)
    for i in range(0, len(Xlayer)):
        colpis = intersec(posizione, traccia, carica, i)
        the_hit = hit(*colpis, j, i)
        ft.write(the_hit.tostr() + '\n')
ft.close()

f = open('dati2.txt', 'w')
#MC per eventi
for jevt, pt in enumerate(pts): #n eventi
    f.write(str(jevt) + '\n')
    #ntrack  = int(rnd.Gaus(100,20))
    ntrack = 100
    for j in range(0, ntrack): #n tracce

        trk   = ROOT.TLorentzVector()
        costh = rnd.Rndm()*2.0-1.0
        zref  = 0
        
        phi   = (rnd.Rndm()*2-1)*0.1
        theta = m.acos(costh)
        
        v = ROOT.TVector3()
        v.SetMagThetaPhi(pt/m.sin(theta),theta,phi)
        trk.SetVectM(v,mPi)
        
        charge = 1
            
        pos = ROOT.TVector3(0,0,zref)
        for i in range(0,len(Xlayer)):
            hits = intersec(pos, trk, charge, i)
            the_hit = hit(*hits, j, i)
            f.write(the_hit.tostr() + '\n')
            
f.close()
