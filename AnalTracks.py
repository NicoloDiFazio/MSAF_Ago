import numpy as np
import matplotlib.pyplot as plt
import ROOT
import math as m

class hit: #aggiunta per comodita' la colonna indice d'evento
    def __init__(self, x, y, z, trkind, layind, eveind):
        self.x = x
        self.y = y
        self.z = z
        self.trkind = trkind
        self.layind = layind
        self.eveind = eveind
    def __call__(self):
        return [self.x, self.y, self.z, self.trkind, self.layind, self.eveind]
    def tostr(self):
        return f"{self.x:.6f} {self.y:.6f} {self.z:.6f} {self.trkind} {self.layind} {self.eveind}"
    def pos(self):
        return self.x, self.y, self.z
    def tocil(self):
        return cil(self.x, self.y, self.z)
    def tosfe(self):
        return sfe(self.x, self.y, self.z)

def cil(x, y, z):
    rho = m.sqrt(x**2 + y**2)
    theta = m.atan2(y, x)
    return rho, theta, z
def sfe(x, y, z):
    rho = m.sqrt(x**2 + y**2 + z**2)
    theta = m.atan2(y, x)
    phi = m.acos(z/rho)
    return rho, theta, phi
    

f = open('dati.txt', 'r')
hits = []
for line in f:
    dati = line.split()
    if len(dati) == 1:
        evento = int(dati[0])
    elif len(dati) == 5:
        the_hit = hit(float(dati[0]), float(dati[1]), float(dati[2]), int(dati[3]), int(dati[4]), evento)
        hits.append(the_hit)
f.close()

plt.figure(figsize=(8, 5))

rhs = np.array([])
ths = np.array([])
phs = np.array([])
cls = np.array([])
colore = ['r', 'g', 'b', 'k', 'm', 'c'] 
for colpo in hits:
    if (colpo.eveind != 0): break
    #if (colpo.eveind != 0 or colpo.trkind !=0): break
    rh, th, ph = colpo.tosfe()
    rhs = np.append(rhs, rh)
    ths = np.append(ths, th)
    phs = np.append(phs, ph)
    cls = np.append(cls, colore[colpo.layind])
    #print(colpo())
"""
plt.scatter(phs, rhs, color=cls)    
plt.title("Grafico di rho in funzione di phi")
plt.xlabel("phi")
plt.ylabel("rho")
#plt.legend()
plt.show()
"""
#plt.hist(dati, bins=numerobin, range=(x_min, x_max), density=True, color='skyblue', edgecolor='black')
nbin = 5*int(len(phs)/6)  #ndati/nstrati
print(nbin)
bin_edges = []
plt.hist(phs, bins=nbin)

plt.title("Conteggi di tracce con angolo phi")
plt.xlabel("phi")
#plt.ylim(0, 7)
plt.ylabel("Conteggi")
plt.grid(True, alpha=1)
plt.show()
#"""
