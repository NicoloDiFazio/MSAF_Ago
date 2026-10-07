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
    def tocam(self):
        return cam(self.x, self.y, self.z)

def cam(x, y, z):
    rho = m.sqrt(x**2 + y**2)
    theta = m.atan2(rho, z)
    return rho, theta
    
def cil(x, y, z):
    rho = m.sqrt(x**2 + y**2)
    theta = m.atan2(y, x)
    return rho, theta, z

def sfe(x, y, z):
    rho = m.sqrt(x**2 + y**2 + z**2)
    theta = m.acos(z/rho)
    phi = m.atan2(y, x)
    return rho, theta, phi

def bin_edges(dati, finezza):
    dati_ordinati = sorted(dati)
    #angolo_0 = dati_ordinati[0]
    bin_edge = [dati_ordinati[0]]
    for angolo in dati_ordinati:
        if (angolo - bin_edge[-1] >= finezza):
            bin_edge.append(angolo)
            #angolo_0 = angolo
    if (dati_ordinati[-1] > bin_edge[-1]):
        bin_edge.append(angolo)
    return bin_edge

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

rhs = []
ths = []
phs = []
cls = []
colore = ['r', 'g', 'b', 'k', 'm', 'c'] 
for colpo in hits:
    if (colpo.eveind != 0): break
    #if (colpo.eveind != 0 or colpo.trkind !=0): break
    #if (colpo.eveind == 0 and colpo.trkind ==0):
    #if (colpo.eveind == 0):
    rh, th = colpo.tocam()
    rhs.append(rh)
    ths.append(th)
    #phs.append(ph)
    cls.append(colore[colpo.layind])
    #print(colpo())

fin = 5e-4
print(f"Finezza: {fin}")

#"""
plt.errorbar(ths, rhs, xerr=2*fin, linestyle='', marker='')
plt.scatter(ths, rhs, color=cls)
plt.title("Grafico di rho in funzione di theta")
plt.xlabel("theta")
plt.ylabel("rho")
#plt.legend()
plt.show()
"""

#plt.hist(dati, bins=numerobin, range=(x_min, x_max), density=True, color='skyblue', edgecolor='black')

plt.hist(ths, bins=bin_edges(ths, fin), edgecolor='blue')
plt.title("Conteggi di tracce con angolo theta")
plt.xlabel("theta")
plt.ylabel("Conteggi")
#plt.ylim(0, 7)
plt.grid(True, alpha=1)
plt.show()
#"""
