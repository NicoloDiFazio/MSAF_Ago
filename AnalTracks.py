import numpy as np
import matplotlib.pyplot as plt
import ROOT
import math as m

class hit:
    def __init__(self, x, y, z, trkind, layind):
        self.x = x
        self.y = y
        self.z = z
        self.trkind = trkind
        self.layind = layind
    def __call__(self):
        return [self.x, self.y, self.z, self.trkind, self.layind]
    def __str__(self):
        return f"{self.x:.6f} {self.y:.6f} {self.z:.6f} {self.trkind} {self.layind}"
    def tocam(self):
        return cam(self.x, self.y, self.z)
  
    
class event:
    def __init__(self, eveind: int, data: list = None, angle: list = None):
        self.eveind = eveind
        self.data = data or []
        self.angle = angle or []
    def __call__(self):
        return [self.eveind, self.data, self.angle]
    def append(self, value):
        self.data.append(value)
        r, t = value.tocam()
        self.angle.append(t)
    def sort(self):
        return sorted(self.angle)
    def sync_sort(self):
        couple = sorted(list(zip(self.data, self.angle)),key=lambda a: a[1])
        dat, ang = map(list, zip(*couple))
        return event(self.eveind, dat, ang)
    def data_print(self):
        for dat in self.data: print(str(dat))
    def ind_print(self):
        print(str(self.eveind))
    def angle_print(self):
        for ang in self.angle: print(str(ang))
    def fin(self):
        ordine = self.sort()
        differenze = np.diff(ordine)
        #return float(np.median(differenze)*0.275)
        #return float(np.quantile(differenze, 0.3))
        return sigma_f
    def len(self):
        return len(self.angle)

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
    theta = m.atan2(y, x)
    phi = m.acos(z/rho)
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
l_evento = 0
evs = []
ev = event(l_evento)     #eventi con rho, theta lungo piano trasverso

for line in f:
    lettura = line.split()
    if len(lettura) == 1:
        if (int(lettura[0]) != l_evento):
            l_evento = int(lettura[0])
            evs.append(ev)
            ev = event(l_evento)
    elif len(lettura) == 5:
        the_hit = hit(float(lettura[0]), float(lettura[1]), float(lettura[2]),
                      int(lettura[3]), int(lettura[4]))
        ev.append(the_hit)

evs.append(ev)
f.close()

plt.figure(figsize=(8, 5))
#fig, axes = plt.subplots(1, 3, figsize=(12, 4))


colore = ['r', 'g', 'b', 'k', 'm']
ev = evs[0]
x = np.array([dato.x for dato in ev.data])
cls = np.array([colore[dato.layind] for dato in ev.data])

plt.scatter(ev.angle, x, color=cls)
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
