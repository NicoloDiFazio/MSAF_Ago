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
        return [self.x, self.y, self.z]

def pol(x, y, z):
    return rho, theta, phi

f = open('dati.txt', 'r')
hits = []
for line in f:
    dati = line.split()
    if len(dati) == 1:
        evento = int(dati[0])
    elif len(dati) == 5:
        the_hit = hit(float(dati[0]), float(dati[1]), float(dati[2]), int(dati[3]), int(dati[4]), evento)
        hits.append(the_hit())
f.close()
