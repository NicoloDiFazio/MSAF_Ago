import numpy as np
import matplotlib.pyplot as plt
import ROOT
import math as m

"""
Obbiettivo: cercare tracce, ovvero i punti con stesso 
angolo polare lungo l'asse del rivelatore.

1)
Leggo i dati scritti come: 

len(line) = 1 -> evento
len(line) = 5 -> X Y Z IndiceTraccia IndiceLayer

2)
Converto i dati posizionali da X, Y, Z a
rho = sqrt(x^2 + y^2): distanza radiale dell'hit
      dall'asse del rivelatore z nel piano trasverso (x, y).
theta = atan(rho/z): angolo polare rispetto all'asse z.

3)
Per ogni evento genero un istogramma che conti quante
volte un determinato angolo e' stato raggiunto.
La larghezza di bin dev'essere variabile per evitare che zone 
troppo densamente popolate portino a conclusioni sbagliate

4)
Eseguo il fit della circonferenza per ottenere l'impulso
(la sua componente ortogonale a B o trasversa in gergo)
"""
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
    def __init__(self, eveind: int, data: list):
        self.eveind = eveind
        self.data = data
    def __call__(self):
        return [self.eveind, self.data]
    def append(self, value):
        self.data.append(value)
    def sort(self):
        return sorted(self.data)
    def data_print(self):
        for dat in self.data: print(str(dat))
    def ind_print(self):
        print(str(self.eveind))
    def fin(self):
        ordine = self.sort()
        differenze = np.diff(ordine)
        return float(np.median(differenze)*1/2)
    #prendiamo la mediana delle differenze degli angoli,
    #cioe' prendiamo la differenza tra angoli che non sia
    #ne troppo grande (tracce lontane) ne troppo piccole
    #(stessa traccia). La consideriamo sigma/2 e prendiamo
    # 3*sigma come finezza
    
def cam(x, y, z):
    rho = m.sqrt(x**2 + y**2)
    theta = m.atan2(rho, z)
    return rho, theta

def bin_edges(evento):
    ordine = evento.sort()
    bin_edge = [ordine[0]]
    for angolo in ordine:
        if (angolo - bin_edge[-1] >= evento.fin()):
            bin_edge.append(angolo)
    if (ordine[-1] > bin_edge[-1]):
        bin_edge.append(ordine[-1])
    return bin_edge

f = open('dati.txt', 'r')
l_evento = 0
#evs_raws = []
evs = []
#ev_raw = event(l_evento, []) #eventi con xyz indtrk indlay
ev = event(l_evento, [])     #eventi con rho, theta lungo piano trasverso

for line in f:
    lettura = line.split()
    if len(lettura) == 1:
        if (int(lettura[0]) != l_evento):
            l_evento = int(lettura[0])
            #evs_raws.append(ev_raw)
            evs.append(ev)
            #ev_raw = event(l_evento, [])
            ev = event(l_evento, [])
    elif len(lettura) == 5:
        the_hit = hit(float(lettura[0]), float(lettura[1]), float(lettura[2]),
                      int(lettura[3]), int(lettura[4]))
        #ev_raw.append(the_hit)
        rho, theta = the_hit.tocam()
        ev.append(theta)

#evs_raws.append(ev_raw)
evs.append(ev)
f.close()
'''
print(len(evs_raws))
for eve in evs_raws:
    eve.data_print()

print(len(evs))
for eve in evs:
    eve.ind_print()
'''
fig, axes = plt.subplots(2, m.ceil(len(evs)/2), figsize=(12, 4))

for i, ax in enumerate(axes.flatten()):
    if i >= len(evs): break
    ev = evs[i]
    ax.hist(ev.data, bins=bin_edges(ev), edgecolor='blue')
    ax.set_title(f"Conteggi $\\theta$ dell'evento {ev.eveind}")
    ax.set_xlabel("theta")
    ax.set_ylabel("Conteggi")
    ax.grid(True, alpha=1)
    
plt.tight_layout()
plt.show()
#'''
