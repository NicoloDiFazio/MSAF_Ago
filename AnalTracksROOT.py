import numpy as np
import matplotlib.pyplot as plt
from ROOT import TCanvas, TH1D, TF1, gApplication
import math as m
from iminuit import Minuit
from iminuit.cost import LeastSquares

"""
Leggo i dati scritti come: 

len(line) = 1 -> evento
len(line) = 5 -> X Y Z IndiceTraccia IndiceLayer
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
        return sigma_f
    def len(self):
        return len(self.angle)
    
def cam(x, y, z):
    rho = m.sqrt(x**2 + y**2)
    theta = m.atan2(rho, z)
    return rho, theta

#Parametri di generazione
B = 1.0
um = 1e-6

Xsizes = [300 * um] * 5
Ysizes = [50 * um, 50 * um, 200 * um, 200 * um, 200 * um]
Zsizes = [100 * um] * 5

#errori sulle x, y e z
def sig(coord):
    return coord/m.sqrt(12)

#Minuit
def fcn(xc, yc, R):
    dx = x - xc
    dy = y - yc
    r = np.sqrt(dx**2 + dy**2)
    sigma = (dx**2 * sx**2 + dy**2 * sy**2)/(dx**2+dy**2)
    return np.sum(((r - R)**2)/sigma)
"""
risposta = input("output verboso? (y/n o t/f) ").strip().lower()
verbose = risposta in ["y", "t"]
"""
verbose = True

f = open('test.txt', 'r')
evento_test = event(0)
for line in f:
    lec = line.split()
    h = hit(float(lec[0]), float(lec[1]), float(lec[2]),
            int(lec[3]), int(lec[4]))
    evento_test.append(h)
f.close()

c0 = TCanvas("c0", "Fit Distribuzione Angolare", 800, 600)

# Costruttore: TH1D(nome, titolo_con_assi, nbins, xlow, xup)
thm, thM = min(evento_test.angle), max(evento_test.angle)
h_test = TH1D("h_test", "Distribuzione degli angoli a theta fissato; Theta; Conteggi", 50, thm, thM)
for ang in evento_test.angle: h_test.Fill(ang)

h_test.SetFillColor(38) # Azzurro
h_test.SetLineColor(4) # Blu
h_test.SetLineWidth(2)

fit_test = TF1("fit_test", "gaus", thm, thM)

A_i     = h_test.GetMaximum()
mu_i    = 0.5
sigma_i = h_test.GetRMS()

fit_test.SetParameters(A_i, mu_i, sigma_i)
h_test.Fit("fit_test")

A_f,  mu_f,  sigma_f  = (fit_test.GetParameter(i) for i in range(3))
eA_f, emu_f, esigma_f = (fit_test.GetParError(i)  for i in range(3))

chi2, ndf, p_value = fit_test.GetChisquare(), fit_test.GetNDF(), fit_test.GetProb()

if verbose:
    print(f"Chi2 / NDF = {chi2:.2f} / {ndf} (p-value: {p_value:.4f})")
    print(f"mu = {mu_f:.2e} +/- {emu_f:.2e} rad")
    print(f"sigma = {sigma_f:.2e} +/- {esigma_f:.2e} rad")

h_test.Draw("E") # Disegna l'istogramma con le barre d'errore
c0.Update() # Mantiene aperta la finestra grafica ROOT

#################################################
#              ANALISI EFFETTIVA                #
#################################################

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

c1 = TCanvas("c1", "Distribuzione Angolare dei Dati divisi per evento", 800, 600)
col = 5
row = 2
c1.Divide(col, row)

lista_histo = []
for i, ev in enumerate(evs):
    c1.cd(i+1)
    thm, thM = min(ev.angle), max(ev.angle)
    #nb = int((thM - thm)/(sigma_f*4))
    nb = 100000
    #print(nb)
    h = TH1D(f"h {i}", f"Evento {i}; Theta; Conteggi", nb, thm, thM)
    for ang in ev.angle: h.Fill(ang)

    h.SetFillColor(38) # Azzurro
    h.SetLineColor(4) # Blu
    h.SetLineWidth(2)
    
    h.Draw("HIST")
    lista_histo.append(h)
    #print(f"Larghezza bin = {h.GetXaxis().GetBinWidth(0)/(sigma_f)}sigma")
    
c1.Update()
c1.Modified()

c2 = TCanvas("c2", "Distribuzione Angolare dei Dati divisi per evento", 800, 600)
lista_histo[0].Draw("HIST")
c2.Update()

lar = 4*sigma_f
listaTracceEventi = []
for i, hist in enumerate(lista_histo):
    listaTracce = []
    ev = evs[i].sync_sort()
    a_0 = ev.angle[0]
    j_0 = hist.FindBin(a_0)
    c_0 = hist.GetBinContent(j_0)

    indice = 0
    traccia = event(indice)    
    
    for angint, angolo in enumerate(ev.angle):
        j = hist.FindBin(angolo)
        c = hist.GetBinContent(j)
        if j != j_0 and ((c + c_0) > 5 or (angolo - a_0) > lar):
            if (traccia.len()>3):listaTracce.append(traccia)
            indice += 1
            traccia = event(indice)
            traccia.append(ev.data[angint])
            a_0 = angolo
        else:
            traccia.append(ev.data[angint])
        c_0 = traccia.len()
        j_0 = j
    listaTracce.append(traccia)
    som = 0
    """
    for tr in listaTracce:
        som += 1
        print(f"angolo numero {som}/{len(listaTracce)}")
        tr.angle_print()
    """
    listaTracceEventi.append(listaTracce)
    print(f"evento numero {i}, tracce trovate {len(listaTracce)} su {ev.len()/5}, efficenza di {(len(listaTracce)/(ev.len()/5)):.2f}")
    
gApplication.Run(True)
