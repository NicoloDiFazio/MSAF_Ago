import numpy as np
import matplotlib.pyplot as plt
from ROOT import TCanvas, TH1D, TF1, gApplication, TGraphErrors
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

pts_reali = [0.5, 1, 1.5, 2, 5, 10]

#errori sulle x, y e z
def sig(coord):
    return coord/m.sqrt(12)

risposta = input("output verboso? (y/n o t/f) ").strip().lower()
verbose = risposta in ["y", "t"]

risposta = input("output grafico? (y/n o t/f) ").strip().lower()
graphics = risposta in ["y", "t"]

f = open('test.txt', 'r')
l_evtest = 0
evstest = []
evento_test = event(l_evtest)

for line in f:
    lec = line.split()
    if len(lec) == 1:
        if (int(lec[0]) != l_evtest):
            l_evtest = int(lec[0])
            evstest.append(evento_test)
            evento_test = event(l_evtest)
    elif len(lec) == 5:
        h = hit(float(lec[0]), float(lec[1]), float(lec[2]),
                int(lec[3]), int(lec[4]))
        evento_test.append(h)
evstest.append(evento_test)
f.close()

sigmas_f = []
for evtest in evstest:
    ass = np.array([a for a in evtest.angle])
    sigmas_f.append(np.std(ass))

for s in sigmas_f: print(s)
#################################################
#              ANALISI EFFETTIVA                #
#################################################

f = open('dati2.txt', 'r')
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

if graphics: 
    c1 = TCanvas("c1", "Distribuzione Angolare dei Dati divisi per evento", 800, 600)
    col = 5
    row = m.ceil(len(evs)/col)
    c1.Divide(col, row)

lista_histo = []
for i, ev in enumerate(evs):
    sigma_f = sigmas_f[-1]
    if graphics: c1.cd(i+1)
    thm, thM = min(ev.angle), max(ev.angle)
    nb = 100000
    h = TH1D(f"h_{i}", f"Evento {i}; Theta; Conteggi", nb, thm, thM)
    for ang in ev.angle: h.Fill(ang)

    if graphics: 
        h.SetFillColor(38) # Azzurro
        h.SetLineColor(4) # Blu
        h.SetLineWidth(2)
        h.Draw("HIST")
        
    lista_histo.append(h)
    if verbose: print(f"Larghezza bin = {h.GetXaxis().GetBinWidth(0)/(sigma_f)}sigma")
if graphics: 
    c1.Modified()
    c1.Update()

#################################################
#              RICERCA VICINANZA                #
#################################################

listaTracceEventi = []
tr_tro = []
tr_tot = []
somma_tot = 0
lista_canvas = []

risposta = input("vuoi visualizzare i grafici specifici? (y/n o t/f) ").strip().lower()
allgraph = risposta in ["y", "t"]

for i, hist in enumerate(lista_histo):
    sigma_f = sigmas_f[-1]
    lar = sigma_f*3
    if allgraph:
        c2 = TCanvas(f"c2_{i}", f"Distribuzione Angolare dei Dati dell' evento {i}", 800, 600)
        hist.Draw("HIST")
        c2.Modified()
        c2.Update()
        lista_canvas.append(c2)
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
        if j != j_0 and ((c + traccia.len()) > 5 or (angolo - a_0) > lar):
            if (traccia.len()>3): listaTracce.append(traccia)
            indice += 1
            traccia = event(indice)
            traccia.append(ev.data[angint])
            a_0 = angolo
        else:
            traccia.append(ev.data[angint])
        j_0 = j
    if traccia.len() > 3: listaTracce.append(traccia)
    listaTracceEventi.append(listaTracce)
    
    tr_tro.append(len(listaTracce))
    tr_tot.append(ev.len()/5)
    somma_tot += ev.len()/5
    
    som = 0
    
    if verbose:
        for tr in listaTracce:
            som += 1
            print(f"angolo numero {som}")
            tr.angle_print()
    if True:
        print(f"evento numero {i}")
        print(f"tracce ricostruite {tr_tro[i]}/{tr_tot[i]}={(tr_tro[i]/tr_tot[i]):.2f}")

listaGrafici = []
Rss = []
eRss = []
pTs = []
epTs = []
somma = 0
count = 0
Ris = []
for pt_re in pts_reali:
    Ris.append(10/3*pt_re/B)


risposta = input("vuoi visualizzare i fit specifici? (y/n o t/f) ").strip().lower()
allgraph = risposta in ["y", "t"]

for eint, e in enumerate(listaTracceEventi):
    Rs = []
    eRs = []
    ps = []
    if allgraph: 
        c3 = TCanvas(f"c_{eint}", "Fit circolare sulle tracce", 800, 600)
        col = 10
        row = m.ceil(tr_tro[eint]/col)
        c3.Divide(col, row)
    for tint, t in enumerate(e):
        if allgraph: c3.cd(tint+1)
        indici_di_traccia = np.array([dato.trkind for dato in t.data])
        x = np.array([dato.x for dato in t.data])
        ex = np.array([sig(Xsizes[dato.layind]) for dato in t.data])
        y = np.array([dato.y for dato in t.data])
        ey = np.array([sig(Ysizes[dato.layind]) for dato in t.data])
        g = TGraphErrors(t.len(), x, y, ex, ey)
        g.SetTitle(f"Fit circolare del raggio di curvatura della traccia {tint}")
        listaGrafici.append(g)
        
        xm, xM = min(x), max(x)
        circ = "[0]+sqrt(pow([1], 2) - pow(x - [2], 2))"
        f_circ = TF1(f"f_circ_{tint}", circ, xm, xM)
        ri   = Ris[eint]
        yci  = -ri
        xci  = 0.0
        
        f_circ.SetParameters(yci, ri, xci)
        g.Fit(f_circ, "RQS")
        if allgraph: g.Draw("AP")
        ycf,  Rf,  xcf  = (f_circ.GetParameter(i) for i in range(3))
        eycf, eRf, exff = (f_circ.GetParError(i)  for i in range(3))
        chi2, ndf, p_value = f_circ.GetChisquare(), f_circ.GetNDF(), f_circ.GetProb()

        pT = 0.3*Rf*B
        epT = pT*eRf/Rf
        
        if p_value >= 0.05 and len(np.unique(indici_di_traccia)) == 1:
            Rs.append(Rf)
            eRs.append(eRf)
            pTs.append(pT)
            epTs.append(epT)
        
        if verbose:
            print(f"Chi2 / NDF = {chi2:.2f} / {ndf} (p-value: {p_value:.4f})")
    if True:
        print(f"evento numero {eint}")
        print(f"tracce trovate {len(Rs)}/{tr_tot[i]}={(len(Rs)/tr_tot[i]):.2f}")

    somma += len(e)
    Rss.append(Rs)
    eRss.append(eRs)
    if graphics: 
        c3.Modified()
        c3.Update()
        lista_canvas.append(c3)

graphics=True

if graphics: 
    c4 = TCanvas(f"c4", "Confronto impulsi reali/calcolati", 800, 600)

pTm, pTM = min(pTs), max(pTs)
nb = len(pts_reali)*8

h_calc = TH1D(f"h_calc", f"Impulsi calcolati; pT; Conteggi", nb, pTm-1, pTM+1)
for pT in pTs: h_calc.Fill(pT)
"""
#fit = TF1(f"fit", "[0]*TMath::Log(x)/TMath::Log([1]) + [2]", pTm, pTM)
fit = TF1(f"fit", "[0]+[1]*atan([2]*x)", 0, pTM)
fit.SetParameters(0, 1, 10)
h_calc.Fit(fit, "R")
"""
if graphics:
    h_calc.SetFillColor(38) # Arancione
    h_calc.SetLineColor(4) # Rosso
    h_calc.SetLineWidth(2)
    h_calc.Draw("HIST")
    #fit.Draw("SAME")

if graphics: 
    c4.Modified()
    c4.Update()

#"""    
if graphics or allgraph: gApplication.Run(True)
