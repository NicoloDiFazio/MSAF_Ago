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

#errori sulle x, y e z
def sig(coord):
    return coord/m.sqrt(12)

risposta = input("output verboso? (y/n o t/f) ").strip().lower()
verbose = risposta in ["y", "t"]

risposta = input("output grafico? (y/n o t/f) ").strip().lower()
graphics = risposta in ["y", "t"]

f = open('test.txt', 'r')
evento_test = event(0)
for line in f:
    lec = line.split()
    h = hit(float(lec[0]), float(lec[1]), float(lec[2]),
            int(lec[3]), int(lec[4]))
    evento_test.append(h)
f.close()

if graphics: c0 = TCanvas("c0", "Fit Distribuzione Angolare", 800, 600)

# Costruttore: TH1D(nome, titolo_con_assi, nbins, xlow, xup)
thm, thM = min(evento_test.angle), max(evento_test.angle)
h_test = TH1D("h_test", "Distribuzione degli angoli a theta fissato; Theta; Conteggi", 20, thm, thM)
for ang in evento_test.angle: h_test.Fill(ang)

h_test.SetFillColor(38) # Azzurro
h_test.SetLineColor(4) # Blu
h_test.SetLineWidth(2)

fit_test = TF1("fit_test", "gaus", thm, thM)

A_i     = h_test.GetMaximum()
mu_i    = 0.5
sigma_i = h_test.GetRMS()

fit_test.SetParameters(A_i, mu_i, sigma_i)
h_test.Fit("fit_test", "R")

A_f,  mu_f,  sigma_f  = (fit_test.GetParameter(i) for i in range(3))
eA_f, emu_f, esigma_f = (fit_test.GetParError(i)  for i in range(3))

chi2, ndf, p_value = fit_test.GetChisquare(), fit_test.GetNDF(), fit_test.GetProb()

if verbose:
    print(f"Chi2 / NDF = {chi2:.2f} / {ndf} (p-value: {p_value:.4f})")
    print(f"mu = {mu_f:.2e} +/- {emu_f:.2e} rad")
    print(f"sigma = {sigma_f:.2e} +/- {esigma_f:.2e} rad")

if graphics: 
    h_test.Draw("HISTO")
    fit_test.Draw("SAME")
    c0.Update()

ass = np.array([a for a in evento_test.angle])
"""
print(f"1. {sigma_f}     2. {np.std(ass)}")
scelta = int(input("quale sigma si vuole usare ? (1/2) "))
if scelta == 2:
    sigma_f = np.std(ass)
    
print(sigma_f)
"""
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
if graphics: 
    c1 = TCanvas("c1", "Distribuzione Angolare dei Dati divisi per evento", 800, 600)
    col = 5
    row = 2
    c1.Divide(col, row)

lista_histo = []
for i, ev in enumerate(evs):
    if graphics: c1.cd(i+1)
    thm, thM = min(ev.angle), max(ev.angle)
    #nb = int((thM - thm)/(sigma_f*4))
    nb = 1000000
    #print(nb)
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

lar = sigma_f*16
listaTracceEventi = []
tr_tro = []
tr_tot = []
somma_tot = 0
lista_canvas = []

risposta = input("vuoi visualizzare i grafici specifici? (y/n o t/f) ").strip().lower()
allgraph = risposta in ["y", "t"]


for i, hist in enumerate(lista_histo):
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
    
    if verbose:
        print(f"evento numero {i}")
        print(f"tracce trovate {tr_tro[i]}/{tr_tot[i]}={(tr_tro[i]/tr_tot[i]):.2f}")

listaGrafici = []
Rs = []
eRs = []
ps = []
somma = 0
count = 0
for eint, e in enumerate(listaTracceEventi):
    if graphics: 
        c3 = TCanvas(f"c_{eint}", "Fit circolare sulle tracce", 800, 600)
        col = 10
        row = m.ceil(tr_tro[eint]/col)
        c3.Divide(col, row)
    for tint, t in enumerate(e):
        if graphics: c3.cd(tint+1)
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
        yci  = -33.0
        ri   = 33.0
        xci  = 0.0
        
        f_circ.SetParameters(yci, ri, xci)
        #f_circ.SetParLimits(1, 30.0, 40.0)
        g.Fit(f_circ, "RQS")
        if graphics: g.Draw("AP")
        ycf,  Rf,  xcf  = (f_circ.GetParameter(i) for i in range(3))
        eycf, eRf, exff = (f_circ.GetParError(i)  for i in range(3))
        chi2, ndf, p_value = f_circ.GetChisquare(), f_circ.GetNDF(), f_circ.GetProb()
        
        if p_value >= 0.05:
            Rs.append(Rf)
            eRs.append(eRf)
            ps.append(p_value)
            if len(np.unique(indici_di_traccia)) == 1:
                count += 1
        
        if verbose:
            print(f"Chi2 / NDF = {chi2:.2f} / {ndf} (p-value: {p_value:.4f})")
    somma += len(e)

    if graphics: 
        c3.Modified()
        c3.Update()
        lista_canvas.append(c3)

pTs = np.array([])
epTs = np.array([])
scorre = np.array([])

for rind, r in enumerate(Rs):
    pT = 0.3*B*r
    epT = pT*eRs[rind]/r
    pTs = np.append(pTs, pT)
    epTs = np.append(epTs, epT)
    scorre = np.append(scorre, float(rind+1))
    if verbose:
        print(f"Raggio R = {r} +- {eRs[rind]} m")
        print(f"Impulso pT = {pT} +- {epT} GeV/c")
        print(f"p_value = {ps[rind]}")

esc = np.zeros(len(scorre))
c4 = TCanvas("c4", "pTs", 800, 600)
gr = TGraphErrors(len(pTs), scorre, pTs, esc, epTs)
gr.SetTitle("Grafico pTs")
gr.SetMarkerSize(10)
gr.Draw("AP")
c4.Modified()
c4.Update()
        
print(f"Ho trovato {len(Rs)} impulsi su {somma} su un totale di {somma_tot}")
print(f"Ho trovato {count} tracce della stessa traccia")
#"""    
if graphics or allgraph: gApplication.Run(True)
gApplication.Run(True)
