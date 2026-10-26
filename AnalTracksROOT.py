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

c = TCanvas("c", "Fit Distribuzione Angolare", 800, 600)

# Costruttore: TH1D(nome, titolo_con_assi, nbins, xlow, xup)
thm, thM = min(evento_test.angle), max(evento_test.angle)
h_test = TH1D("h_test", "Distribuzione degli angoli a theta fissato; Theta; Conteggi", 15, thm, thM)
for ang in evento_test.angle: h_test.Fill(ang)

# Personalizzazione estetica
h_test.SetFillColor(38) # Azzurro tenue
h_test.SetLineColor(4) # Blu
h_test.SetLineWidth(2)

fit_test = TF1("fit_test", "[0] * exp(-0.5 *pow((x - [1]) / [2], 2))", thm, thM)

A_i     = 1.1e3
mu_i    = 0.5
sigma_i = 2e-6

fit_test.SetParameters(A_i, mu_i, sigma_i)
h_test.Fit("fit_test")

A_f = fit_test.GetParameter(0)
mu_f = fit_test.GetParameter(1)
sigma_f = fit_test.GetParameter(2)

eA_f = fit_test.GetParError(0)
emu_f = fit_test.GetParError(1)
esigma_f = fit_test.GetParError(2)

chi2 = fit_test.GetChisquare()
ndf = fit_test.GetNDF()
p_value = fit_test.GetProb()

if verbose:
    print(f"Chi2 / NDF = {chi2:.2f} / {ndf} (p-value: {p_value:.4f})")
    print(f"$$mu$$ = {mu_f:.3f} +/- {emu_f:.3f} keV")
    print(f"$$sigma$$ = {sigma_f:.3f} +/- {esigma_f:.3f} keV")

h_test.Draw("E") # Disegna l'istogramma con le barre d'errore
c.Update() # Mantiene aperta la finestra grafica ROOT
gApplication.Run(True)
"""
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

tracce_eventi = []
l = 4*sigma_f

for evento in evs:
    tracce = []
    indice = 0
    traccia = event(indice)
    ev = evento.sync_sort()
    a_prec = ev.angle[0]
    for i, angolo in enumerate(ev.angle):
        if angolo - a_prec <= l and traccia.len() < 5:
            traccia.append(ev.data[i])
            if verbose: print("messo hit in traccia")
        elif angolo - a_prec > l or traccia.len() == 5:
            if traccia.len() >= 3:
                tracce.append(traccia)
                if verbose: print("messa traccia in lista tracce")
            elif traccia.len() < 3:                
                if verbose: print("traccia skippata")
            indice += 1
            traccia = event(indice)
            if i+1 < ev.len():
                a_prec = ev.angle[i+1]
    tracce_eventi.append(tracce)
    if verbose: print("messe liste tracce in liste eventi")

ptss = []
sptss = []
for ind, traccias in enumerate(tracce_eventi):
    pts = []
    spts = []
    fig, axes = plt.subplots(10, m.ceil(len(traccias)/10), figsize=(12, 4))
    plt.suptitle(f"Tracce dell'evento {ind}")
    if verbose: print(f"Evento {ind}")
    for i, ax in enumerate(axes.flatten()):
        if i >= len(traccias): break
        tr = traccias[i]
        x = np.array([dato.x for dato in tr.data])
        sx = np.array([sig(Xsizes[dato.layind]) for dato in tr.data])
        y = np.array([dato.y for dato in tr.data])
        sy = np.array([sig(Ysizes[dato.layind]) for dato in tr.data])
        ax.errorbar(x, y, xerr=sx, yerr=sy, fmt = '.', label='Hit traccia')

        xci = np.mean(x)
        yci = np.mean(y)
        Ri = 33
        mi = Minuit(fcn, xc = xci, yc = yci, R = Ri)
        mi.errordef = Minuit.LEAST_SQUARES
        mi.migrad()

        theta = np.linspace(0, 2*np.pi, 200)
        xcf = mi.values['xc']
        ycf = mi.values['yc']
        Rf = mi.values['R']
        sRf = mi.errors['R']
        ax.plot(xcf + Rf*np.cos(theta), ycf + Rf*np.sin(theta), 'r--', label='Fit cerchio')
        ax.set_aspect('equal', adjustable='datalim')

        if verbose:
            if mi.valid:
                print("Il fit e' converso correttamente.")
            else:
                print("Attenzione: il fit non e' valido o non e' converso!")
            chi2 = mi.fval
            ndof = 2*len(x)-3  #2 coordinate per ogni hit meno 3 parametri liberi (xc, yc, R)
            print(f"Traccia {i}: Chi2/ndof = {chi2/ndof}")

        pT = 0.3*B*Rf
        spT = pT*sRf/Rf
        if verbose:
            print(f"Raggio R: {Rf} +- {sRf} m")
            print(f"Impulso pT: {pT} +- {spT} GeV/c")
        pts.append(pT)
        spts.append(spT)
    ptss.append(pts)
    sptss.append(spts)
    plt.tight_layout()
    plt.show()
    plt.close(fig)

for i, p in enumerate(ptss):
    sp = sptss[i]
    x = np.array([j for j, pt in enumerate(p)])
    plt.figure(figsize=(8, 5))
    plt.errorbar(x, p, sp, fmt='.', linewidth=2)
    plt.show()
#"""
