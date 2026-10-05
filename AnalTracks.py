import numpy as np
import matplotlib.pyplot as plt
import ROOT
import math as m

f = open('dati.txt', 'r')

for line in f:
    if len(line) == 1:
        evento = line
        print('evento')
    while len(line) == 5:
        np.append(hit, line)
        print('hit')
print(hit)
f.close()
