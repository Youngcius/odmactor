import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as C

fc = np.loadtxt('contrast-by-single.txt')
plt.plot(*fc)
plt.xlim(2.855*C.giga, 2.885*C.giga)
plt.ylabel('contrast')
plt.show()

counts = np.loadtxt('counts.txt')
freqs = fc[0][5:]
plt.plot(freqs, counts[:,0][5:], color='r', label='MW on')
plt.plot(freqs, counts[:,1][5:], color='k', label='MW off')
plt.ylabel('Counts')
plt.legend()
plt.show()
