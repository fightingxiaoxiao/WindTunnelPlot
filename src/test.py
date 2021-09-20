import numpy as np
import matplotlib.pyplot as plt
a = np.random.normal(1, 1, 1201)
b = np.linspace(0, 60, 1201)
#t = np.arange(0.0, 2.0, 0.01)
#s = 1 + np.sin(2 * np.pi * t)

fig, ax = plt.subplots()
ax.plot(b, a)

ax.set(xlabel='time (s)', ylabel='voltage (mV)',
       title='About as simple as it gets, folks')
ax.grid()

plt.show()
