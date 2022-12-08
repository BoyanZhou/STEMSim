

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from matplotlib import style
style.use('ggplot')
params = [0.5, 1, 2, 3]
x = np.linspace(0, 1, 100)
f, ax = plt.subplots(len(params), len(params), sharex=True, sharey=True)
for i in range(4):
  for j in range(4):
    alpha = params[i]
    beta = params[j]
    pdf = stats.beta(alpha, beta).pdf(x)
    ax[i, j].plot(x, pdf)
    ax[i, j].plot(0, 0, label='alpha={:3.2f}\nbeta={:3.2f}'.format(alpha, beta), alpha=0)
    plt.setp(ax[i, j], xticks=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0], yticks=[0, 2, 4, 6, 8, 10])
    ax[i, j].legend(fontsize=10)
ax[3, 0].set_xlabel('theta', fontsize=16)
ax[0, 0].set_ylabel('pdf(theta)', fontsize=16)
plt.suptitle('Beta PDF', fontsize=16)
plt.tight_layout()
plt.show()
