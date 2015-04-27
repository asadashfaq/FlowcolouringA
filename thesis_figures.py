from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from aurespf.tools import get_q

# Flow histogram with quantiles
F = np.load('./results/linear-flows.npy')
f, bins = np.histogram(F[0], bins=400, normed=True)

qs = [get_q(F[0], q) for q in [0.005, 0.01, 0.99, 0.995]]
topBot = [0, 0.00015]

textHeight = 0.000153

plt.figure(figsize=(11, 8))
plt.subplot(211)
plt.fill_between(bins[:-1], f, color='#0000aa', edgecolor='#0000aa')
#plt.plot([m, m], [0, 1], '-k', lw=2)
for q in qs:
    plt.plot([q, q], topBot, '-', lw=2, color='#aa0000')
plt.text(qs[0] - 1080, textHeight, '0.5%')
plt.text(qs[1] + 0, textHeight, '1%')
plt.text(qs[2] - 950, textHeight, '99%')
plt.text(qs[3] + 0, textHeight, '99.5%')
plt.text(bins[0] + 400, 0.9 * 0.0002, 'a)', fontsize=16)
plt.axis([bins[0], bins[-2], 0, 0.0002])
plt.xlabel(r'$F_l$ [MW]')
plt.ylabel(r'$p(F_l)$')


f, bins = np.histogram(abs(F[0]), bins=400, normed=True)
qs = [get_q(abs(F[0]), q) for q in [0.99, 0.995]]
topBot = [0, 0.00015]

plt.subplot(212)
plt.fill_between(bins[:-1], f, color='#0000aa', edgecolor='#0000aa')
#plt.plot([m, m], [0, 1], '-k', lw=2)
for q in qs:
    plt.plot([q, q], topBot, '-', lw=2, color='#aa0000')
plt.plot([max(abs(F[0])), max(abs(F[0]))], topBot, '-', lw=2, color='#aa0000')
plt.text(qs[0] - 550, textHeight, '99%')
plt.text(qs[1] + 0, textHeight, '99.5%')
plt.text(max(abs(F[0])) - 700, textHeight, '100%')
plt.text(bins[0] + 200, 0.9 * 0.0002, 'b)', color='#f6f6f6', fontsize=16)
plt.axis([bins[0], 14000, 0, 0.0002])
plt.xlabel(r'$F_l$ [MW]')
plt.ylabel(r'$p(F_l)$')

plt.savefig('./figures/thesis/flowHist.pdf', bbox_inches='tight')
