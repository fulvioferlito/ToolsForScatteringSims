import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.lines as mlines
import seaborn as sns
from ff_plot_options import *

################################ SELECT OPTIONS ################################

simnames = ['lcdm', 'q_0', 'q_1', 'q_2', 'q_3', 'p_0', 'p_1', 'p_2','p_3']
#simnames = ['lcdm', 'q_2','q_3']
part_type = "GAS+CDM"                   #possible options are GAS, CDM, GAS+CDM
z         = ["0.0","0.5","1.0"]     #list of redshift as a string in #.# format

################################# LOAD SPECTRA #################################

pk = []                     #list of powerspectra
n_sims = len(simnames)
for i in range(n_sims):
    for z_j in z:
        pk.append(np.loadtxt('Pk/Pk_'+simnames[i]+'_'+part_type+'_z='+z_j+'00.dat'))

##################################### PLOT #####################################
matplotlib.rcParams.update({'font.size': 13})
fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, sharex=True, figsize=(10,12))
axes = (ax1,ax2,ax3,ax4,ax5,ax6)
plt.subplots_adjust(wspace=0.22,  hspace=0.07)
#fig.tight_layout()
for ax in (ax1,ax2,ax3,ax4,ax5,ax6):
    ax.set_xlim([0.1,12])
    ax.set_ylim([0,2.75])
    ax.set_xscale('log')
    #ax.set_aspect(0.4)
for ax in (ax1,ax3,ax5): ax.set_ylabel(r'$P(k)/P(k)_{\Lambda \mathrm{CDM}}$')
for ax in (ax2,ax4,ax6): ax.set_ylabel(r'$P(k)/P(k)_{\alpha_{xb}=0}$')
for ax in (ax5,ax6): ax.set_xlabel('$k[h\mathrm{Mpc}^{-1}]$')

for j_z in range(3):
    for i in range(n_sims-1,-1,-1):
        axes[2*j_z].plot(pk[i*3+j_z][:,0],pk[i*3+j_z][:,1]/pk[0+j_z][:,1],
                 color = colors[i], linestyle=linestyles[i], linewidth=linewidths[i])
        if z[j_z][2]=='5':
            txt=r'$z='+z[j_z]+"$"
        else:
            txt=r'$z='+z[j_z][0]+"$"
        axes[2*j_z].text(3,2.4,txt)

for j_z in range(3):
    for i in range(n_sims-1,0,-1):
        axes[2*j_z+1].plot(pk[i*3+j_z][:,0],pk[i*3+j_z][:,1]/pk[den[i]*3+j_z][:,1],
                 color = colors[i], linestyle=linestyles[i], linewidth=linewidths[i-1])
        if z[j_z][2]=='5':
            txt=r'$z='+z[j_z]+"$"
        else:
            txt=r'$z='+z[j_z][0]+"$"
        axes[2*j_z+1].text(3,2.4,txt)




first_legend_lines=[]
first_label = ['$\Lambda \mathrm{CDM}$',r'$\alpha_{xb}=0$',r'$\alpha_{xb}=1$',r'$\alpha_{xb}=10$',r'$\alpha_{xb}=100$']
for i in range(5):
    first_legend_lines.append(mlines.Line2D([], [], color=pal[i], label=first_label[i]))

second_legend_lines=[]
second_label = ['$w_x=-0.9$','$w_x=-1.1$']
second_linestyle = ['-','--']
for i in range(2):
    second_legend_lines.append(mlines.Line2D([], [], color='grey',
                               label=second_label[i],
                               linestyle=second_linestyle[i]))

leg1 = ax1.legend(handles=first_legend_lines,loc="upper left")
#ax7 = plt.gca().add_artist(first_legend)
leg2 = ax1.legend(handles=second_legend_lines, loc="lower left")

ax1.add_artist(leg1)
ax1.add_artist(leg2)

#plt.legend(loc="upper left", ncol=2)



plt.savefig('Figures/Pk_ratio_'+part_type+'.pdf', bbox_inches='tight')
plt.show()
plt.close()
