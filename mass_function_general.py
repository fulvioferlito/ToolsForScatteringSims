import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import readfof
import readsubf
import readgadget
import mass_function_library as MFL
import matplotlib.lines as mlines
import seaborn as sns
from ff_plot_options import *
from matplotlib import rc


################################## INPUT FILE ##################################

snapdir = '/Volumes/ELEMENTS/FULVIO/Simulations/baryon_512_'  #folder hosting the catalogue
simnames = ['lcdm','q_0', 'q_1', 'q_2', 'q_3', 'p_0', 'p_1', 'p_2','p_3']
snapnum = [4,3,2]    #redshift 0, 0.5, 1
massrange = [12, 13.98]
nbins = 10

############## READ BOX SIZE IN Mpc/h AND DETERMINE THE REDSHIFT ###############

snapshot = snapdir+simnames[0]+"/snapdir_00"+str(snapnum[0])+"/snap_00"+str(snapnum[0])
BoxSize  = readgadget.header(snapshot).boxsize/1e3
z_dict = {4:'0.0', 3:'0.5', 2:'1.0', 1:'1.5', 0:'2.0'}
z = [z_dict[n] for n in snapnum]

##################### GET GROUP MASSES FOR EACH SIMULATION #####################

groupmasses = []

for i in range(len(simnames)):
    for snap in snapnum:
        SubF = readsubf.subfind_catalog(snapdir+simnames[i] , snap,
                                        group_veldisp = True, masstab = True)
        print("Simulation: " + simnames[i])
        print("number of groups =", SubF.ngroups)
        print("number of subgroups =", SubF.nsubs)
        groupmasses.append(sorted(SubF.group_m_crit200*1e10))

###################### FILL THE BINS FOR EACH SIMULATION #######################
# The bins list contains a list for each simulation each one containing how    #
# many objects fall into the n-th bin according to the masses array, which     #
# defines the boundaries of the bins                                           #
################################################################################

masses, bins = np.logspace(massrange[0], massrange[1], nbins+1), []
for i in range(len(simnames)):
    for j_z in range(3):
        bins.append([0])
        mass_idx = 0
        for m in groupmasses[i*3+j_z]:
            if m > masses[-1]:
                break
            if m >= masses[mass_idx] and m < masses[mass_idx+1]:
                #print(bins[i][-1])
                bins[i*3+j_z][-1] += 1
            elif m >= masses[mass_idx+1]:
                while m >= masses[mass_idx+1]:
                    bins[i*3+j_z].append(0)
                    mass_idx += 1
                bins[i*3+j_z][-1] += 1
        bins[i*3+j_z] = bins[i*3+j_z] + [0]*(nbins-len(bins[i*3+j_z]))
        print("Done with " + simnames[i])

print("***")
print(bins)

################################ PLOT THE RATIO ################################
matplotlib.rcParams.update({'font.size': 13})
fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, sharex=True,
                                                         figsize=(10,12))
axes = (ax1,ax2,ax3,ax4,ax5,ax6)
plt.subplots_adjust(wspace=0.22,  hspace=0.07)

for ax in (ax1,ax2,ax3,ax4,ax5,ax6):
    ax.set_ylim([0.75,1.3])
    ax.set_xscale('log')
    #ax.set_aspect(0.4)
for ax in (ax1,ax3,ax5): ax.set_ylabel(r'$N(M)/N(M)_{\Lambda \mathrm{CDM}}$')
for ax in (ax2,ax4,ax6): ax.set_ylabel(r'$N(M)/N(M)_{\alpha_{xb}=0}$')
for ax in (ax5,ax6): ax.set_xlabel(r'$M_{200} [h^{-1}M_{\odot}]$')

for j_z in range(3):
    for i in range(len(simnames)-1,-1,-1):
        print(simnames[i],np.array(bins[i*3+j_z]))
        axes[2*j_z].plot(masses[:-1],
                         np.array(bins[i*3+j_z])/np.array(bins[0+j_z]),
                         label=labels[i], color=colors[i], linewidth=linewidths[i],
                         linestyle=linestyles[i])
        if z[j_z][2]=='5':
            txt=r'$z='+z[j_z]+"$"
        else:
            txt=r'$z='+z[j_z][0]+"$"
        axes[2*j_z].text(3e13,1.2,txt)

for j_z in range(3):
    for i in range(len(simnames)-1,0,-1):
        axes[2*j_z+1].plot(masses[:-1],
                           np.array(bins[i*3+j_z])/np.array(bins[den[i]*3+j_z]),
                           label=labels[i], color=colors[i], linewidth=linewidths[i-1],
                           linestyle=linestyles[i])
        if z[j_z][2]=='5':
            txt=r'$z='+z[j_z]+"$"
        else:
            txt=r'$z='+z[j_z][0]+"$"
        axes[2*j_z+1].text(3e13,1.2,txt)


first_legend_lines=[]
first_label = [r'$\Lambda\mathrm{CDM}$',r'$\alpha_{xb}=0$',r'$\alpha_{xb}=1$',r'$\alpha_{xb}=10$',r'$\alpha_{xb}=100$']
for i in range(5):
    first_legend_lines.append(mlines.Line2D([], [], color=pal[i], label=first_label[i]))

second_legend_lines=[]
second_label = ['$w_x=-0.9$','$w_x=-1.1$']
second_linestyle = ['-','--']
for i in range(2):
    second_legend_lines.append(mlines.Line2D([], [], color='grey',
                               label=second_label[i],
                               linestyle=second_linestyle[i]))

leg1 = ax1.legend(handles=first_legend_lines,loc="upper left",fontsize=12,ncol=2)
#ax7 = plt.gca().add_artist(first_legend)
leg2 = ax1.legend(handles=second_legend_lines, loc="lower left",fontsize=12,ncol=2)

ax1.add_artist(leg1)
ax1.add_artist(leg2)

#plt.title("Mass Function at z=" + str(redshift))


#plt.legend(ncol=2)
plt.savefig('Figures/Mass_Function.pdf', bbox_inches='tight')
plt.show()
plt.close()
