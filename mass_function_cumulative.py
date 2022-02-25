import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import readfof
import readgadget
import readsubf
import mass_function_library as MFL
from scipy.interpolate import interp1d
import matplotlib.lines as mlines
from ff_plot_options import *

# input file
snapdir = '/Volumes/ELEMENTS/FULVIO/Simulations/baryon_512_'  #folder hosting the catalogue
simnames = ['lcdm', 'q_0', 'q_1', 'q_2', 'q_3', 'p_0', 'p_1', 'p_2','p_3']
snapnum = 4    #redshift 0

#Read box size in Mpc/h
snapshot = snapdir+simnames[0]+"/snapdir_00"+str(snapnum)+"/snap_00"+str(snapnum)
BoxSize  = readgadget.header(snapshot).boxsize/1e3

# determine the redshift of the catalogue
z_dict = {4:0.0, 3:0.5, 2:1.0, 1:1.5, 0:1.0}
redshift = z_dict[snapnum]
mass = []

for i in range(len(simnames)):
    SubF = readsubf.subfind_catalog(snapdir+simnames[i] , snapnum,
                                    group_veldisp = True, masstab = True)
    print ("number of groups =", SubF.ngroups)
    print ("number of subgroups =", SubF.nsubs)
    mass.append(sorted(i for i in SubF.group_mass*1e10 if i > 1e11))

# theoretical halo mass function parameters
print("Computing theoretical mass function ...")
f_Pk   = 'CAMB/camb_matterpower_0.dat'  #file with linear Pk
OmegaM = 0.27825
Masses = np.logspace(11, 16, 100) #array with halo masses
author = 'Jenkins'   #Sheth-Tormen halo mass function
bins   = 100  #number of bins to use for Pk
z      = 0.0    #redshift; only used for Tinker, Tinker10 and Crocce
delta  = 200.0  #overdensity; only for Tinker and Tinker10

# read linear matter Pk
k, Pk = np.loadtxt(f_Pk, unpack=True)

# compute theoretical halo mass function
HMF = MFL.MF_theory(k, Pk, OmegaM, Masses, author, bins, z, delta)

# cumulative halo mass function
CHMF = [np.trapz(HMF[i:], x=Masses[i:]) for i in range(len(HMF)-1)]

#Set plot characteristics
matplotlib.rcParams.update({'font.size': 13})
fig, (ax1, ax2) = plt.subplots(2, sharex=True,gridspec_kw={'height_ratios': [1.5, 1]},figsize=(7,9))
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlim(2e11,2e14)
ax1.set_ylim(7e-6,0.02)
ax2.set_ylim(0.8,1.2)
plt.subplots_adjust(hspace=0.1)

pal[7]=pal[5]
# Compute our cumulative mass function
print("computing cumlative mass functions ...")
for i in range(len(simnames)-1,-1,-1):
    N_tot = len(mass[i])
    N = [n/BoxSize**3 for n in range(N_tot, 0, -1)]
    x, y = mass[i], N
    f = interp1d(x, y, fill_value="extrapolate")
    ax1.plot(Masses[:-1], f(Masses[:-1]), label=labels[i],
             color = colors[i], linestyle=linestyles[i], linewidth=linewidths[i])
    ax2.plot(Masses[:-1], f(Masses[:-1])/CHMF, label=labels[i],
             color = colors[i], linestyle=linestyles[i], linewidth=linewidths[i])
    print("Done with " + simnames[i])

ax1.plot(Masses[:-1], CHMF, label="Jenkins", color=pal[7],linewidth=linewidths[0]*0.9,linestyle="-")
ax2.plot(Masses[:-1], [a/a for a in CHMF], label="Jenkins", color=pal[7],linewidth=linewidths[0]*0.9,linestyle="-")

first_legend_lines=[mlines.Line2D([], [], color=pal[7], label='$\mathrm{Jenkins}$', linestyle='-')]
first_label = ['$\Lambda\mathrm{CDM}$',r'$\alpha_{xb}=0$',r'$\alpha_{xb}=1$',r'$\alpha_{xb}=10$',r'$\alpha_{xb}=100$']
for i in range(5):
    first_legend_lines.append(mlines.Line2D([], [], color=pal[i], label=first_label[i]))

second_legend_lines=[]
second_label = ['$w_x=-0.9$','$w_x=-1.1$']
second_linestyle = ['-','--']
for i in range(2):
    second_legend_lines.append(mlines.Line2D([], [], color='grey',
                               label=second_label[i],
                               linestyle=second_linestyle[i]))


leg1 = ax1.legend(handles=first_legend_lines,loc="upper right")
#ax7 = plt.gca().add_artist(first_legend)
leg2 = ax1.legend(handles=second_legend_lines, loc="lower left")
ax1.add_artist(leg1)
ax1.add_artist(leg2)
#ax1.set_title("Cumulative Halo Mass Function at z=0")
ax1.set_ylabel(r'$N(>M)[(h^{-1}\mathrm{Mpc})^{-3}]$')
ax2.set_ylabel(r'$N(>M)/N(>M)_{\mathrm{Jenkins}}$')
ax2.set_xlabel('$M [h^{-1}M_{\odot}]$')
plt.savefig('Figures/Cumulative_Mass_Function.pdf')
plt.show()
