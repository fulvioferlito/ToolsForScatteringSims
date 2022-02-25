import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import seaborn as sns
from ff_plot_options import *

################################ SELECT OPTIONS ################################

basename   = "/Volumes/ELEMENTS/FULVIO/Simulations"
simnames    =   ['lcdm','q_0','q_1','q_2','q_3','p_0','p_1','p_2','p_3']
snapnum     =   4
nbin        =   4
softening   =   12     #kpc  (we cut the profiles at radii r < 3*softening)
part_type   =   [0,1]

################################# LOAD GROUPS ##################################
# Groups is a dictionary of dictionaries. We have a key for each simulation:   #
# the key of the subdictionary is the number of the group, the value is a list #
# of array, each array is the density profile of a component (in order).       #
################################################################################

groups = {}
r_min  = 1e10  #minimum radius (initialized with a deliberately high value)
for simname in simnames:
    groups[simname] = {}
    gr_file = (basename + "/baryon_512_"+ simname
               + "/profiles/groups_properties_bin_" + "{:02d}".format(nbin) + "_"
               + "{:03d}".format(snapnum) + ".txt")
    gr_nm = np.loadtxt(gr_file)[:,0]
    groups_names = gr_nm.astype('int')
    r = min(np.loadtxt(gr_file)[:,8])
    if r < r_min:
        r_min = r
    for group_name in groups_names:
        groups[simname][group_name] = []
        for ptype in part_type:
            gr = np.loadtxt(basename + "/baryon_512_"+ simname
                            + "/profiles/snap_" + "{:03d}".format(snapnum)
                            + "_group_" + "{:03d}".format(group_name)
                            + "_density_type_" + str(ptype) + ".txt")[:,:2]
            groups[simname][group_name].append(gr)

########################### CREATE STACKED PROFILES ############################
# The stacked_profiles dictionary contains a list of arrays for each           #
# simulation. Each list contain one array for every particle type.             #
################################################################################

stacked_profiles = {}
for simname in simnames:
    stacked_profiles[simname] = []
    for ptype in part_type:
        st_pr = []
        for group in groups[simname].keys():
            if len(st_pr) == 0:
                st_pr = groups[simname][group][ptype]
            else:
                st_pr[:,1] += groups[simname][group][ptype][:,1]
        st_pr[:,1] /= len(groups[simname].keys())
        stacked_profiles[simname].append(st_pr)

################################### IDX CUT ####################################
# We compute the minimum radius to plot so that is exceeds three times the     #
# softening lenght                                                             #
################################################################################

r_2 = np.loadtxt(basename + "/baryon_512_"+ simname
                + "/profiles/snap_" + "{:03d}".format(snapnum)
                + "_group_" + "{:03d}".format(group_name)
                + "_density_type_0.txt")[:,2]
idx_min = 0
for i in range(len(r_2)):
    print(r_min*r_2[i] , 3*softening)
    if r_min*r_2[i] < 3*softening:
        idx_min = i
    else:
        break
print(idx_min)

############################### SET PLOT OPTIONS ###############################

matplotlib.rcParams.update({'font.size': 12})
fig, ax = plt.subplots()
ax.set_ylim([-0.01,0.31])
#ax.set_xlim([0.01,2.2])
ax.set_xscale('log')

##################################### PLOT #####################################

for i in range(len(simnames)-1,-1,-1):
    r_2 = np.loadtxt(basename + "/baryon_512_"+ simname
                    + "/profiles/snap_" + "{:03d}".format(snapnum)
                    + "_group_" + "{:03d}".format(group_name)
                    + "_density_type_0.txt")[:,2]
    rho_b = stacked_profiles[simnames[i]][0][:,1]
    rho_c = stacked_profiles[simnames[i]][1][:,1]
    cum_m_b = [4*np.pi*r_2[0]**3*rho_b[0]/3]
    cum_m_c = [4*np.pi*r_2[0]**3*rho_c[0]/3]
    for j in range(1,len(r_2)):
        cum_m_b.append(cum_m_b[-1]+4*np.pi*(r_2[j]**3-r_2[j-1]**3)*rho_b[j]/3)
        cum_m_c.append(cum_m_c[-1]+4*np.pi*(r_2[j]**3-r_2[j-1]**3)*rho_c[j]/3)
    baryon_fraction = [cum_m_b[i]/cum_m_c[i] for i in range(len(r_2))]
    plt.plot(r_2[idx_min:], baryon_fraction[idx_min:], color=colors[i], linestyle=linestyles[i],
             alpha=alphas[0], label=labels[i], linewidth=linewidths[i])

plt.title("$\mathrm{Stacked \; baryon \: fraction \; profile \; of \; 100 \; halos \; at \;} z=0$"+"\n"+"$ \mathrm{Mass \; bin \; n. \;} "+str(nbin)+"$")
plt.ylabel(r'$f_b(R)$')
plt.xlabel('$R/R_{200}$')

#################################### LEGEND ####################################

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

leg1 = ax.legend(handles=first_legend_lines,loc="lower right")
leg2 = ax.legend(handles=second_legend_lines, loc="upper right")

ax.add_artist(leg1)
ax.add_artist(leg2)

##################################### TEXT #####################################

m = [r"5 \times 10^{12}", r"1 \times 10^{13}", r"5 \times 10^{13}", r"1 \times 10^{14}", r"5 \times 10^{14}"]
txt = "$"+m[nbin-1]+"\, h^{-1}M_{\odot} < M_{200} < "+m[nbin]+"\, h^{-1}M_{\odot}$"
plt.text(r_2[idx_min], 0.01, txt, fontsize=12)


plt.savefig('Figures/Baryon_fraction_stacked_profile_'+str(nbin)+'.pdf', bbox_inches='tight')
plt.show()
plt.close()
