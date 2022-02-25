import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import seaborn as sns
from ff_plot_options import *

################################ SELECT OPTIONS ################################

basename   = "/Volumes/ELEMENTS/FULVIO/Simulations"
simnames   = ['lcdm','q_0','q_1','q_2','q_3','p_0','p_1','p_2','p_3']
snapnum    = 4
nbin       = 1
part_type  = [0,1]
rho_crit   = 7.72283e-09

################################# LOAD GROUPS ##################################
# Groups is a dictionary of dictionaries. We have a key for each simulation:   #
# the key of the subdictionary is the number of the group, the value is a list #
# of array, each array is the density profile of a component (in order).       #
################################################################################

groups = {}
for simname in simnames:
    groups[simname] = {}
    gr_info = np.loadtxt(basename + "/baryon_512_"+ simname
                         + "/profiles/groups_properties_bin_"
                         + "{:02d}".format(nbin) + "_"
                         + "{:03d}".format(snapnum) + ".txt")
    groups_names = gr_info[:,0].astype('int')
    groups_r200 = gr_info[:,8]
    for i in range(len(groups_names)):
        #print("group:",groups_names[i],"r200",groups_r200[i])
        groups[simname][groups_names[i]] = []
        for ptype in part_type:
            gr = np.loadtxt(basename + "/baryon_512_"+ simname
                            + "/profiles/snap_" + "{:03d}".format(snapnum)
                            + "_group_" + "{:03d}".format(groups_names[i])
                            + "_density_type_" + str(ptype) + ".txt")[:,:2]
            gr[:,1] /= rho_crit*groups_r200[i]**3
            groups[simname][groups_names[i]].append(gr)

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


############################### SET PLOT OPTIONS ###############################

matplotlib.rcParams.update({'font.size': 12})
fig, ax = plt.subplots(2, 1, gridspec_kw={'height_ratios': [2, 1]})
fig.subplots_adjust(hspace=0)
fig.set_figheight(7)
#plt.title("Stacked profile")
#axes.set_xlim([0.1,10.1])
ax[0].set_ylim([2,7e5])
ax[1].set_xscale('log')
ax[0].set_xscale('log')
ax[0].set_yscale('log')

##################################### PLOT #####################################

for i in range(len(simnames)-1,-1,-1):
    st_pr_lcdm_tot = stacked_profiles[simnames[0]][0][:,1] + stacked_profiles[simnames[0]][1][:,1]
    st_pr_tot = stacked_profiles[simnames[i]][0][:,1] + stacked_profiles[simnames[i]][1][:,1]
    for ptype in part_type:
        st_pr = stacked_profiles[simnames[i]][ptype]
        ax[0].plot(st_pr[:,0][1:],st_pr[:,1][1:], label = (1-ptype)*labels[i],
                   color=colors[i], linestyle=linestyles[i], alpha=alphas[ptype],
                   linewidth=linewidths[i])
    ptype = 0
    ax[1].plot(st_pr[:,0][1:],st_pr_tot[1:]/st_pr_lcdm_tot[1:], label = (1-ptype)*labels[i],
               color=colors[i], linestyle=linestyles[i], alpha=alphas[ptype],
               linewidth=linewidths[i])

ax[0].set_title("$\mathrm{Stacked \; density \; profile \; of \; 100 \; halos \; at \;} z=0$"+"\n"+"$ \mathrm{Mass \; bin \; n. \;} "+str(nbin)+"$")
ax[0].set_ylabel(r'$\rho/\rho_{\mathrm{crit}}$')
ax[1].set_ylabel(r'$(\rho/\rho_{\mathrm{crit}})/(\rho/\rho_{\mathrm{crit}})_{\Lambda \mathrm{CDM}}$')
ax[1].set_xlabel('$R/R_{200}$')

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

third_legend_lines=[]
third_label = ['$\mathrm{Dark \; Matter}$','$\mathrm{Baryons}$']
third_alpha = alphas[::-1]
for i in range(2):
    third_legend_lines.append(mlines.Line2D([], [], color=pal[0],
                               label=third_label[i],
                               alpha=third_alpha[i]))


leg1 = ax[0].legend(handles=first_legend_lines,loc="upper right")
leg2 = ax[0].legend(handles=second_legend_lines, loc="lower left")
leg3 = ax[0].legend(handles=third_legend_lines, loc="upper center")

ax[0].add_artist(leg1)
ax[0].add_artist(leg2)
ax[0].add_artist(leg3)


##################################### TEXT #####################################

x0, xmax = ax[0].get_xlim()
y0, ymax = ax[0].get_ylim()
data_width = xmax - x0
data_height = ymax - y0
m = [r"5 \times 10^{12}", r"1 \times 10^{13}", r"5 \times 10^{13}", r"1 \times 10^{14}", r"5 \times 10^{14}"]
txt = "$"+m[nbin-1]+"\, h^{-1}M_{\odot} < M_{200} < "+m[nbin]+"\, h^{-1}M_{\odot}$"
ax[0].text(x0 + data_width * 0.0007, y0 + data_height * 0.00007, txt, fontsize=12)

plt.savefig('Figures/Stacked_profile_alt_tot_'+str(nbin)+'.pdf', bbox_inches='tight')
plt.show()
plt.close()
