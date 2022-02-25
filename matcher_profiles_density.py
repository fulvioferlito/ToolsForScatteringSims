import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib
import seaborn as sns
from math import frexp
from decimal import Decimal
from ff_plot_options import *
import sys

############## FUNCTION TO CHECK IF TWO GROUPS ARE THE SAME GROUP ##############

def check(g_1, g_2):
    # check position (2.5*r_200)
    delta = 0.8*((g_1[8] + g_2[8])/2)
    dist = np.sqrt((g_1[5]-g_2[5])**2 + (g_1[6]-g_2[6])**2 + (g_1[7]-g_2[7])**2)
    if dist > delta:
        return False
    # check mass and radius (20%)
    for i in range(8, 10):
        if abs(g_1[8] - g_2[8])/(max(g_1[8], g_2[8])) > 0.1:
            return False
    return True

######################## FUNCTION TO REMOVE NULL VALUES ########################

def clean_gr(gr):
    bad = []
    for i in range(gr.shape[0]):
        if gr[i,1] < 1 or gr[i,1] > 5e5:
            bad.append(i)
    gr = np.delete(gr, bad, axis=0)
    return gr

################################# INPUT FILES ##################################

snapdir = '/Volumes/ELEMENTS/FULVIO/Simulations/baryon_512_'  #folder hosting the catalogue
simnames    = ['lcdm','q_0','q_1','q_2','q_3','p_0','p_1','p_2','p_3']
part_type   = [0,1]
snapnum     = 4
nbin        = 4      #starting from 1
rho_crit   = 7.72283e-09

############################### LOAD GROUPS INFO ###############################
# groupdata is a list containing one array for each simulation in simnames     #
# every array contains the data of the relative groups_properties file.        #
# The columns are the following:                                               #
# GrNr, GroupLength, SubGroupLength, GroupMass, SubGroupMass, GroupPos[0],     #
# GroupPos[1], GroupPos[2], R_200, M_200                                       #
################################################################################

n_sims, groupdata = len(simnames), []
for i in range(n_sims):
    fname = (snapdir + simnames[i] + "/profiles/groups_properties_bin_"
             + "{:02d}".format(nbin) + "_"
             + "{:03d}".format(snapnum) + ".txt")
    arr = np.loadtxt(fname,
                     dtype='int,'*3+'float,'*6+'float')
    groupdata.append(arr)

############################### LOOK FOR A MATCH ###############################
# If we find a match, we store the identities in the ids list, corresponding   #
# the  groups identity in the various simulations in the same order.           #
# The matches list will contain a sublist for every  match detected. Each      #
# sublist contains the ids list of the relative group                          #
################################################################################

matches, matches_radius, matches_mass = [], [], []
for gr_0 in groupdata[0]:
    ids = [gr_0[0]]
    gr_radius = [gr_0[8]]
    gr_mass = [gr_0[9]]
    for other_sim in groupdata[1:]:
        found = False
        for gr_i in other_sim:
            if check(gr_0, gr_i):
                ids.append(gr_i[0])
                gr_radius.append(gr_i[8])
                gr_mass.append(gr_i[9])
                found = True
                break
        if found == False:
            break
    if found == False:
        continue
    else:
        matches.append(ids)
        matches_radius.append(gr_radius)
        matches_mass.append(gr_mass)
        print("********** Match n. "+str(len(matches)-1)+" **********\n")
        print("   Simulation  -  Group ID\n")
        for i in range(n_sims):
            print(" "*5,simnames[i].ljust(6), "  -  ","{:5d}".format(ids[i]))
        print("\n")

############################ DO YOU WANT TO PLOT?? #############################

if len(matches) > 0:
    chosen = int(input("Which group do you want plot? "))
else:
    print("No matches found :(")
    sys.exit()

############################### SET PLOT OPTIONS ###############################

matplotlib.rcParams.update({'font.size': 12})
fig, ax = plt.subplots()
#plt.title("Stacked profile")
#axes.set_xlim([0.1,10.1])
ax.set_ylim([2,7e5])
ax.set_xscale('log')
ax.set_yscale('log')

##################################### PLOT #####################################

for i in range(len(simnames)):
    for ptype in part_type:
        gr = np.loadtxt(snapdir + simnames[i] + "/profiles/snap_"
                        + "{:03d}".format(snapnum) + "_group_"
                        + "{:03d}".format(matches[chosen][i])
                        + "_density_type_" + str(ptype) + ".txt")[:,:2]
        gr[:,1] /= rho_crit*matches_radius[chosen][i]**3
        gr[:,0] *= matches_radius[chosen][i]
        gr = clean_gr(gr)
        plt.plot(gr[:,0],gr[:,1], label = (1-ptype)*labels[i],
                 color=colors[i], linestyle=linestyles[i], alpha=alphas[ptype])

plt.title(r"$\mathrm{Density \; profile \; at \; } z=0 \mathrm{, \; mass \; bin \; n. \;}"+str(nbin)+"$")
plt.ylabel(r'$\rho/\rho_{\mathrm{crit}}$')
plt.xlabel('$R \,[\mathrm{kpc}/h]$')


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


leg1 = ax.legend(handles=first_legend_lines,loc="upper right")
leg2 = ax.legend(handles=second_legend_lines, loc="lower left")
leg3 = ax.legend(handles=third_legend_lines, loc="upper center")

ax.add_artist(leg1)
ax.add_artist(leg2)
ax.add_artist(leg3)

##################################### TEXT #####################################

x0, xmax = plt.xlim()
y0, ymax = plt.ylim()
data_width = xmax - x0
data_height = ymax - y0
mass = '%.2E' % Decimal(str(matches_mass[chosen][0]*1e10))
txt ="$M_{200}(\Lambda \mathrm{CDM}) = "+ mass[:4] +r" \times 10^{"+ mass[6:] +"}\, h^{-1}M_{\odot}$"
plt.text(x0 + data_width * 0.0007, y0 + data_height * 0.00007, txt, fontsize=12)



plt.savefig('Figures/Halo_profile_' + str(nbin) + '.pdf',
            bbox_inches='tight')
plt.show()
plt.close()
