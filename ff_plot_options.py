import seaborn as sns
from matplotlib import rc
############################### SET PLOT OPTIONS ###############################

#set labels
labels = ['LCDM', 'Q0', 'Q1' ,'Q10', 'Q100', 'P0', 'P1', 'P10','P100']
#set colors etc..
pal = sns.color_palette("bright")
colors = [pal[0],pal[1],pal[2],pal[3],pal[4],pal[1],pal[2],pal[3],pal[4]]
linestyles = ["-","-","-","-","-","--","--","--","--"]
linewidths = [2.8] + [1.5]*8
den = [0,1,1,1,1,5,5,5,5]
alphas = [1, 0.6]

rc('font',**{'family':'sans-serif'})
rc('text', usetex=True)
