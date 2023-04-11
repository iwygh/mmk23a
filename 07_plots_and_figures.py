def circular_hist(ax, x, edgecolor,bins=16, density=True, offset=0, gaps=True):
    """
    https://stackoverflow.com/questions/22562364/circular-polar-histogram-in-python
    Produce a circular histogram of angles on ax.
    Parameters
    ----------
    ax : matplotlib.axes._subplots.PolarAxesSubplot
        axis instance created with subplot_kw=dict(projection='polar').
    x : array
        Angles to plot, expected in units of radians.
    bins : int, optional
        Defines the number of equal-width bins in the range. The default is 16.
    density : bool, optional
        If True plot frequency proportional to area. If False plot frequency
        proportional to radius. The default is True.
    offset : float, optional
        Sets the offset for the location of the 0 direction in units of
        radians. The default is 0.
    gaps : bool, optional
        Whether to allow gaps between bins. When gaps = False the bins are
        forced to partition the entire [-pi, pi] range. The default is True.
    Returns
    -------
    n : array or list of arrays
        The number of values in each bin.
    bins : array
        The edges of the bins.
    patches : `.BarContainer` or list of a single `.Polygon`
        Container of individual artists used to create the histogram
        or list of such containers if there are multiple input datasets.
    """
    # Wrap angles to [-pi, pi)
    # x = (x+np.pi) % (2*np.pi) - np.pi
    x = np.mod(x+np.pi,2*np.pi) - np.pi
    # Force bins to partition entire circle
    if not gaps:
        bins = np.linspace(-np.pi, np.pi, num=bins+1)
    # Bin data and record counts
    n, bins = np.histogram(x, bins=bins)
    # Compute width of each bin
    widths = np.diff(bins)
    # By default plot frequency proportional to area
    if density:
        # Area to assign each bin
        area = n / x.size
        # Calculate corresponding bin radius
        radius = (area/np.pi) ** .5
    # Otherwise plot frequency proportional to radius
    else:
        radius = n
    # Plot data on ax
    patches = ax.bar(bins[:-1], radius, zorder=1, align='edge', width=widths,
                     edgecolor=edgecolor, fill=False, linewidth=1)
    # Set the direction of the zero angle
    ax.set_theta_offset(offset)
    # Remove ylabels for area plots (they are mostly obstructive)
    if density:
        ax.set_yticks([])
    return n, bins, patches
#%%
def JDfun(year,month,day,hour,minute,second):
    # from Vallado pg 183, valid for yrs 1900 to 2100
    import numpy as np
    JD = 367*year - np.floor(7/4*(year+np.floor(1/12*(month+9)))) + \
        np.floor(275*month/9) + day + 1721013.5 + 1/24*(hour+1/60*(second/60+minute))
    return JD
#%%
# def truncated_rayleigh_pdf(x, sigma):
def truncated_rayleigh_pdf(x,kappa):
    import numpy as np
    sigma = 1/np.sqrt(kappa)
    C_R = 1/(1-np.exp(-np.pi**2/2/sigma**2))
    return C_R/sigma**2 * x * np.exp(-x**2/2/sigma**2)
#%%
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from astroquery.jplhorizons import Horizons
#%%
df = pd.read_csv('calculations_dictionary.csv')
df_mean = pd.read_csv('vm17_mean_plane_clones_noclones.csv',delim_whitespace=True)
theta_mid_deg_clones = df_mean['theta_mid_deg'][0]
phi_mid_deg_clones = df_mean['phi_mid_deg'][0]
theta_mid_deg_noclones = df_mean['theta_mid_deg'][1]
phi_mid_deg_noclones = df_mean['phi_mid_deg'][1]
Rbar_clones = df['Rbar_clones'][0]
Rbar_noclones = df['Rbar_noclones'][0]
alphaconfidence = df['alphaconfidence'][0]
xhat_clones = df['xhat_clones'][0]
yhat_clones = df['yhat_clones'][0]
zhat_clones = df['zhat_clones'][0]
xhat_noclones = df['xhat_noclones'][0]
yhat_noclones = df['yhat_noclones'][0]
zhat_noclones = df['zhat_noclones'][0]
theta_mid_clones = np.radians(theta_mid_deg_clones)
phi_mid_clones = np.radians(phi_mid_deg_clones)
theta_mid_noclones = np.radians(theta_mid_deg_noclones)
phi_mid_noclones = np.radians(phi_mid_deg_noclones)
xhat_vm17_clones = np.sin(theta_mid_clones)*np.cos(phi_mid_clones)
yhat_vm17_clones = np.sin(theta_mid_clones)*np.sin(phi_mid_clones)
zhat_vm17_clones = np.cos(theta_mid_clones)
xhat_vm17_noclones = np.sin(theta_mid_noclones)*np.cos(phi_mid_noclones)
yhat_vm17_noclones = np.sin(theta_mid_noclones)*np.sin(phi_mid_noclones)
zhat_vm17_noclones = np.cos(theta_mid_noclones)
plutinos_file_clones = 'plutinos_for_mnras_clones.csv'
df_clones = pd.read_csv(plutinos_file_clones)
n_clones = df_clones.shape[0]
des_clones = df_clones['Packed MPC designation'].tolist()
aPb_au_clones = df_clones['Semimajor axis au barycentric'].tolist()
ePb_clones = df_clones['Eccentricity barycentric'].tolist()
iPb_deg_clones = df_clones['Inclination ecliptic J2000 barycentric degrees'].tolist()
WPb_deg_clones = df_clones['Longitude of ascending node ecliptic J2000 barycentric degrees'].tolist()
wPb_deg_clones = df_clones['Argument of perihelion ecliptic J2000 barycentric degrees'].tolist()
MPb_deg_clones = df_clones['Mean anomaly ecliptic J2000 barycentric degrees'].tolist()
aPb_au_clones = np.array(aPb_au_clones)
ePb_clones = np.array(ePb_clones)
iPb_rad_clones = np.radians(np.array(iPb_deg_clones))
WPb_rad_clones = np.radians(np.array(WPb_deg_clones))
wPb_rad_clones = np.radians(np.array(wPb_deg_clones))
MPb_rad_clones = np.radians(np.array(MPb_deg_clones))
x_clones =  np.sin(iPb_rad_clones)*np.sin(WPb_rad_clones) # x component of orbit normal vector
y_clones = -np.sin(iPb_rad_clones)*np.cos(WPb_rad_clones) # y component of orbit normal vector
z_clones =  np.cos(iPb_rad_clones) # z component of orbit normal vector
plutinos_file_noclones = 'plutinos_for_mnras_noclones.csv'
df_noclones = pd.read_csv(plutinos_file_noclones)
n_noclones = df_noclones.shape[0]
des_noclones = df_noclones['Packed MPC designation'].tolist()
aPb_au_noclones = df_noclones['Semimajor axis au barycentric'].tolist()
ePb_noclones = df_noclones['Eccentricity barycentric'].tolist()
iPb_deg_noclones = df_noclones['Inclination ecliptic J2000 barycentric degrees'].tolist()
WPb_deg_noclones = df_noclones['Longitude of ascending node ecliptic J2000 barycentric degrees'].tolist()
wPb_deg_noclones = df_noclones['Argument of perihelion ecliptic J2000 barycentric degrees'].tolist()
MPb_deg_noclones = df_noclones['Mean anomaly ecliptic J2000 barycentric degrees'].tolist()
aPb_au_noclones = np.array(aPb_au_noclones)
ePb_noclones = np.array(ePb_noclones)
iPb_rad_noclones = np.radians(np.array(iPb_deg_noclones))
WPb_rad_noclones = np.radians(np.array(WPb_deg_noclones))
wPb_rad_noclones = np.radians(np.array(wPb_deg_noclones))
MPb_rad_noclones = np.radians(np.array(MPb_deg_noclones))
x_noclones =  np.sin(iPb_rad_noclones)*np.sin(WPb_rad_noclones) # x component of orbit normal vector
y_noclones = -np.sin(iPb_rad_noclones)*np.cos(WPb_rad_noclones) # y component of orbit normal vector
z_noclones =  np.cos(iPb_rad_noclones) # z component of orbit normal vector
#%% section 2.5 estimating a confidence region for the mean direction
d_vm17_clones = 1 - 1/n_clones * np.sum( (x_clones*xhat_vm17_clones + \
                    y_clones*yhat_vm17_clones + z_clones*zhat_vm17_clones)**2 )
sigmahat_vm17_clones = np.sqrt(d_vm17_clones/(n_clones*Rbar_clones**2))
sigmahat_vm17_clones = sigmahat_vm17_clones
sigmahat_degrees_vm17_clones = np.degrees(sigmahat_vm17_clones)
q_vm17_clones = np.arcsin(sigmahat_vm17_clones*np.sqrt(-np.log(alphaconfidence)))
q_vm17_clones = q_vm17_clones
q_vm17_degrees_clones = np.degrees(q_vm17_clones)
d_vm17_noclones = 1 - 1/n_noclones * np.sum( (x_noclones*xhat_noclones + \
                    y_noclones*yhat_noclones + z_noclones*zhat_noclones)**2 )
sigmahat_vm17_noclones = np.sqrt(d_vm17_noclones/(n_noclones*Rbar_noclones**2))
sigmahat_vm17_noclones = sigmahat_vm17_noclones
sigmahat_degrees_vm17_noclones = np.degrees(sigmahat_vm17_noclones)
q_vm17_noclones = np.arcsin(sigmahat_vm17_noclones*np.sqrt(-np.log(alphaconfidence)))
q_vm17_noclones = q_vm17_noclones
q_vm17_degrees_noclones = np.degrees(q_vm17_noclones)
#%% FIGURE 1 plot normalization parameter C_F(kappa)
fig = plt.figure(figsize=(2,2))
ax1  = fig.add_subplot(111)
ax1.set_xlabel('$\kappa$')
ax1.set_ylabel('$C_F(\kappa)$')
plt.rcParams['font.size'] = 6
kappavec = np.linspace(start=0.1,stop=50,num=51,endpoint=True)
CF_vec = kappavec / (2*np.pi*(np.exp(kappavec)-np.exp(-kappavec)) )
g = ax1.semilogy(kappavec,CF_vec,color='blue')
titlestr = 'fig1_CF_kappa'
plt.savefig(titlestr + '.eps',dpi=300,bbox_inches='tight',pad_inches=0)
plt.show()
#%% FIGURE 2 plot comparison of vmf inclination and truncated rayleigh inclination
kappavec = np.array([30,15,4,1])
colors = ['green','brown','magenta','blue']
radians_vec = np.linspace(start=0,stop=np.pi,num=1000,endpoint=False)
fig = plt.figure(figsize=(3,2))
ax1  = fig.add_subplot(111)
ax1.set_xlabel('Inclination (radians)')
ax1.set_ylabel('pdf')
ax1.set_xticks([0,np.pi/4,np.pi/2,3*np.pi/4,np.pi])
ax1.set_xticklabels(['0','$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$'])
plt.rcParams['font.size'] = 6
axins = inset_axes(ax1,width="50%",height="50%",borderpad=1)
axins.set_xlabel('Inclination (radians)')
axins.set_ylabel('Rayleigh - vMF')
axins.set_xticks([0,np.pi/4,np.pi/2,3*np.pi/4,np.pi])
axins.set_xticklabels(['0','$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$'])
mindiff = []
maxdiff = []
maxy = []
for ik in range(4):
    kappa_here = kappavec[ik]
    color = colors[ik]
    sigma_here = 1/np.sqrt(kappa_here)
    curve_vec_exact = kappa_here/(np.exp(kappa_here)-np.exp(-kappa_here))*\
        np.exp(kappa_here*np.cos(radians_vec))*np.sin(radians_vec)
    C_R = (1-np.exp(-np.pi**2/(2*sigma_here**2)))**(-1)
    curve_vec_truncated = C_R/sigma_here**2 * radians_vec * \
        np.exp(-radians_vec**2/2/sigma_here**2)
    h = ax1.plot(radians_vec[0:-1:10],curve_vec_exact[0:-1:10],color=color,zorder=1,lw=0.5,linestyle='solid')
    h2 = ax1.plot(radians_vec[0:-1:10],curve_vec_truncated[0:-1:10],color=color,zorder=2,lw=0.5,linestyle='dotted')
    h3 = axins.plot(radians_vec[0:-1:10],curve_vec_truncated[0:-1:10]-curve_vec_exact[0:-1:10],color=color,\
                    lw=0.5,linestyle='solid')
    mindiff.append(np.min(curve_vec_truncated-curve_vec_exact))
    maxdiff.append(np.max(curve_vec_truncated-curve_vec_exact))
    maxy.append(1.1*np.max([np.max(curve_vec_exact),np.max(curve_vec_truncated)]))
ax1.set_xlim([0,np.pi])
ax1.set_ylim([0,np.max(maxy)])
axins.set_xlim([0,np.pi])
titlestr = 'fig2_fourcurveswithinset'
plt.savefig(titlestr + '.eps',dpi=300,bbox_inches='tight',pad_inches=0)
plt.show()
#%% FIGURE 3 3d rendering of orbit poles on the unit sphere
u, v = np.mgrid[0:np.radians(360):40j, 0:np.radians(90):40j]
x = np.cos(u) * np.sin(v)
y = np.sin(u) * np.sin(v)
z = np.cos(v)
scale = 1
x = x * scale
y = y * scale
z = z * scale
fig = plt.figure(figsize=(2.5,2.5))
plt.rcParams['font.size'] = 8
ax = fig.add_subplot(111,projection='3d')
ax.plot_surface(x,y,z,rstride=1,cstride=1,color='whitesmoke',alpha=0.1,edgecolor='gray',linewidth=0.25) # unit sphere
ax.scatter(x_clones,y_clones,z_clones,color='tomato',s=2) # plutino poles
ax.set_box_aspect((1,1,0.5))
ax.set_xticks(ticks=[-0.5,0.5], minor=False)
ax.set_yticks(ticks=[-0.5,0,0.5], minor=False)
ax.set_zticks(ticks=[0,0.5,1], minor=False)
ax.view_init(60, 45)
plt.tight_layout()
ax.set_xlabel('$h_x$ = sin(i) sin(Ω)')
ax.set_ylabel('$h_y$ = -sin(i) cos(Ω)')
ax.set_zlabel('$h_z$ = cos(i)')
titlestr = 'fig3unitsphere_clones'
plt.tight_layout()
plt.savefig(titlestr + '.pdf',dpi=300,bbox_inches='tight',pad_inches=0)
plt.show()
fig = plt.figure(figsize=(2.5,2.5))
plt.rcParams['font.size'] = 8
ax = fig.add_subplot(111,projection='3d')
ax.plot_surface(x,y,z,rstride=1,cstride=1,color='whitesmoke',alpha=0.1,edgecolor='gray',linewidth=0.25) # unit sphere
ax.scatter(x_noclones,y_noclones,z_noclones,color='tomato',s=2) # plutino poles
ax.set_box_aspect((1,1,0.5))
ax.set_xticks(ticks=[-0.5,0.5], minor=False)
ax.set_yticks(ticks=[-0.5,0,0.5], minor=False)
ax.set_zticks(ticks=[0,0.5,1], minor=False)
ax.view_init(60, 45)
plt.tight_layout()
ax.set_xlabel('$h_x$ = sin(i) sin(Ω)')
ax.set_ylabel('$h_y$ = -sin(i) cos(Ω)')
ax.set_zlabel('$h_z$ = cos(i)')
titlestr = 'fig3unitsphere_noclones'
plt.tight_layout()
plt.savefig(titlestr + '.pdf',dpi=300,bbox_inches='tight',pad_inches=0)
plt.show()
#%% FIGURE 4 plot orbit poles in (hx,hy) plane
fig = plt.figure(figsize=(2.5,2.5))
plt.rcParams['font.size'] = 8
ax = fig.add_subplot(111)
th = np.linspace(start=0,stop=2*np.pi,num=100,endpoint=True)
costh = np.cos(th)
sinth=  np.sin(th)
ax.plot(np.cos(np.radians(90))*costh,np.cos(np.radians(90))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(80))*costh,np.cos(np.radians(80))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(70))*costh,np.cos(np.radians(70))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(60))*costh,np.cos(np.radians(60))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(50))*costh,np.cos(np.radians(50))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(40))*costh,np.cos(np.radians(40))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(30))*costh,np.cos(np.radians(30))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(20))*costh,np.cos(np.radians(20))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(10))*costh,np.cos(np.radians(10))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(0))*costh,np.cos(np.radians(0))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(0*30))],[0,np.sin(np.radians(0*30))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(1*30))],[0,np.sin(np.radians(1*30))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(2*30))],[0,np.sin(np.radians(2*30))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(3*30))],[0,np.sin(np.radians(3*30))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(4*30))],[0,np.sin(np.radians(4*30))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(5*30))],[0,np.sin(np.radians(5*30))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(6*30))],[0,np.sin(np.radians(6*30))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(7*30))],[0,np.sin(np.radians(7*30))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(8*30))],[0,np.sin(np.radians(8*30))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(9*30))],[0,np.sin(np.radians(9*30))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(10*30))],[0,np.sin(np.radians(10*30))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(11*30))],[0,np.sin(np.radians(11*30))],color='gray',linestyle='-',linewidth=0.25)
ax.axhline(color='black',linestyle='-',linewidth=1)
ax.axvline(color='black',linestyle='-',linewidth=1)
ax.scatter(x_clones,y_clones,color='tomato',s=1) # plutino poles
plt.tight_layout()
plt.axis('equal')
ax.set_xlabel('$h_x$')
ax.set_ylabel('$h_y$')
ax.set_xlim([-1,1])
ax.set_ylim([-1,1])
ax.set_box_aspect(1)
titlestr = 'fig4plothxhy_clones'
plt.tight_layout()
plt.savefig(titlestr + '.eps',dpi=300,bbox_inches='tight',pad_inches=0)
plt.show()
fig = plt.figure(figsize=(2.5,2.5))
plt.rcParams['font.size'] = 8
ax = fig.add_subplot(111)
th = np.linspace(start=0,stop=2*np.pi,num=100,endpoint=True)
costh = np.cos(th)
sinth=  np.sin(th)
ax.plot(np.cos(np.radians(90))*costh,np.cos(np.radians(90))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(80))*costh,np.cos(np.radians(80))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(70))*costh,np.cos(np.radians(70))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(60))*costh,np.cos(np.radians(60))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(50))*costh,np.cos(np.radians(50))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(40))*costh,np.cos(np.radians(40))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(30))*costh,np.cos(np.radians(30))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(20))*costh,np.cos(np.radians(20))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(10))*costh,np.cos(np.radians(10))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(0))*costh,np.cos(np.radians(0))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(0*30))],[0,np.sin(np.radians(0*30))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(1*30))],[0,np.sin(np.radians(1*30))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(2*30))],[0,np.sin(np.radians(2*30))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(3*30))],[0,np.sin(np.radians(3*30))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(4*30))],[0,np.sin(np.radians(4*30))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(5*30))],[0,np.sin(np.radians(5*30))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(6*30))],[0,np.sin(np.radians(6*30))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(7*30))],[0,np.sin(np.radians(7*30))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(8*30))],[0,np.sin(np.radians(8*30))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(9*30))],[0,np.sin(np.radians(9*30))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(10*30))],[0,np.sin(np.radians(10*30))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(11*30))],[0,np.sin(np.radians(11*30))],color='gray',linestyle='-',linewidth=0.25)
ax.axhline(color='black',linestyle='-',linewidth=1)
ax.axvline(color='black',linestyle='-',linewidth=1)
ax.scatter(x_noclones,y_noclones,color='tomato',s=1) # plutino poles
plt.tight_layout()
plt.axis('equal')
ax.set_xlabel('$h_x$')
ax.set_ylabel('$h_y$')
ax.set_xlim([-1,1])
ax.set_ylim([-1,1])
ax.set_box_aspect(1)
titlestr = 'fig4plothxhy_noclones'
plt.tight_layout()
plt.savefig(titlestr + '.eps',dpi=300,bbox_inches='tight',pad_inches=0)
plt.show()
#%% FIGURE 5 detail near the origin
fig = plt.figure(figsize=(2.5,2.5))
plt.rcParams['font.size'] = 8
ax = fig.add_subplot(111)
th = np.linspace(start=0,stop=2*np.pi,num=100,endpoint=True)
costh = np.cos(th)
sinth=  np.sin(th)
ax.plot(np.cos(np.radians(80))*costh,np.cos(np.radians(80))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(70))*costh,np.cos(np.radians(70))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(60))*costh,np.cos(np.radians(60))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(50))*costh,np.cos(np.radians(50))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(40))*costh,np.cos(np.radians(40))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(30))*costh,np.cos(np.radians(30))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(20))*costh,np.cos(np.radians(20))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(10))*costh,np.cos(np.radians(10))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(90-9))*costh,np.cos(np.radians(90-9))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(90-8))*costh,np.cos(np.radians(90-8))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(90-7))*costh,np.cos(np.radians(90-7))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(90-6))*costh,np.cos(np.radians(90-6))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(90-5))*costh,np.cos(np.radians(90-5))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(90-4))*costh,np.cos(np.radians(90-4))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(90-3))*costh,np.cos(np.radians(90-3))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(90-2))*costh,np.cos(np.radians(90-2))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(90-1))*costh,np.cos(np.radians(90-1))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(10))],[0,np.sin(np.radians(10))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(20))],[0,np.sin(np.radians(20))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(30))],[0,np.sin(np.radians(30))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(40))],[0,np.sin(np.radians(40))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(50))],[0,np.sin(np.radians(50))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(60))],[0,np.sin(np.radians(60))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(70))],[0,np.sin(np.radians(70))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(80))],[0,np.sin(np.radians(80))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(90))],[0,np.sin(np.radians(90))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(100))],[0,np.sin(np.radians(100))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(110))],[0,np.sin(np.radians(110))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(120))],[0,np.sin(np.radians(120))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(130))],[0,np.sin(np.radians(130))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(140))],[0,np.sin(np.radians(140))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(150))],[0,np.sin(np.radians(150))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(160))],[0,np.sin(np.radians(160))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(170))],[0,np.sin(np.radians(170))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(180))],[0,np.sin(np.radians(180))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(190))],[0,np.sin(np.radians(190))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(200))],[0,np.sin(np.radians(200))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(210))],[0,np.sin(np.radians(210))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(220))],[0,np.sin(np.radians(220))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(230))],[0,np.sin(np.radians(230))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(240))],[0,np.sin(np.radians(240))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(250))],[0,np.sin(np.radians(250))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(260))],[0,np.sin(np.radians(260))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(270))],[0,np.sin(np.radians(270))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(280))],[0,np.sin(np.radians(280))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(290))],[0,np.sin(np.radians(290))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(300))],[0,np.sin(np.radians(300))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(310))],[0,np.sin(np.radians(310))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(320))],[0,np.sin(np.radians(320))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(330))],[0,np.sin(np.radians(330))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(340))],[0,np.sin(np.radians(340))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(350))],[0,np.sin(np.radians(350))],color='gray',linestyle='-',linewidth=0.25)
ax.axhline(color='black',linestyle='-',linewidth=1)
ax.axvline(color='black',linestyle='-',linewidth=1)
ax.scatter(x_clones,y_clones,color='tomato',s=1) # plutino poles
Npts_circle = 100
circle_clock_angles = np.linspace(start=0,stop=2*np.pi,num=Npts_circle,endpoint=True)
q_clones = df['q_clones'][0]
radius_circle = np.sin(q_clones)
hx_circle_pre = np.cos(circle_clock_angles)*radius_circle
hy_circle_pre = np.sin(circle_clock_angles)*radius_circle
hz_circle_pre = np.cos(q_clones)*np.ones(Npts_circle)
# rotate the confidence circle around to center on the vMF mean pole
alphahat_clones = df['alphahat_clones'][0]
betahat_clones = df['betahat_clones'][0]
kappahat_clones = df['kappahat_clones'][0]
A_vmf_clones = np.array([[np.cos(alphahat_clones)*np.cos(betahat_clones),np.cos(alphahat_clones)*np.sin(betahat_clones),-np.sin(alphahat_clones)],\
                          [-np.sin(betahat_clones),np.cos(betahat_clones),0],\
                          [np.sin(alphahat_clones)*np.cos(betahat_clones),np.sin(alphahat_clones)*np.sin(betahat_clones),np.cos(alphahat_clones)]])
hvec_pre = np.array([hx_circle_pre,hy_circle_pre,hz_circle_pre])
hvec_post = np.matmul(np.transpose(A_vmf_clones),hvec_pre)
hx_circle_post = hvec_post[0,:]
hy_circle_post = hvec_post[1,:]
hz_circle_post = hvec_post[2,:]
ax.plot(hx_circle_post,hy_circle_post,color='blue',linestyle='-',linewidth=1) # confidence circle
ax.scatter(xhat_clones,yhat_clones,color='blue',s=20,marker='o') # vmf midplane of plutinos
A_vm17_clones = np.array([[np.cos(theta_mid_clones)*np.cos(phi_mid_clones),np.cos(theta_mid_clones)*np.sin(phi_mid_clones),-np.sin(theta_mid_clones)],\
                          [-np.sin(phi_mid_clones),np.cos(phi_mid_clones),0],\
                          [np.sin(theta_mid_clones)*np.cos(phi_mid_clones),np.sin(theta_mid_clones)*np.sin(phi_mid_clones),np.cos(theta_mid_clones)]])
hvec_post = np.matmul(np.transpose(A_vm17_clones),hvec_pre)
hx_circle_post = hvec_post[0,:]
hy_circle_post = hvec_post[1,:]
hz_circle_post = hvec_post[2,:]
ax.plot(hx_circle_post,hy_circle_post,color='forestgreen',linestyle='-',linewidth=1) # confidence circle
ax.scatter(xhat_vm17_clones,yhat_vm17_clones,color='forestgreen',s=20,marker='D') # vm17 midplane of plutinos
i_invariable = np.radians(1.578694)
W_invariable = np.radians(107.582222)
hx_invariable = np.sin(i_invariable)*np.sin(W_invariable)
hy_invariable = -np.sin(i_invariable)*np.cos(W_invariable)
hz_invariable = np.cos(i_invariable)
ax.scatter(hx_invariable,hy_invariable,color='red',s=20,marker='v') # invariable pole of the solar system
date = '20221012'
year = 2022
month = 1
day = 1
hour = 0
minute = 0
second = 0
JD = JDfun(year,month,day,hour,minute,second)
JD = str(JD)
name = '8' # neptune barycenter
center = '500@0' # solar system barycenter
obj = Horizons(id=name,location=center,epochs=JD)
el = obj.elements()
i_neptune = np.radians(float(el['incl']))
W_neptune = np.radians(float(el['Omega']))
# equation 1 in the paper
hx_neptune = np.sin(i_neptune)*np.sin(W_neptune)
hy_neptune = -np.sin(i_neptune)*np.cos(W_neptune)
hz_neptune = np.cos(i_neptune)
ax.scatter(hx_neptune,hy_neptune,color='blue',s=20,marker='>') # orbit pole of neptune
plt.tight_layout()
plt.axis('equal')
ax.set_xlabel('$h_x$')
ax.set_ylabel('$h_y$')
ax.set_xlim([-0.07,0.09])
ax.set_ylim([-0.07,0.09])
ax.set_box_aspect(1)
titlestr = 'fig5plothxhydetail_clones'
plt.tight_layout()
plt.savefig(titlestr + '.eps',dpi=300,bbox_inches='tight',pad_inches=0)
plt.show()
fig = plt.figure(figsize=(2.5,2.5))
plt.rcParams['font.size'] = 8
ax = fig.add_subplot(111)
th = np.linspace(start=0,stop=2*np.pi,num=100,endpoint=True)
costh = np.cos(th)
sinth=  np.sin(th)
ax.plot(np.cos(np.radians(80))*costh,np.cos(np.radians(80))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(70))*costh,np.cos(np.radians(70))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(60))*costh,np.cos(np.radians(60))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(50))*costh,np.cos(np.radians(50))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(40))*costh,np.cos(np.radians(40))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(30))*costh,np.cos(np.radians(30))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(20))*costh,np.cos(np.radians(20))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(10))*costh,np.cos(np.radians(10))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(90-9))*costh,np.cos(np.radians(90-9))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(90-8))*costh,np.cos(np.radians(90-8))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(90-7))*costh,np.cos(np.radians(90-7))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(90-6))*costh,np.cos(np.radians(90-6))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(90-5))*costh,np.cos(np.radians(90-5))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(90-4))*costh,np.cos(np.radians(90-4))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(90-3))*costh,np.cos(np.radians(90-3))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(90-2))*costh,np.cos(np.radians(90-2))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(90-1))*costh,np.cos(np.radians(90-1))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(10))],[0,np.sin(np.radians(10))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(20))],[0,np.sin(np.radians(20))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(30))],[0,np.sin(np.radians(30))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(40))],[0,np.sin(np.radians(40))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(50))],[0,np.sin(np.radians(50))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(60))],[0,np.sin(np.radians(60))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(70))],[0,np.sin(np.radians(70))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(80))],[0,np.sin(np.radians(80))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(90))],[0,np.sin(np.radians(90))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(100))],[0,np.sin(np.radians(100))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(110))],[0,np.sin(np.radians(110))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(120))],[0,np.sin(np.radians(120))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(130))],[0,np.sin(np.radians(130))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(140))],[0,np.sin(np.radians(140))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(150))],[0,np.sin(np.radians(150))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(160))],[0,np.sin(np.radians(160))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(170))],[0,np.sin(np.radians(170))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(180))],[0,np.sin(np.radians(180))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(190))],[0,np.sin(np.radians(190))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(200))],[0,np.sin(np.radians(200))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(210))],[0,np.sin(np.radians(210))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(220))],[0,np.sin(np.radians(220))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(230))],[0,np.sin(np.radians(230))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(240))],[0,np.sin(np.radians(240))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(250))],[0,np.sin(np.radians(250))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(260))],[0,np.sin(np.radians(260))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(270))],[0,np.sin(np.radians(270))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(280))],[0,np.sin(np.radians(280))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(290))],[0,np.sin(np.radians(290))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(300))],[0,np.sin(np.radians(300))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(310))],[0,np.sin(np.radians(310))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(320))],[0,np.sin(np.radians(320))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(330))],[0,np.sin(np.radians(330))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(340))],[0,np.sin(np.radians(340))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(350))],[0,np.sin(np.radians(350))],color='gray',linestyle='-',linewidth=0.25)
ax.axhline(color='black',linestyle='-',linewidth=1)
ax.axvline(color='black',linestyle='-',linewidth=1)
ax.scatter(x_noclones,y_noclones,color='tomato',s=1) # plutino poles
Npts_circle = 100
circle_clock_angles = np.linspace(start=0,stop=2*np.pi,num=Npts_circle,endpoint=True)
q_noclones = df['q_noclones'][0]
radius_circle = np.sin(q_noclones)
hx_circle_pre = np.cos(circle_clock_angles)*radius_circle
hy_circle_pre = np.sin(circle_clock_angles)*radius_circle
hz_circle_pre = np.cos(q_noclones)*np.ones(Npts_circle)
# rotate the confidence circle around to center on the vMF mean pole
alphahat_noclones = df['alphahat_noclones'][0]
betahat_noclones = df['betahat_noclones'][0]
kappahat_noclones = df['kappahat_noclones'][0]
A_vmf_noclones = np.array([[np.cos(alphahat_noclones)*np.cos(betahat_noclones),np.cos(alphahat_noclones)*np.sin(betahat_noclones),-np.sin(alphahat_noclones)],\
                          [-np.sin(betahat_noclones),np.cos(betahat_noclones),0],\
                          [np.sin(alphahat_noclones)*np.cos(betahat_noclones),np.sin(alphahat_noclones)*np.sin(betahat_noclones),np.cos(alphahat_noclones)]])
hvec_pre = np.array([hx_circle_pre,hy_circle_pre,hz_circle_pre])
hvec_post = np.matmul(np.transpose(A_vmf_noclones),hvec_pre)
hx_circle_post = hvec_post[0,:]
hy_circle_post = hvec_post[1,:]
hz_circle_post = hvec_post[2,:]
ax.plot(hx_circle_post,hy_circle_post,color='blue',linestyle='-',linewidth=1) # confidence circle
ax.scatter(xhat_noclones,yhat_noclones,color='blue',s=20,marker='o') # vmf midplane of plutinos
A_vm17_noclones = np.array([[np.cos(theta_mid_noclones)*np.cos(phi_mid_noclones),np.cos(theta_mid_noclones)*np.sin(phi_mid_noclones),-np.sin(theta_mid_noclones)],\
                          [-np.sin(phi_mid_noclones),np.cos(phi_mid_noclones),0],\
                          [np.sin(theta_mid_noclones)*np.cos(phi_mid_noclones),np.sin(theta_mid_noclones)*np.sin(phi_mid_noclones),np.cos(theta_mid_noclones)]])
hvec_post = np.matmul(np.transpose(A_vm17_noclones),hvec_pre)
hx_circle_post = hvec_post[0,:]
hy_circle_post = hvec_post[1,:]
hz_circle_post = hvec_post[2,:]
ax.plot(hx_circle_post,hy_circle_post,color='forestgreen',linestyle='-',linewidth=1) # confidence circle
ax.scatter(xhat_vm17_noclones,yhat_vm17_noclones,color='forestgreen',s=20,marker='D') # vm17 midplane of plutinos
i_invariable = np.radians(1.578694)
W_invariable = np.radians(107.582222)
hx_invariable = np.sin(i_invariable)*np.sin(W_invariable)
hy_invariable = -np.sin(i_invariable)*np.cos(W_invariable)
hz_invariable = np.cos(i_invariable)
ax.scatter(hx_invariable,hy_invariable,color='red',s=20,marker='v') # invariable pole of the solar system
date = '20221012'
year = 2022
month = 1
day = 1
hour = 0
minute = 0
second = 0
JD = JDfun(year,month,day,hour,minute,second)
JD = str(JD)
name = '8' # neptune barycenter
center = '500@0' # solar system barycenter
obj = Horizons(id=name,location=center,epochs=JD)
el = obj.elements()
i_neptune = np.radians(float(el['incl']))
W_neptune = np.radians(float(el['Omega']))
# equation 1 in the paper
hx_neptune = np.sin(i_neptune)*np.sin(W_neptune)
hy_neptune = -np.sin(i_neptune)*np.cos(W_neptune)
ax.scatter(hx_neptune,hy_neptune,color='blue',s=20,marker='>') # orbit pole of neptune
plt.tight_layout()
plt.axis('equal')
ax.set_xlabel('$h_x$')
ax.set_ylabel('$h_y$')
ax.set_xlim([-0.07,0.09])
ax.set_ylim([-0.07,0.09])
ax.set_box_aspect(1)
titlestr = 'fig5plothxhydetail_noclones'
plt.tight_layout()
plt.savefig(titlestr + '.eps',dpi=300,bbox_inches='tight',pad_inches=0)
plt.show()
#%% FIGURE 6 area-weighted circular histogram of ecliptic longitudes only
fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
n,bins,patches=circular_hist(ax,np.array(WPb_rad_clones),edgecolor='black',bins=16)
titlestr = 'fig6cirhistogram_clones'
plt.savefig(titlestr + '.eps',dpi=300,bbox_inches='tight',pad_inches=0)
plt.show()
fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
n,bins,patches=circular_hist(ax,np.array(WPb_rad_noclones),edgecolor='black',bins=16)
titlestr = 'fig6cirhistogram_noclones'
plt.savefig(titlestr + '.eps',dpi=300,bbox_inches='tight',pad_inches=0)
plt.show()
#%% FIGURE 7 histogram of relative inclinations
A_vmf_clones = np.array([[np.cos(alphahat_clones)*np.cos(betahat_clones),np.cos(alphahat_clones)*np.sin(betahat_clones),-np.sin(alphahat_clones)],\
                          [-np.sin(betahat_clones),np.cos(betahat_clones),0],\
                          [np.sin(alphahat_clones)*np.cos(betahat_clones),np.sin(alphahat_clones)*np.sin(betahat_clones),np.cos(alphahat_clones)]])
xvec_clones = np.array([x_clones,y_clones,z_clones])
xivec_clones = np.matmul(A_vmf_clones,xvec_clones)
xi_clones = xivec_clones[0,:]
yi_clones = xivec_clones[1,:]
zi_clones = xivec_clones[2,:]
thetai_clones = np.arccos(zi_clones)
thetai_degrees_clones = np.degrees(thetai_clones)
radians_vec = np.linspace(0,np.max(thetai_clones),num=1000000,endpoint=True)
curve_vec = kappahat_clones/(1-np.exp(-2*kappahat_clones))*np.sin(radians_vec) * \
    np.exp(-2*kappahat_clones*np.sin(radians_vec/2)**2)
fig = plt.figure(figsize=(2,2))
ax1  = fig.add_subplot(111)
ax1.set_xlabel('Relative inclination (radians)')
ax1.set_xticks([0,np.pi/4,np.pi/2,3*np.pi/4,np.pi])
ax1.set_xticklabels(['0','$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$'])
plt.rcParams['font.size'] = 6
data = thetai_clones
myHist = np.histogram(data,bins=20,density=True)
myHist_binheights = myHist[0]
myHist_binboundaries = myHist[1]
myHist_binwidth = myHist_binboundaries[1]-myHist_binboundaries[0]
myHist_plot = ax1.stairs(values=myHist_binheights,edges=myHist_binboundaries,ec='blue',label='myHist')
data2 = iPb_rad_clones
myHist2 = np.histogram(data2,bins=20,density=True)
myHist2_binheights = myHist2[0]
myHist2_binboundaries = myHist2[1]
myHist2_binwidth = myHist2_binboundaries[1]-myHist2_binboundaries[0]
myHist2_plot = ax1.stairs(values=myHist2_binheights,edges=myHist2_binboundaries,ec='gray',label='myHist2')
A_vm17_clones = np.array([[np.cos(theta_mid_clones)*np.cos(phi_mid_clones),np.cos(theta_mid_clones)*np.sin(phi_mid_clones),-np.sin(theta_mid_clones)],\
                          [-np.sin(phi_mid_clones),np.cos(phi_mid_clones),0],\
                          [np.sin(theta_mid_clones)*np.cos(phi_mid_clones),np.sin(theta_mid_clones)*np.sin(phi_mid_clones),np.cos(theta_mid_clones)]])
xivec_vm17_clones = np.matmul(A_vm17_clones,xvec_clones)
xi_vm17_clones = xivec_vm17_clones[0,:]
yi_vm17_clones = xivec_vm17_clones[1,:]
zi_vm17_clones = xivec_vm17_clones[2,:]
thetai_vm17_clones = np.arccos(zi_vm17_clones)
thetai_vm17_degrees_clones = np.degrees(thetai_vm17_clones)
data3 = thetai_vm17_clones
myHist3 = np.histogram(data3,bins=20,density=True)
myHist3_binheights = myHist3[0]
myHist3_binboundaries = myHist3[1]
myHist3_binwidth = myHist3_binboundaries[1]-myHist3_binboundaries[0]
myHist3_plot = ax1.stairs(values=myHist3_binheights,edges=myHist3_binboundaries,ec='forestgreen',label='myHist3')
sigmahat_mle_clones = df['sigmahat_mle_clones'][0]
g = ax1.plot(radians_vec[0:-1:100],curve_vec[0:-1:100],color='gold',zorder=3,lw=1,label='g')
curve_vec_3 = truncated_rayleigh_pdf(radians_vec,1/sigmahat_mle_clones**2)
h = ax1.plot(radians_vec[0:-1:40000],curve_vec_3[0:-1:40000],color='red',linestyle='none',\
              marker='D',markersize=1,zorder=15,lw=1,label='h')
maxx = np.max([np.max(thetai_clones),np.max(thetai_vm17_clones)])
maxy = np.max([np.max(curve_vec),np.max(curve_vec_3),np.max(myHist[0]),\
               np.max(myHist2[0]),np.max(myHist3[0])])
maxx = 1.1*maxx
maxy = 1.1*maxy
ax1.set_xlim([0,maxx])
ax1.set_ylim([0,maxy])
import matplotlib.lines as mlines
red_dots = mlines.Line2D([], [], color='red', marker='D',
                          markersize=1, label='Eq. 50, $\hat{\sigma}_{MLE}=0.18$')
gold_line = mlines.Line2D([],[],color='gold',label='Eq. 46, $\kappa=32$')
blue_hist = mlines.Line2D([],[],color='blue',label='$i_{rel}$ (vMF)')
green_hist = mlines.Line2D([],[],color='forestgreen',label='$i_{rel}$ (debiased)')
black_hist = mlines.Line2D([],[],color='gray',label='$i$ (ecliptic)')
ax1.legend(handles=[gold_line,red_dots,black_hist,blue_hist,green_hist])
titlestr = 'fig7inclinationhistogram_clones'
plt.savefig(titlestr + '.eps',dpi=300,bbox_inches='tight',pad_inches=0)
plt.show()
alphahat_noclones = df['alphahat_noclones'][0]
betahat_noclones = df['betahat_noclones'][0]
kappahat_noclones = df['kappahat_noclones'][0]
A_vmf_noclones = np.array([[np.cos(alphahat_noclones)*np.cos(betahat_noclones),np.cos(alphahat_noclones)*np.sin(betahat_noclones),-np.sin(alphahat_noclones)],\
                          [-np.sin(betahat_noclones),np.cos(betahat_noclones),0],\
                          [np.sin(alphahat_noclones)*np.cos(betahat_noclones),np.sin(alphahat_noclones)*np.sin(betahat_noclones),np.cos(alphahat_noclones)]])
xvec_noclones = np.array([x_noclones,y_noclones,z_noclones])
xivec_noclones = np.matmul(A_vmf_noclones,xvec_noclones)
xi_noclones = xivec_noclones[0,:]
yi_noclones = xivec_noclones[1,:]
zi_noclones = xivec_noclones[2,:]
thetai_noclones = np.arccos(zi_noclones)
thetai_degrees_noclones = np.degrees(thetai_noclones)
radians_vec = np.linspace(0,np.max(thetai_noclones),num=1000000,endpoint=True)
curve_vec = kappahat_noclones/(1-np.exp(-2*kappahat_noclones))*np.sin(radians_vec) * \
    np.exp(-2*kappahat_noclones*np.sin(radians_vec/2)**2)
fig = plt.figure(figsize=(2,2))
ax1  = fig.add_subplot(111)
ax1.set_xlabel('Relative inclination (radians)')
ax1.set_xticks([0,np.pi/4,np.pi/2,3*np.pi/4,np.pi])
ax1.set_xticklabels(['0','$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$'])
plt.rcParams['font.size'] = 6
data = thetai_noclones
myHist = np.histogram(data,bins=20,density=True)
myHist_binheights = myHist[0]
myHist_binboundaries = myHist[1]
myHist_binwidth = myHist_binboundaries[1]-myHist_binboundaries[0]
myHist_plot = ax1.stairs(values=myHist_binheights,edges=myHist_binboundaries,ec='blue',label='myHist')
data2 = iPb_rad_noclones
myHist2 = np.histogram(data2,bins=20,density=True)
myHist2_binheights = myHist2[0]
myHist2_binboundaries = myHist2[1]
myHist2_binwidth = myHist2_binboundaries[1]-myHist2_binboundaries[0]
myHist2_plot = ax1.stairs(values=myHist2_binheights,edges=myHist2_binboundaries,ec='gray',label='myHist2')
A_vm17_noclones = np.array([[np.cos(theta_mid_noclones)*np.cos(phi_mid_noclones),np.cos(theta_mid_noclones)*np.sin(phi_mid_noclones),-np.sin(theta_mid_noclones)],\
                          [-np.sin(phi_mid_noclones),np.cos(phi_mid_noclones),0],\
                          [np.sin(theta_mid_noclones)*np.cos(phi_mid_noclones),np.sin(theta_mid_noclones)*np.sin(phi_mid_noclones),np.cos(theta_mid_noclones)]])
xivec_vm17_noclones = np.matmul(A_vm17_noclones,xvec_noclones)
xi_vm17_noclones = xivec_vm17_noclones[0,:]
yi_vm17_noclones = xivec_vm17_noclones[1,:]
zi_vm17_noclones = xivec_vm17_noclones[2,:]
thetai_vm17_noclones = np.arccos(zi_vm17_noclones)
thetai_vm17_degrees_noclones = np.degrees(thetai_vm17_noclones)
data3 = thetai_vm17_noclones
myHist3 = np.histogram(data3,bins=20,density=True)
myHist3_binheights = myHist3[0]
myHist3_binboundaries = myHist3[1]
myHist3_binwidth = myHist3_binboundaries[1]-myHist3_binboundaries[0]
myHist3_plot = ax1.stairs(values=myHist3_binheights,edges=myHist3_binboundaries,ec='forestgreen',label='myHist3')
sigmahat_mle_noclones = df['sigmahat_mle_noclones'][0]
g = ax1.plot(radians_vec[0:-1:100],curve_vec[0:-1:100],color='gold',zorder=3,lw=1,label='g')
curve_vec_3 = truncated_rayleigh_pdf(radians_vec,1/sigmahat_mle_noclones**2)
h = ax1.plot(radians_vec[0:-1:40000],curve_vec_3[0:-1:40000],color='red',linestyle='none',\
              marker='D',markersize=1,zorder=15,lw=1,label='h')
maxx = np.max([np.max(thetai_noclones),np.max(thetai_vm17_noclones)])
maxy = np.max([np.max(curve_vec),np.max(curve_vec_3),np.max(myHist[0]),\
               np.max(myHist2[0]),np.max(myHist3[0])])
maxx = 1.1*maxx
maxy = 1.1*maxy
ax1.set_xlim([0,maxx])
ax1.set_ylim([0,maxy])
import matplotlib.lines as mlines
red_dots = mlines.Line2D([], [], color='red', marker='D',
                          markersize=1, label='Eq. 50, $\hat{\sigma}_{MLE}=0.18$')
gold_line = mlines.Line2D([],[],color='gold',label='Eq. 46, $\kappa=32$')
blue_hist = mlines.Line2D([],[],color='blue',label='$i_{rel}$ (vMF)')
green_hist = mlines.Line2D([],[],color='forestgreen',label='$i_{rel}$ (debiased)')
black_hist = mlines.Line2D([],[],color='gray',label='$i$ (ecliptic)')
ax1.legend(handles=[gold_line,red_dots,black_hist,blue_hist,green_hist])
titlestr = 'fig7inclinationhistogram_noclones'
plt.savefig(titlestr + '.eps',dpi=300,bbox_inches='tight',pad_inches=0)
plt.show()
#%% some more numbers
d1 = np.degrees(np.arccos(np.dot([hx_invariable,hy_invariable,hz_invariable],\
                                 [xhat_noclones,yhat_noclones,zhat_noclones])))
d2 = np.degrees(np.arccos(np.dot([hx_invariable,hy_invariable,hz_invariable],\
                                 [xhat_vm17_noclones,yhat_vm17_noclones,zhat_vm17_noclones])))
d3 = np.degrees(np.arccos(np.dot([hx_neptune,hy_neptune,hz_neptune],\
                                 [xhat_noclones,yhat_noclones,zhat_noclones])))
d4 = np.degrees(np.arccos(np.dot([hx_neptune,hy_neptune,hz_neptune],\
                                 [xhat_vm17_noclones,yhat_vm17_noclones,zhat_vm17_noclones])))
d5 = np.degrees(np.arccos(np.dot([xhat_vm17_noclones,yhat_vm17_noclones,zhat_vm17_noclones],\
                                 [xhat_noclones,yhat_noclones,zhat_noclones])))
