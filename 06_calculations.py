#%%
def dl_dsigma(sigma,uvec):
    n = len(uvec)
    term1 = -2*n*sigma**2
    term2 = np.sum(uvec**2)
    term3coefficient = np.exp(-np.pi**2/2/sigma**2)
    term3parenthesis = term1 - term2 + n*np.pi**2
    val = term1 + term2 + term3coefficient*term3parenthesis
    return val
#%%
import pandas as pd
import numpy as np
from scipy import stats
from scipy.optimize import fsolve
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
# section 2.2 testing for uniformity
Sx_clones = np.sum(x_clones)
Sy_clones = np.sum(y_clones)
Sz_clones = np.sum(z_clones)
R_clones = np.sqrt(Sx_clones**2 + Sy_clones**2 + Sz_clones**2)
xhat_clones = Sx_clones/R_clones
yhat_clones = Sy_clones/R_clones
zhat_clones = Sz_clones/R_clones
ihat_clones = np.arccos(zhat_clones)
ihat_degrees_clones = np.degrees(ihat_clones)
What_clones = np.arctan2(xhat_clones,-yhat_clones)
What_degrees_clones = np.degrees(What_clones)
threeRsquaredovern_clones = 3*R_clones**2 / n_clones
# section 2.3 testing for rotational symmetry
T11_clones = np.sum(x_clones**2)
T12_clones = np.sum(x_clones*y_clones)
T13_clones = np.sum(x_clones*z_clones)
T21_clones = T12_clones
T22_clones = np.sum(y_clones**2)
T23_clones = np.sum(y_clones*z_clones)
T31_clones = T13_clones
T32_clones = T23_clones
T33_clones = np.sum(z_clones**2)
T_clones = np.array([[T11_clones,T12_clones,T13_clones],\
                     [T21_clones,T22_clones,T23_clones],\
                     [T31_clones,T32_clones,T33_clones]])
eigs_clones = np.linalg.eig(T_clones)
eigvals_clones = eigs_clones[0]
eigvecs_clones = eigs_clones[1]
tau1hat_clones = np.min(eigvals_clones)
tau3hat_clones = np.max(eigvals_clones)
tau2hat_clones = eigvals_clones[eigvals_clones!=tau1hat_clones]
tau2hat_clones = tau2hat_clones[tau2hat_clones!=tau3hat_clones]
tau2hat_clones = tau2hat_clones[0]
tau1bar_clones = tau1hat_clones/np.sum(eigvals_clones)
tau2bar_clones = tau2hat_clones/np.sum(eigvals_clones)
tau3bar_clones = tau3hat_clones/np.sum(eigvals_clones)
if tau3hat_clones == eigvals_clones[0]:
    u3hat_clones = eigvecs_clones[0]
if tau3hat_clones == eigvals_clones[1]:
    u3hat_clones = eigvecs_clones[1]
if tau3hat_clones == eigvals_clones[2]:
    u3hat_clones = eigvecs_clones[2]
if u3hat_clones[0] < 0:
    u3hat_clones = -u3hat_clones
lambdahat_clones = u3hat_clones[0]
muhat_clones = u3hat_clones[1]
nuhat_clones = u3hat_clones[2]
Gamma_clones = 1/n_clones * np.sum( (lambdahat_clones*x_clones+\
                    muhat_clones*y_clones+nuhat_clones*z_clones)**4 )
Pn_clones = 2 * n_clones * (tau2bar_clones-tau1bar_clones)**2 / (1-2*tau3bar_clones+Gamma_clones)
p_clones = np.exp(-1/2 * Pn_clones)
# section 2.4 estimating the mean direction
thetahat_clones = np.arccos(zhat_clones)
phihat_clones = np.arctan2(yhat_clones,xhat_clones)
alphahat_clones = thetahat_clones
betahat_clones = phihat_clones
alphahat_degrees_clones = np.degrees(alphahat_clones)
betahat_degrees_clones = np.degrees(betahat_clones)
# section 2.5 estimating a confidence region for the mean direction
Rbar_clones = R_clones/n_clones
d_clones = 1 - 1/n_clones * np.sum( (x_clones*xhat_clones + \
                    y_clones*yhat_clones + z_clones*zhat_clones)**2 )
sigmahat_clones = np.sqrt(d_clones/(n_clones*Rbar_clones**2))
alphaconfidence = 0.003 # 99.7 per cent confidence region
q_clones = np.arcsin(sigmahat_clones*np.sqrt(-np.log(alphaconfidence)))
q_degrees_clones = np.degrees(q_clones)
# section 2.6 estimating the concentration parameter
kappahat_clones = (n_clones-1)/(n_clones-R_clones)
dof_clones = 2*n_clones - 2
kappaL_clones = 1/2 * stats.chi2.isf(1-alphaconfidence/2,dof_clones) / (n_clones-R_clones)
kappaU_clones = 1/2 * stats.chi2.isf(alphaconfidence/2,dof_clones) / (n_clones-R_clones)
kappaL_16_clones = 1/2 * stats.chi2.isf(1-0.32/2,dof_clones) / (n_clones-R_clones)
kappaU_84_clones = 1/2 * stats.chi2.isf(0.32/2,dof_clones) / (n_clones-R_clones)
qkappa_clones = np.arccos(1-(n_clones-R_clones)/R_clones*(alphaconfidence**(-1/(n_clones-1))-1))
qkappa_degrees_clones = np.degrees(qkappa_clones)
# section 2.7 checking goodness of fit
A_vmf_clones = np.array([[np.cos(alphahat_clones)*np.cos(betahat_clones),np.cos(alphahat_clones)*np.sin(betahat_clones),-np.sin(alphahat_clones)],\
                         [-np.sin(betahat_clones),np.cos(betahat_clones),0],\
                         [np.sin(alphahat_clones)*np.cos(betahat_clones),np.sin(alphahat_clones)*np.sin(betahat_clones),np.cos(alphahat_clones)]])
xvec_clones = np.array([x_clones,y_clones,z_clones])
xivec_clones = np.matmul(A_vmf_clones,xvec_clones)
xi_clones = xivec_clones[0,:]
yi_clones = xivec_clones[1,:]
zi_clones = xivec_clones[2,:]
thetai_clones = np.arccos(zi_clones)
phii_clones = np.arctan2(yi_clones,xi_clones)
# section 2.7.1 colatitudes
Xi_clones = np.sort(1-np.cos(thetai_clones))
kappahat_colatitudes_clones = (n_clones-1)/np.sum(Xi_clones)
FXi_clones = 1 - np.exp(-kappahat_colatitudes_clones*Xi_clones)
Dnminus_clones = np.max(FXi_clones-(np.linspace(start=1,stop=n_clones,num=n_clones,endpoint=True)-1)/n_clones)
Dnplus_clones = np.max(np.linspace(start=1,stop=n_clones,num=n_clones,endpoint=True)/n_clones-FXi_clones)
Dn_clones = np.max(np.array([Dnminus_clones,Dnplus_clones]))
ME_clones = (Dn_clones-0.2/np.sqrt(n_clones))*(np.sqrt(n_clones)+0.26+0.5/np.sqrt(n_clones))
# section 2.7.2 longitudes
Xi_clones = np.sort(phii_clones/(2*np.pi))
FXi_clones = Xi_clones
Dnminus_clones = np.max(FXi_clones-(np.linspace(start=1,stop=n_clones,num=n_clones,endpoint=True)-1)/n_clones)
Dnplus_clones = np.max(np.linspace(start=1,stop=n_clones,num=n_clones,endpoint=True)/n_clones-FXi_clones)
Vn_clones = Dnminus_clones + Dnplus_clones
MU_clones = Vn_clones * (np.sqrt(n_clones) - 0.467 + 1.623/np.sqrt(n_clones))
# section 2.7.3 two-variable test
A_vmf_clones = np.array([[np.cos(3*np.pi/2-alphahat_clones)*np.cos(betahat_clones-np.pi),np.cos(3*np.pi/2-alphahat_clones)*np.sin(betahat_clones-np.pi),-np.sin(3*np.pi/2-alphahat_clones)],\
                         [-np.sin(betahat_clones-np.pi),np.cos(betahat_clones-np.pi),0],\
                         [np.sin(3*np.pi/2-alphahat_clones)*np.cos(betahat_clones-np.pi),np.sin(3*np.pi/2-alphahat_clones)*np.sin(betahat_clones-np.pi),np.cos(3*np.pi/2-alphahat_clones)]])
xiivec_clones = np.matmul(A_vmf_clones,xvec_clones)
xii_clones = xiivec_clones[0,:]
yii_clones = xiivec_clones[1,:]
zii_clones = xiivec_clones[2,:]
thetaii_clones = np.arccos(zii_clones)
phiii_clones = np.arctan2(yii_clones,xii_clones)
Xi_clones = np.sort( (phiii_clones-np.pi)*np.sqrt(np.sin(thetaii_clones)) )
s = np.sqrt( 1/n_clones * np.sum(Xi_clones**2) )
xj_clones = np.sort(Xi_clones/s)
Fxj_clones = stats.norm.cdf(xj_clones)
Dnminus_clones = np.max(Fxj_clones-(np.linspace(start=1,stop=n_clones,num=n_clones,endpoint=True)-1)/n_clones)
Dnplus_clones = np.max(np.linspace(start=1,stop=n_clones,num=n_clones,endpoint=True)/n_clones-Fxj_clones)
Dn_clones = np.max(np.array([Dnminus_clones,Dnplus_clones]))
MN_clones = Dn_clones * ( np.sqrt(n_clones) - 0.01 + 0.85/np.sqrt(n_clones) )
# section 2.8 relationship of the vMF distribution to the univariate Rayleigh distribution
sigma_clones = 1/np.sqrt(kappahat_clones)
sigma_degrees_clones = np.degrees(sigma_clones)
CR_clones = 1/( 1-np.exp(-np.pi**2/(2*sigma_clones**2)) )
uvec = thetai_clones
func = lambda sigma : dl_dsigma(sigma,uvec)
sigmahat_mle_clones = fsolve(func,sigma_clones)
sigmahat_mle_degrees_clones = np.degrees(sigmahat_mle_clones[0])
# section 3 mitigating observational bias
EPb_rad_clones = [] # eccentric anomaly
fPb_rad_clones = [] # true anomaly
Etol = 1e-13
ratio = 1
for iobj in range(n_clones): # Curtis algorithm 3.1, really just Newton-Raphson
    M = MPb_rad_clones[iobj]
    ecc = ePb_clones[iobj]
    if M < np.pi:
        E = M + ecc/2
    else:
        E = M - ecc/2
    while np.abs(ratio) > Etol:
        top = E - ecc*np.sin(E) - M
        bottom = 1 - ecc*np.cos(E)
        ratio = top/bottom
        E = E - ratio
    f = 2 * np.arctan(np.sqrt(1+ecc)/np.sqrt(1-ecc)*np.tan(E/2)) # Curtis eq 3.13a
    EPb_rad_clones.append(E)
    fPb_rad_clones.append(f)
fPb_rad_clones = np.array(fPb_rad_clones)
thetaPb_rad_clones = wPb_rad_clones + fPb_rad_clones
hx_clones = x_clones # orbit normal vectors (unit angular momentum vectors)
hy_clones = y_clones
hz_clones = z_clones
x_clones = np.cos(WPb_rad_clones)*np.cos(thetaPb_rad_clones) - \
    np.sin(WPb_rad_clones)*np.sin(thetaPb_rad_clones)*np.cos(iPb_rad_clones) # unit position vectors, Schaub & Junkins eq 9.164
y_clones = np.sin(WPb_rad_clones)*np.cos(thetaPb_rad_clones) + \
    np.cos(WPb_rad_clones)*np.sin(thetaPb_rad_clones)*np.cos(iPb_rad_clones)
z_clones = np.sin(thetaPb_rad_clones)*np.sin(iPb_rad_clones)
vtx_clones = []
vty_clones = []
vtz_clones = []
for iobj in range(n_clones):
    x = x_clones[iobj]
    y = y_clones[iobj]
    z = z_clones[iobj]
    xvec = np.array([x,y,z])
    hx = hx_clones[iobj]
    hy = hy_clones[iobj]
    hz = hz_clones[iobj]
    hvec = np.array([hx,hy,hz])
    vtvec = np.cross(hvec,xvec)
    vtx_clones.append(vtvec[0])
    vty_clones.append(vtvec[1])
    vtz_clones.append(vtvec[2])
vtx_clones = np.array(vtx_clones)
vty_clones = np.array(vty_clones)
vtz_clones = np.array(vtz_clones)
dictionary_vtvec = {'packed_designation':des_clones,'vtx':vtx_clones,'vty':vty_clones,'vtz':vtz_clones}
df_vtvec = pd.DataFrame.from_dict(dictionary_vtvec)
df_vtvec.to_csv('vtvec_clones.csv',index=False)

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
# section 2.2 testing for uniformity
Sx_noclones = np.sum(x_noclones)
Sy_noclones = np.sum(y_noclones)
Sz_noclones = np.sum(z_noclones)
R_noclones = np.sqrt(Sx_noclones**2 + Sy_noclones**2 + Sz_noclones**2)
xhat_noclones = Sx_noclones/R_noclones
yhat_noclones = Sy_noclones/R_noclones
zhat_noclones = Sz_noclones/R_noclones
ihat_noclones = np.arccos(zhat_noclones)
ihat_degrees_noclones = np.degrees(ihat_noclones)
What_noclones = np.arctan2(xhat_noclones,-yhat_noclones)
What_degrees_noclones = np.degrees(What_noclones)
threeRsquaredovern_noclones = 3*R_noclones**2 / n_noclones
# section 2.3 testing for rotational symmetry
T11_noclones = np.sum(x_noclones**2)
T12_noclones = np.sum(x_noclones*y_noclones)
T13_noclones = np.sum(x_noclones*z_noclones)
T21_noclones = T12_noclones
T22_noclones = np.sum(y_noclones**2)
T23_noclones = np.sum(y_noclones*z_noclones)
T31_noclones = T13_noclones
T32_noclones = T23_noclones
T33_noclones = np.sum(z_noclones**2)
T_noclones = np.array([[T11_noclones,T12_noclones,T13_noclones],\
                     [T21_noclones,T22_noclones,T23_noclones],\
                     [T31_noclones,T32_noclones,T33_noclones]])
eigs_noclones = np.linalg.eig(T_noclones)
eigvals_noclones = eigs_noclones[0]
eigvecs_noclones = eigs_noclones[1]
tau1hat_noclones = np.min(eigvals_noclones)
tau3hat_noclones = np.max(eigvals_noclones)
tau2hat_noclones = eigvals_noclones[eigvals_noclones!=tau1hat_noclones]
tau2hat_noclones = tau2hat_noclones[tau2hat_noclones!=tau3hat_noclones]
tau2hat_noclones = tau2hat_noclones[0]
tau1bar_noclones = tau1hat_noclones/np.sum(eigvals_noclones)
tau2bar_noclones = tau2hat_noclones/np.sum(eigvals_noclones)
tau3bar_noclones = tau3hat_noclones/np.sum(eigvals_noclones)
if tau3hat_noclones == eigvals_noclones[0]:
    u3hat_noclones = eigvecs_noclones[0]
if tau3hat_noclones == eigvals_noclones[1]:
    u3hat_noclones = eigvecs_noclones[1]
if tau3hat_noclones == eigvals_noclones[2]:
    u3hat_noclones = eigvecs_noclones[2]
if u3hat_noclones[0] < 0:
    u3hat_noclones = -u3hat_noclones
lambdahat_noclones = u3hat_noclones[0]
muhat_noclones = u3hat_noclones[1]
nuhat_noclones = u3hat_noclones[2]
Gamma_noclones = 1/n_noclones * np.sum( (lambdahat_noclones*x_noclones+\
                    muhat_noclones*y_noclones+nuhat_noclones*z_noclones)**4 )
Pn_noclones = 2 * n_noclones * (tau2bar_noclones-tau1bar_noclones)**2 / (1-2*tau3bar_noclones+Gamma_noclones)
p_noclones = np.exp(-1/2 * Pn_noclones)
# section 2.4 estimating the mean direction
thetahat_noclones = np.arccos(zhat_noclones)
phihat_noclones = np.arctan2(yhat_noclones,xhat_noclones)
alphahat_noclones = thetahat_noclones
betahat_noclones = phihat_noclones
alphahat_degrees_noclones = np.degrees(alphahat_noclones)
betahat_degrees_noclones = np.degrees(betahat_noclones)
# section 2.5 estimating a confidence region for the mean direction
Rbar_noclones = R_noclones/n_noclones
d_noclones = 1 - 1/n_noclones * np.sum( (x_noclones*xhat_noclones + \
                    y_noclones*yhat_noclones + z_noclones*zhat_noclones)**2 )
sigmahat_noclones = np.sqrt(d_noclones/(n_noclones*Rbar_noclones**2))
alphaconfidence = 0.003 # 99.7 per cent confidence region
q_noclones = np.arcsin(sigmahat_noclones*np.sqrt(-np.log(alphaconfidence)))
q_degrees_noclones = np.degrees(q_noclones)
# section 2.6 estimating the concentration parameter
kappahat_noclones = (n_noclones-1)/(n_noclones-R_noclones)
dof_noclones = 2*n_noclones - 2
kappaL_noclones = 1/2 * stats.chi2.isf(1-alphaconfidence/2,dof_noclones) / (n_noclones-R_noclones)
kappaU_noclones = 1/2 * stats.chi2.isf(alphaconfidence/2,dof_noclones) / (n_noclones-R_noclones)
kappaL_16_noclones = 1/2 * stats.chi2.isf(1-0.32/2,dof_noclones) / (n_noclones-R_noclones)
kappaU_84_noclones = 1/2 * stats.chi2.isf(0.32/2,dof_noclones) / (n_noclones-R_noclones)
qkappa_noclones = np.arccos(1-(n_noclones-R_noclones)/R_noclones*(alphaconfidence**(-1/(n_noclones-1))-1))
qkappa_degrees_noclones = np.degrees(qkappa_noclones)
# section 2.7 checking goodness of fit
A_vmf_noclones = np.array([[np.cos(alphahat_noclones)*np.cos(betahat_noclones),np.cos(alphahat_noclones)*np.sin(betahat_noclones),-np.sin(alphahat_noclones)],\
                         [-np.sin(betahat_noclones),np.cos(betahat_noclones),0],\
                         [np.sin(alphahat_noclones)*np.cos(betahat_noclones),np.sin(alphahat_noclones)*np.sin(betahat_noclones),np.cos(alphahat_noclones)]])
xvec_noclones = np.array([x_noclones,y_noclones,z_noclones])
xivec_noclones = np.matmul(A_vmf_noclones,xvec_noclones)
xi_noclones = xivec_noclones[0,:]
yi_noclones = xivec_noclones[1,:]
zi_noclones = xivec_noclones[2,:]
thetai_noclones = np.arccos(zi_noclones)
phii_noclones = np.arctan2(yi_noclones,xi_noclones)
# section 2.7.1 colatitudes
Xi_noclones = np.sort(1-np.cos(thetai_noclones))
kappahat_colatitudes_noclones = (n_noclones-1)/np.sum(Xi_noclones)
FXi_noclones = 1 - np.exp(-kappahat_colatitudes_noclones*Xi_noclones)
Dnminus_noclones = np.max(FXi_noclones-(np.linspace(start=1,stop=n_noclones,num=n_noclones,endpoint=True)-1)/n_noclones)
Dnplus_noclones = np.max(np.linspace(start=1,stop=n_noclones,num=n_noclones,endpoint=True)/n_noclones-FXi_noclones)
Dn_noclones = np.max(np.array([Dnminus_noclones,Dnplus_noclones]))
ME_noclones = (Dn_noclones-0.2/np.sqrt(n_noclones))*(np.sqrt(n_noclones)+0.26+0.5/np.sqrt(n_noclones))
# section 2.7.2 longitudes
Xi_noclones = np.sort(phii_noclones/(2*np.pi))
FXi_noclones = Xi_noclones
Dnminus_noclones = np.max(FXi_noclones-(np.linspace(start=1,stop=n_noclones,num=n_noclones,endpoint=True)-1)/n_noclones)
Dnplus_noclones = np.max(np.linspace(start=1,stop=n_noclones,num=n_noclones,endpoint=True)/n_noclones-FXi_noclones)
Vn_noclones = Dnminus_noclones + Dnplus_noclones
MU_noclones = Vn_noclones * (np.sqrt(n_noclones) - 0.467 + 1.623/np.sqrt(n_noclones))
# section 2.7.3 two-variable test
A_vmf_noclones = np.array([[np.cos(3*np.pi/2-alphahat_noclones)*np.cos(betahat_noclones-np.pi),np.cos(3*np.pi/2-alphahat_noclones)*np.sin(betahat_noclones-np.pi),-np.sin(3*np.pi/2-alphahat_noclones)],\
                         [-np.sin(betahat_noclones-np.pi),np.cos(betahat_noclones-np.pi),0],\
                         [np.sin(3*np.pi/2-alphahat_noclones)*np.cos(betahat_noclones-np.pi),np.sin(3*np.pi/2-alphahat_noclones)*np.sin(betahat_noclones-np.pi),np.cos(3*np.pi/2-alphahat_noclones)]])
xiivec_noclones = np.matmul(A_vmf_noclones,xvec_noclones)
xii_noclones = xiivec_noclones[0,:]
yii_noclones = xiivec_noclones[1,:]
zii_noclones = xiivec_noclones[2,:]
thetaii_noclones = np.arccos(zii_noclones)
phiii_noclones = np.arctan2(yii_noclones,xii_noclones)
Xi_noclones = np.sort( (phiii_noclones-np.pi)*np.sqrt(np.sin(thetaii_noclones)) )
s = np.sqrt( 1/n_noclones * np.sum(Xi_noclones**2) )
xj_noclones = np.sort(Xi_noclones/s)
Fxj_noclones = stats.norm.cdf(xj_noclones)
Dnminus_noclones = np.max(Fxj_noclones-(np.linspace(start=1,stop=n_noclones,num=n_noclones,endpoint=True)-1)/n_noclones)
Dnplus_noclones = np.max(np.linspace(start=1,stop=n_noclones,num=n_noclones,endpoint=True)/n_noclones-Fxj_noclones)
Dn_noclones = np.max(np.array([Dnminus_noclones,Dnplus_noclones]))
MN_noclones = Dn_noclones * ( np.sqrt(n_noclones) - 0.01 + 0.85/np.sqrt(n_noclones) )
# section 2.8 relationship of the vMF distribution to the univariate Rayleigh distribution
sigma_noclones = 1/np.sqrt(kappahat_noclones)
sigma_degrees_noclones = np.degrees(sigma_noclones)
CR_noclones = 1/( 1-np.exp(-np.pi**2/(2*sigma_noclones**2)) )
uvec = thetai_noclones
func = lambda sigma : dl_dsigma(sigma,uvec)
sigmahat_mle_noclones = fsolve(func,sigma_noclones)
sigmahat_mle_degrees_noclones = np.degrees(sigmahat_mle_noclones[0])
# section 3 mitigating observational bias
EPb_rad_noclones = [] # eccentric anomaly
fPb_rad_noclones = [] # true anomaly
Etol = 1e-13
ratio = 1
for iobj in range(n_noclones): # Curtis algorithm 3.1, really just Newton-Raphson
    M = MPb_rad_noclones[iobj]
    ecc = ePb_noclones[iobj]
    if M < np.pi:
        E = M + ecc/2
    else:
        E = M - ecc/2
    while np.abs(ratio) > Etol:
        top = E - ecc*np.sin(E) - M
        bottom = 1 - ecc*np.cos(E)
        ratio = top/bottom
        E = E - ratio
    f = 2 * np.arctan(np.sqrt(1+ecc)/np.sqrt(1-ecc)*np.tan(E/2)) # Curtis eq 3.13a
    EPb_rad_noclones.append(E)
    fPb_rad_noclones.append(f)
fPb_rad_noclones = np.array(fPb_rad_noclones)
thetaPb_rad_noclones = wPb_rad_noclones + fPb_rad_noclones
hx_noclones = x_noclones # orbit normal vectors (unit angular momentum vectors)
hy_noclones = y_noclones
hz_noclones = z_noclones
x_noclones = np.cos(WPb_rad_noclones)*np.cos(thetaPb_rad_noclones) - \
    np.sin(WPb_rad_noclones)*np.sin(thetaPb_rad_noclones)*np.cos(iPb_rad_noclones) # unit position vectors, Schaub & Junkins eq 9.164
y_noclones = np.sin(WPb_rad_noclones)*np.cos(thetaPb_rad_noclones) + \
    np.cos(WPb_rad_noclones)*np.sin(thetaPb_rad_noclones)*np.cos(iPb_rad_noclones)
z_noclones = np.sin(thetaPb_rad_noclones)*np.sin(iPb_rad_noclones)
vtx_noclones = []
vty_noclones = []
vtz_noclones = []
for iobj in range(n_noclones):
    x = x_noclones[iobj]
    y = y_noclones[iobj]
    z = z_noclones[iobj]
    xvec = np.array([x,y,z])
    hx = hx_noclones[iobj]
    hy = hy_noclones[iobj]
    hz = hz_noclones[iobj]
    hvec = np.array([hx,hy,hz])
    vtvec = np.cross(hvec,xvec)
    vtx_noclones.append(vtvec[0])
    vty_noclones.append(vtvec[1])
    vtz_noclones.append(vtvec[2])
vtx_noclones = np.array(vtx_noclones)
vty_noclones = np.array(vty_noclones)
vtz_noclones = np.array(vtz_noclones)
dictionary_vtvec = {'packed_designation':des_noclones,'vtx':vtx_noclones,'vty':vty_noclones,'vtz':vtz_noclones}
df_vtvec = pd.DataFrame.from_dict(dictionary_vtvec)
df_vtvec.to_csv('vtvec_noclones.csv',index=False)
dictionary = {\
  'n_clones':n_clones,'xhat_clones':xhat_clones,'yhat_clones':yhat_clones,'R_clones':R_clones,\
  'zhat_clones':zhat_clones,'ihat_clones':ihat_clones,'ihat_degrees_clones':ihat_degrees_clones,\
  'What_clones':What_clones,'What_degrees_clones':What_degrees_clones,'threeRsquaredovern_clones':\
  threeRsquaredovern_clones,'T11_clones':T11_clones,'T12_clones':T12_clones,'T13_clones':T13_clones,\
  'T21_clones':T21_clones,'T22_clones':T22_clones,'T23_clones':T23_clones,'T31_clones':T31_clones,\
  'T32_clones':T32_clones,'T33_clones':T33_clones,'tau1hat_clones':tau1hat_clones,\
  'tau2hat_clones':tau2hat_clones,'tau3hat_clones':tau3hat_clones,'tau1bar_clones':tau1bar_clones,\
  'tau2bar_clones':tau2bar_clones,'tau3bar_clones':tau3bar_clones,'lambdahat_clones':lambdahat_clones,\
  'muhat_clones':muhat_clones,'nuhat_clones':nuhat_clones,'Gamma_clones':Gamma_clones,\
  'Pn_clones':Pn_clones,'p_clones':p_clones,'thetahat_clones':thetahat_clones,'phihat_clones':phihat_clones,\
  'alphahat_clones':alphahat_clones,'betahat_clones':betahat_clones,'alphahat_degrees_clones':\
  alphahat_degrees_clones,'betahat_degrees_clones':betahat_degrees_clones,'Rbar_clones':Rbar_clones,\
  'd_clones':d_clones,'sigmahat_clones':sigmahat_clones,'alphaconfidence':alphaconfidence,\
  'q_clones':q_clones,'q_degrees_clones':q_degrees_clones,'kappahat_clones':kappahat_clones,\
  'dof_clones':dof_clones,'kappaL_clones':kappaL_clones,'kappaU_clones':kappaU_clones,\
  'qkappa_clones':qkappa_clones,'qkappa_degrees_clones':qkappa_degrees_clones,'ME_clones':ME_clones,\
  'MU_clones':MU_clones,'MN_clones':MN_clones,'sigma_clones':sigma_clones,'sigma_degrees_clones':\
  sigma_degrees_clones,'CR_clones':CR_clones,'sigmahat_mle_clones':sigmahat_mle_clones,\
  'sigmahat_mle_degrees_clones':sigmahat_mle_degrees_clones,\
  'n_noclones':n_noclones,'xhat_noclones':xhat_noclones,'yhat_noclones':yhat_noclones,'R_noclones':R_noclones,\
  'zhat_noclones':zhat_noclones,'ihat_noclones':ihat_noclones,'ihat_degrees_noclones':ihat_degrees_noclones,\
  'What_noclones':What_noclones,'What_degrees_noclones':What_degrees_noclones,'threeRsquaredovern_noclones':\
  threeRsquaredovern_noclones,'T11_noclones':T11_noclones,'T12_noclones':T12_noclones,'T13_noclones':T13_noclones,\
  'T21_noclones':T21_noclones,'T22_noclones':T22_noclones,'T23_noclones':T23_noclones,'T31_noclones':T31_noclones,\
  'T32_noclones':T32_noclones,'T33_noclones':T33_noclones,'tau1hat_noclones':tau1hat_noclones,\
  'tau2hat_noclones':tau2hat_noclones,'tau3hat_noclones':tau3hat_noclones,'tau1bar_noclones':tau1bar_noclones,\
  'tau2bar_noclones':tau2bar_noclones,'tau3bar_noclones':tau3bar_noclones,'lambdahat_noclones':lambdahat_noclones,\
  'muhat_noclones':muhat_noclones,'nuhat_noclones':nuhat_noclones,'Gamma_noclones':Gamma_noclones,\
  'Pn_noclones':Pn_noclones,'p_noclones':p_noclones,'thetahat_noclones':thetahat_noclones,'phihat_noclones':phihat_noclones,\
  'alphahat_noclones':alphahat_noclones,'betahat_noclones':betahat_noclones,'alphahat_degrees_noclones':\
  alphahat_degrees_noclones,'betahat_degrees_noclones':betahat_degrees_noclones,'Rbar_noclones':Rbar_noclones,\
  'd_noclones':d_noclones,'sigmahat_noclones':sigmahat_noclones,'alphaconfidence':alphaconfidence,\
  'q_noclones':q_noclones,'q_degrees_noclones':q_degrees_noclones,'kappahat_noclones':kappahat_noclones,\
  'dof_noclones':dof_noclones,'kappaL_noclones':kappaL_noclones,'kappaU_noclones':kappaU_noclones,\
  'qkappa_noclones':qkappa_noclones,'qkappa_degrees_noclones':qkappa_degrees_noclones,'ME_noclones':ME_noclones,\
  'MU_noclones':MU_noclones,'MN_noclones':MN_noclones,'sigma_noclones':sigma_noclones,'sigma_degrees_noclones':\
  sigma_degrees_noclones,'CR_noclones':CR_noclones,'sigmahat_mle_noclones':sigmahat_mle_noclones,\
  'sigmahat_mle_degrees_noclones':sigmahat_mle_degrees_noclones,\
  }
df = pd.DataFrame.from_dict(dictionary)
df.to_csv('calculations_dictionary.csv',index=False)
