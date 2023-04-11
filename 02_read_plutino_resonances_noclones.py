#%%
def datestr(datevec):
    # datevec should be [yyyy,mm,dd,hh,mm,ss]
    datestr = ''
    for i in range(len(datevec)-1):
        thisstr = str(datevec[i])
        if len(thisstr) == 1:
            thisstr = '0' + thisstr
        if i == 0:
            datestr = datestr + thisstr
        else:
            datestr = datestr + '_' + thisstr
    thisstr = str(datevec[-1])
    if datevec[-1] < 10:
        thisstr = '0' + thisstr
    datestr = datestr + '_' + thisstr
    return datestr
#%%
def JDfun(year,month,day,hour,minute,second):
    # from Vallado pg 183, valid for yrs 1900 to 2100
    JD = 367*year - np.floor(7/4*(year+np.floor(1/12*(month+9)))) + \
        np.floor(275*month/9) + day + 1721013.5 + 1/24*(hour+1/60*(second/60+minute))
    return JD
#%% K04K19H
# K14Se4J
import numpy as np
import pandas as pd
import os
date = '20221012'
year = 2022
month = 1
day = 1
hour = 0
minute = 0
second = 0
JD = JDfun(year,month,day,hour,minute,second)
datestr = datestr([year,month,day,hour,minute,second])
horizons_file = 'horizons_barycentric_nom_noclones_' + date + '_' + datestr + '.csv'
resfile = 'resonant_barycentric_nom_noclones' + date + '_' + datestr + '.csv'
path = os.getcwd()
df = pd.read_csv(horizons_file)
# packed_designation,a_au,e,i_deg,w_deg,W_deg,M_deg
Nobj = df.shape[0]
resdes = []
res_a = []
res_e = []
res_i = []
res_w = []
res_W = []
res_M = []
false_negative_list = ['K04K19H'] # objects that look resonant in 3:2 but are tagged False
false_positive_list = ['K14Se4J','l0308'] # objects that don't look resonant in 3:2 but are tagged True
for iobj in range(Nobj):
    des = df['packed_designation'][iobj]
    if des in false_negative_list:
        resdes.append(des)
        res_a.append(df['a_au'][iobj])
        res_e.append(df['e'][iobj])
        res_i.append(df['i_deg'][iobj])
        res_w.append(df['w_deg'][iobj])
        res_W.append(df['W_deg'][iobj])
        res_M.append(df['M_deg'][iobj])
    elif des not in false_positive_list:
        for i in os.listdir(path):
            if os.path.isfile(os.path.join(path,i)) and ('long_True_plutino ' + des) in i:
                resdes.append(des)
                res_a.append(df['a_au'][iobj])
                res_e.append(df['e'][iobj])
                res_i.append(df['i_deg'][iobj])
                res_w.append(df['w_deg'][iobj])
                res_W.append(df['W_deg'][iobj])
                res_M.append(df['M_deg'][iobj])
Nres = len(resdes)
print('number of source objects =',Nobj,'number of plutinos =',Nres)
dictionary = {'packed_designation':resdes,'a_au':res_a,'e':res_e,'i_deg':res_i,\
              'w_deg':res_w,'W_deg':res_W,'M_deg':res_M}
df = pd.DataFrame.from_dict(dictionary)
df.to_csv(resfile,index=False)
df = pd.read_csv(resfile)
des_list = df['packed_designation'].tolist()
a_list = df['a_au'].to_list()
e_list = df['e'].to_list()
i_list = df['i_deg'].to_list()
w_list = df['w_deg'].to_list()
W_list = df['W_deg'].to_list()
M_list = df['M_deg'].to_list()
# # make new dataframe with reordered columns
df2 = pd.DataFrame()
df2['Packed MPC designation'] = df['packed_designation']
df2['Semimajor axis au barycentric'] = a_list
df2['Eccentricity barycentric'] = e_list
df2['Inclination ecliptic J2000 barycentric degrees'] = i_list
df2['Longitude of ascending node ecliptic J2000 barycentric degrees'] = W_list
df2['Argument of perihelion ecliptic J2000 barycentric degrees'] = w_list
df2['Mean anomaly ecliptic J2000 barycentric degrees'] = M_list
df2['Epoch JD'] = [JD for i in range(len(a_list))]
df2['Barycentric element source'] = ['JPL Horizons via Astroquery' for i in range(len(a_list))]
df2.to_csv('plutinos_for_mnras_noclones.csv',index=False)
thisdir = os.getcwd()
files_in_dir = os.listdir(thisdir)
for this_file in files_in_dir:
    if this_file.endswith('.out'):
        os.remove(os.path.join(thisdir,this_file))
