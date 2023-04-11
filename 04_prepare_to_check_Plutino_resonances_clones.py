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
    import numpy as np
    # from Vallado pg 183, valid for yrs 1900 to 2100
    JD = 367*year - np.floor(7/4*(year+np.floor(1/12*(month+9)))) + \
        np.floor(275*month/9) + day + 1721013.5 + 1/24*(hour+1/60*(second/60+minute))
    return JD
#%%
# This file is intended to be run after '01_prepare_to_classify_clones.py' has
# been run locally, 'sbatch slurm_classify_clones.slurm' has been run on the Puma
# cluster, and the resulting folder has been downloaded from Puma. This file
# reads the clone classifications, and saves a list of the objects for which more than
# half of the clones are Resonant, along with the initial conditions of the first
# clone of each object to be classified as Resonant. Next, it prepares a set of
# array jobs to check whether those initial conditions are Resonant in the 3:2 MMR.
# After running this file, upload the entire folder to Puma and run
# 'sbatch slurm_check_plutino_resonances.slurm'.
import pandas as pd
import os
# First, delete the .out files cluttering up the folder from the previous run on Puma.
thisdir = os.getcwd()
files_in_dir = os.listdir(thisdir)
for this_file in files_in_dir:
    if this_file.endswith('.out'):
        os.remove(os.path.join(thisdir,this_file))
# Get list of objects to process.
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
df = pd.read_csv(horizons_file)
des_list = df['packed_designation'].tolist()
Nobj = df.shape[0]
# Initialize lists of Resonant objects and their first Resonant initial conditions.
resonant_des = []
ePh     = [] # eccentricity, Plutino, heliocentric
qPh_au  = [] # perihelion distance, Plutino, heliocentric, au
tpPh_jd = [] # time of perihelion passage, Plutino, heliocentric, Julian date TDB
WPh_deg = [] # longitude of ascending node, Plutino, heliocentric, degrees
wPh_deg = [] # argument of perihelion, Plutino, heliocentric, degrees
iPh_deg = [] # inclination, Plutino, heliocentric, degrees
# Go through the objects and add them to the lists if more than half of the clones are Resonant.
for iobj in range(Nobj):
    added_flag = 0 # Flags whether the first Resonant initial conditions have been added to the list yet.
    des = des_list[iobj]
    classes_file = 'class_lists_' + des + '.csv'
    df = pd.read_csv(classes_file)
    Nclones = df.shape[0]
    # Count Resonant clones.
    resonant_count = 0
    for iclone in range(Nclones):
        cat = df['category'][iclone]
        cat = cat.lstrip()
        cat = cat.rstrip()
        if cat == 'Resonant':
            resonant_count = resonant_count + 1
    # Add object to Resonant list if necessary.
    if resonant_count == Nclones:
        resonant_des.append(des)
        # Find first Resonant initial conditions and add to list.
        for iclone in range(Nclones):
            cat = df['category'][iclone]
            cat = cat.lstrip()
            cat = cat.rstrip()
            if cat == 'Resonant' and added_flag == 0:
                clones_file = 'clones_' + des + '.csv'
                df2 = pd.read_csv(clones_file)
                ePh.append(df2['ePh'][iclone])
                qPh_au.append(df2['qPh_au'][iclone])
                tpPh_jd.append(df2['tpPh_jd'][iclone])
                WPh_deg.append(df2['WPh_deg'][iclone])
                wPh_deg.append(df2['wPh_deg'][iclone])
                iPh_deg.append(df2['iPh_deg'][iclone])
                added_flag = 1
# Save file of Resonant objects and their first Resonant initial conditions.
df3 = pd.DataFrame()
df3['resonant_des'] = resonant_des
df3['ePh'] = ePh
df3['qPh_au'] = qPh_au
df3['tpPh_jd'] = tpPh_jd
df3['WPh_deg'] = WPh_deg
df3['wPh_deg'] = wPh_deg
df3['iPh_deg'] = iPh_deg
resonant_file = 'resonant_objects_clones.csv'
df3.to_csv(resonant_file,index=False)
Nres = df3.shape[0] # Number of Resonant objects
print(Nobj,Nres)
# Generate as many instances of check_plutino_resonances_#.py as there are objects to
# test for membership in the 3:2 MMR. Each instance will classify a separate object on the cluster.
template_file = '04a_check_Plutino_resonances_template_clones.py'
with open(template_file,'r') as file:
    template_data = file.read()
for i in range(Nres):
    ct = i + 1
    ctstr = str(ct)
    outfile = 'check_Plutino_resonances_clones_' + ctstr + '.py'
    filedata = template_data
    filedata = filedata.replace('THIS_INSTANCE = 1','THIS_INSTANCE = '+ctstr)
    with open(outfile,'w') as file:
        file.write(filedata)
# Generate a batch submission file for the cluster to classify the clones with.
# On the cluster, run "sbatch slurm_classify_clones.slurm".
template_file = '04b_slurm_check_Plutino_resonances_clones_template.slurm'
with open(template_file,'r') as file:
    template_data = file.read()
filedata = template_data
filedata = filedata.replace('Nobj',str(Nres))
outfile = 'slurm_check_plutino_resonances_clones.slurm'
with open(outfile,'w') as file:
    file.write(filedata)
