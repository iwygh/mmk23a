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
def generate_clones(JD,datestr,Nclones,infile):
# Read in covariance for each object.
# Generate clones for each object from the covariance.
# For each clone, save heliocentric elements of the object and the planets.
# These will be used later in building a simulation in Rebound for classifying
# the clones. It is easiest to build a simulation from heliocentric elements.
    import numpy as np
    import pandas as pd
    import json
    import urllib
    from astroquery.jplhorizons import Horizons
    # Start with the list of objects left as a rough cut by sbdb_reduce.
    # infile = 'sbdb_reduce_output_clones.csv'
    # Start a data frame with the objects left as a rough cut by sbdb_reduce.
    df = pd.read_csv(infile)
    Nobj = df.shape[0]
    center = '500@10' # easiest to build a simulation in heliocentric elements
    for iobj in range(Nobj):
        # initialize lists of heliocentric elements of the object and the planets
        epochP = []
        ePh = [] # eccentricity, Plutino, heliocentric
        qPh_au = [] # perihelion distance, Plutino, heliocentric, au
        tpPh_jd = [] # time of perihelion passage, Plutino, heliocentric, Julian date TDB
        WPh_deg = [] # longitude of ascending node, Plutino, heliocentric, degrees
        wPh_deg = [] # argument of perihelion, Plutino, heliocentric, degrees
        iPh_deg = [] # inclination, Plutino, heliocentric, degrees
        eJh = [] # Jupiter
        qJh_au = []
        tpJh_jd = []
        WJh_deg = []
        wJh_deg = []
        iJh_deg = []
        eSh = [] # Saturn
        qSh_au = []
        tpSh_jd = []
        WSh_deg = []
        wSh_deg = []
        iSh_deg = []
        eUh = [] # Uranus
        qUh_au = []
        tpUh_jd = []
        WUh_deg = []
        wUh_deg = []
        iUh_deg = []
        eNh = [] # Neptune
        qNh_au = []
        tpNh_jd = []
        WNh_deg = []
        wNh_deg = []
        iNh_deg = []
        # Get packed MPC designation of object.
        des = df['packed_designation'][iobj]
        # Contact JPL Solar System Dynamics Small Body Database API to get covariance.
        # Returns heliocentric orbels and covariance
        url = 'https://ssd-api.jpl.nasa.gov/sbdb.api?sstr='+des+'&full-prec=True&cov=mat'
        response = urllib.request.urlopen(url)
        data = json.loads(response.read())
        # Save covariance and other information to file.
        jsonfile = 'ssd_json_' + des + '.json'
        with open(jsonfile,'w') as of:
            json.dump(data,of)
        # Save heliocentric orbital elements of the Plutino and the planets.
        clones_file = 'clones_' + des + '.csv'
        planets_file = 'planets_for_' + des + '.csv'
        if des == 'D4340': # Pluto
            # Because Pluto's covariance is essentially zero, we save heliocentric orbels at
            # the arbitrary Julian date used as an input.
            # Save heliocentric orbital elements of the planets at the arbitrary epoch.
            for J in range(Nclones+1):
                epochP.append(JD)
                obj = Horizons(id='5',location=center,epochs=JD) # Jupiter barycenter
                el = obj.elements()
                eJh.append(float(el['e']))
                qJh_au.append(float(el['q']))
                tpJh_jd.append(float(el['Tp_jd']))
                WJh_deg.append(np.mod(float(el['Omega']),360))
                wJh_deg.append(np.mod(float(el['w']),360))
                iJh_deg.append(float(el['incl']))
                obj = Horizons(id='6',location=center,epochs=JD) # Saturn barycenter
                el = obj.elements()
                eSh.append(float(el['e']))
                qSh_au.append(float(el['q']))
                tpSh_jd.append(float(el['Tp_jd']))
                WSh_deg.append(np.mod(float(el['Omega']),360))
                wSh_deg.append(np.mod(float(el['w']),360))
                iSh_deg.append(float(el['incl']))
                obj = Horizons(id='7',location=center,epochs=JD) # Uranus barycenter
                el = obj.elements()
                eUh.append(float(el['e']))
                qUh_au.append(float(el['q']))
                tpUh_jd.append(float(el['Tp_jd']))
                WUh_deg.append(np.mod(float(el['Omega']),360))
                wUh_deg.append(np.mod(float(el['w']),360))
                iUh_deg.append(float(el['incl']))
                obj = Horizons(id='8',location=center,epochs=JD) # Neptune barycenter
                el = obj.elements()
                eNh.append(float(el['e']))
                qNh_au.append(float(el['q']))
                tpNh_jd.append(float(el['Tp_jd']))
                WNh_deg.append(np.mod(float(el['Omega']),360))
                wNh_deg.append(np.mod(float(el['w']),360))
                iNh_deg.append(float(el['incl']))
                dictionary = {'epochP_jd':epochP,
                      'eJh':eJh,'qJh_au':qJh_au,'tpJh_jd':tpJh_jd,'WJh_deg':WJh_deg,'wJh_deg':wJh_deg,'iJh_deg':iJh_deg,\
                      'eSh':eSh,'qSh_au':qSh_au,'tpSh_jd':tpSh_jd,'WSh_deg':WSh_deg,'wSh_deg':wSh_deg,'iSh_deg':iSh_deg,\
                      'eUh':eUh,'qUh_au':qUh_au,'tpUh_jd':tpUh_jd,'WUh_deg':WUh_deg,'wUh_deg':wUh_deg,'iUh_deg':iUh_deg,\
                      'eNh':eNh,'qNh_au':qNh_au,'tpNh_jd':tpNh_jd,'WNh_deg':WNh_deg,'wNh_deg':wNh_deg,'iNh_deg':iNh_deg}
                df_out = pd.DataFrame(dictionary)
                df_out.to_csv(planets_file,index=False)
                # The first line of the clones_[des].csv file is always the nominal orbit at epoch.
                obj = Horizons(id='9',location=center,epochs=JD) # Pluto-Charon barycenter
                el = obj.elements()
                ePh.append(float(el['e']))
                qPh_au.append(float(el['q'])) # au
                tpPh_jd.append(float(el['Tp_jd'])) # time of pericenter passage, Julian date
                WPh_deg.append(np.mod(float(el['Omega']),360)) # degrees
                wPh_deg.append(np.mod(float(el['w']),360)) # degrees
                iPh_deg.append(float(el['incl'])) # degrees
            # all Pluto clones are exactly the same
            dictionary = {'ePh':ePh,'qPh_au':qPh_au,'tpPh_jd':tpPh_jd,'WPh_deg':WPh_deg,'wPh_deg':wPh_deg,'iPh_deg':iPh_deg}
            df_out = pd.DataFrame(dictionary)
            df_out.to_csv(clones_file,index=False)
        else:
            # Now parse JSON files for covariance information.
            infile = 'ssd_json_' + des + '.json'
            with open(infile) as json_file:
                data = json.load(json_file)
            cov_array = data['orbit']['covariance']['data'] # orbels in order e,q,tp,node,peri,i
            # Row-major or column-major interpretation of covariance matrix doesn't matter,
            # because covariance matrix is symmetric.
            cov_array_2 = np.array([\
                [float(cov_array[0][0]),float(cov_array[0][1]),float(cov_array[0][2]),\
                     float(cov_array[0][3]),float(cov_array[0][4]),float(cov_array[0][5])],\
                [float(cov_array[1][0]),float(cov_array[1][1]),float(cov_array[1][2]),\
                     float(cov_array[1][3]),float(cov_array[1][4]),float(cov_array[1][5])],\
                [float(cov_array[2][0]),float(cov_array[2][1]),float(cov_array[2][2]),\
                     float(cov_array[2][3]),float(cov_array[2][4]),float(cov_array[2][5])],\
                [float(cov_array[3][0]),float(cov_array[3][1]),float(cov_array[3][2]),\
                     float(cov_array[3][3]),float(cov_array[3][4]),float(cov_array[3][5])],\
                [float(cov_array[4][0]),float(cov_array[4][1]),float(cov_array[4][2]),\
                     float(cov_array[4][3]),float(cov_array[4][4]),float(cov_array[4][5])],\
                [float(cov_array[5][0]),float(cov_array[5][1]),float(cov_array[5][2]),\
                     float(cov_array[5][3]),float(cov_array[5][4]),float(cov_array[5][5])] ])
            # Covariance and elements are specified at this epoch (jd).
            cov_epoch = data['orbit']['covariance']['epoch']
            cov_epoch = float(cov_epoch)
            # Save heliocentric orbital elements of the planets at the covariance epoch.
            epochP.append(cov_epoch)
            obj = Horizons(id='5',location=center,epochs=cov_epoch) # Jupiter barycenter
            el = obj.elements()
            eJh.append(float(el['e']))
            qJh_au.append(float(el['q']))
            tpJh_jd.append(float(el['Tp_jd']))
            WJh_deg.append(np.mod(float(el['Omega']),360))
            wJh_deg.append(np.mod(float(el['w']),360))
            iJh_deg.append(float(el['incl']))
            obj = Horizons(id='6',location=center,epochs=cov_epoch) # Saturn barycenter
            el = obj.elements()
            eSh.append(float(el['e']))
            qSh_au.append(float(el['q']))
            tpSh_jd.append(float(el['Tp_jd']))
            WSh_deg.append(np.mod(float(el['Omega']),360))
            wSh_deg.append(np.mod(float(el['w']),360))
            iSh_deg.append(float(el['incl']))
            obj = Horizons(id='7',location=center,epochs=cov_epoch) # Uranus barycenter
            el = obj.elements()
            eUh.append(float(el['e']))
            qUh_au.append(float(el['q']))
            tpUh_jd.append(float(el['Tp_jd']))
            WUh_deg.append(np.mod(float(el['Omega']),360))
            wUh_deg.append(np.mod(float(el['w']),360))
            iUh_deg.append(float(el['incl']))
            obj = Horizons(id='8',location=center,epochs=cov_epoch) # Neptune barycenter
            el = obj.elements()
            eNh.append(float(el['e']))
            qNh_au.append(float(el['q']))
            tpNh_jd.append(float(el['Tp_jd']))
            WNh_deg.append(np.mod(float(el['Omega']),360))
            wNh_deg.append(np.mod(float(el['w']),360))
            iNh_deg.append(float(el['incl']))
            dictionary = {'epochP_jd':epochP,
                  'eJh':eJh,'qJh_au':qJh_au,'tpJh_jd':tpJh_jd,'WJh_deg':WJh_deg,'wJh_deg':wJh_deg,'iJh_deg':iJh_deg,\
                  'eSh':eSh,'qSh_au':qSh_au,'tpSh_jd':tpSh_jd,'WSh_deg':WSh_deg,'wSh_deg':wSh_deg,'iSh_deg':iSh_deg,\
                  'eUh':eUh,'qUh_au':qUh_au,'tpUh_jd':tpUh_jd,'WUh_deg':WUh_deg,'wUh_deg':wUh_deg,'iUh_deg':iUh_deg,\
                  'eNh':eNh,'qNh_au':qNh_au,'tpNh_jd':tpNh_jd,'WNh_deg':WNh_deg,'wNh_deg':wNh_deg,'iNh_deg':iNh_deg}
            df_out = pd.DataFrame(dictionary)
            df_out.to_csv(planets_file,index=False)
            # Now retrieve the nominal orbital elements at epoch.
            e_nom  = data['orbit']['covariance']['elements'][0]['value']
            q_nom  = data['orbit']['covariance']['elements'][1]['value']
            tp_nom = data['orbit']['covariance']['elements'][2]['value']
            W_nom  = data['orbit']['covariance']['elements'][3]['value']
            w_nom  = data['orbit']['covariance']['elements'][4]['value']
            i_nom  = data['orbit']['covariance']['elements'][5]['value']
            e_nom  = float(e_nom)
            q_nom  = float(q_nom) # au
            tp_nom = float(tp_nom) # time of pericenter passage, Julian date
            W_nom  = float(W_nom) # degrees
            w_nom  = float(w_nom) # degrees
            i_nom  = float(i_nom) # degrees
            # The first line of clones_[des].csv is always the heliocentric orbels
            # of the Plutino and the planets at the epoch of the covariance matrix.
            ePh.append(e_nom)
            qPh_au.append(q_nom)
            tpPh_jd.append(tp_nom)
            WPh_deg.append(W_nom)
            wPh_deg.append(w_nom)
            iPh_deg.append(i_nom)
            epochP.append(cov_epoch)
            # Make nominal state vector before generating clones.
            state_nom = np.array([e_nom,q_nom,tp_nom,W_nom,w_nom,i_nom])
            # Generate clones.
            for j in range(Nclones):
                L = np.linalg.cholesky(cov_array_2)
                clone_1 = np.random.normal(0,1,6) # 6-element random Gaussian mean 0 stdev 1
                clone_2 = np.dot(L,clone_1) # Gaussian perturbations from zero correlated according to cov matrix
                clone_3 = clone_2 + state_nom # add perturbation to nominal orbit
                e_clone = clone_3[0]
                q_clone = clone_3[1]
                tp_clone = clone_3[2]
                W_clone = clone_3[3]
                w_clone = clone_3[4]
                i_clone = clone_3[5]
                if (e_clone<=0 or q_clone<=0): # don't want unphysical orbits
                    print('clone problems with ',iobj,des)
                W_clone = np.mod(W_clone,360)
                w_clone = np.mod(w_clone,360)
                # Record heliocentric orbels of Plutino.
                ePh.append(e_clone)
                qPh_au.append(q_clone)
                tpPh_jd.append(tp_clone)
                WPh_deg.append(W_clone)
                wPh_deg.append(w_clone)
                iPh_deg.append(i_clone)
        dictionary = {'ePh':ePh,'qPh_au':qPh_au,'tpPh_jd':tpPh_jd,'WPh_deg':WPh_deg,'wPh_deg':wPh_deg,'iPh_deg':iPh_deg}
        df_out = pd.DataFrame(dictionary)
        df_out.to_csv(clones_file,index=False)
    return
#%%
def JDfun(year,month,day,hour,minute,second):
    # from Vallado pg 183, valid for yrs 1900 to 2100
    import numpy as np
    JD = 367*year - np.floor(7/4*(year+np.floor(1/12*(month+9)))) + \
        np.floor(275*month/9) + day + 1721013.5 + 1/24*(hour+1/60*(second/60+minute))
    return JD
#%%
import pandas as pd
import os
date = '20221012'
year = 2022
month = 1
day = 1
hour = 0
minute = 0
second = 0
Nclones = 300
# Save settings for epoch and number of clones.
dictionary = {'year':[year],'month':[month],'day':[day],'hour':[hour],'minute':[minute],\
              'second':[second],'Nclones':[Nclones]}
df_out = pd.DataFrame(dictionary)
df_out.to_csv('epoch_settings_clones.csv',index=False)
# Make a string of the date, for use in filenames.
datestr = datestr([year,month,day,hour,minute,second])
# We will need the Julian date to look up orbital elements at the epoch
# chosen for making population-level figures and statistics.
JD = JDfun(year,month,day,hour,minute,second)
elements_file = 'horizons_barycentric_nom_noclones_' + date + '_' + datestr + '.csv'
df_elements = pd.read_csv(elements_file)
Njobs = df_elements.shape[0]
generate_clones(JD,datestr,Nclones,elements_file)
# Generate as many instances of classify_clones_#.py as there are objects to
# test. Each instance will classify a separate object on the cluster.
template_file = '03a_classify_clones_template.py'
with open(template_file,'r') as file:
    template_data = file.read()
for i in range(Njobs):
    ct = i + 1
    ctstr = str(ct)
    outfile = 'classify_clones_' + ctstr + '.py'
    filedata = template_data
    filedata = filedata.replace('THIS_INSTANCE = 1','THIS_INSTANCE = '+str(int(i+1)))
    with open(outfile,'w') as file:
        file.write(filedata)
# Generate a batch submission file for the cluster to classify the clones with.
# On the cluster, run "sbatch slurm_classify_clones.slurm".
template_file = '03b_slurm_classify_clones_template.slurm'
with open(template_file,'r') as file:
    template_data = file.read()
filedata = template_data
filedata = filedata.replace('NJOBS',str(Njobs))
outfile = 'slurm_classify_clones.slurm'
with open(outfile,'w') as file:
    file.write(filedata)
thisdir = os.getcwd()
files_in_dir = os.listdir(thisdir)
for this_file in files_in_dir:
    if this_file.endswith('.out'):
        os.remove(os.path.join(thisdir,this_file))
