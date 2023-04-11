#%%
def add_category_lists_sv20classifer(int_dict,prediction,Detached_list,Resonant_list,Scattering_list,Classical_list):
# Takes in a dictionary of class names.
# Takes in a prediction (probabilities that an object belongs to each class).
# Takes in lists of probabilities assigned to membership of each class, for each prior object.
# Sorts out the probabilities and assigns them to the correct list,
# ie, Classical probability to the Classical list, Resonant probability to the Resonant list, etc.
    if int_dict[0] == 'Detached':
        Detached_list.append(prediction[0][0])
    elif int_dict[0] == 'Resonant':
        Resonant_list.append(prediction[0][0])
    elif int_dict[0] == 'Scattering':
        Scattering_list.append(prediction[0][0])
    elif int_dict[0] == 'Classical':
        Classical_list.append(prediction[0][0])
    else:
        raise Exception('Unknown format in int_dict.',int_dict)
    if int_dict[1] == 'Detached':
        Detached_list.append(prediction[0][1])
    elif int_dict[1] == 'Resonant':
        Resonant_list.append(prediction[0][1])
    elif int_dict[1] == 'Scattering':
        Scattering_list.append(prediction[0][1])
    elif int_dict[1] == 'Classical':
        Classical_list.append(prediction[0][1])
    else:
        raise Exception('Unknown format in int_dict.',int_dict)
    if int_dict[2] == 'Detached':
        Detached_list.append(prediction[0][2])
    elif int_dict[2] == 'Resonant':
        Resonant_list.append(prediction[0][2])
    elif int_dict[2] == 'Scattering':
        Scattering_list.append(prediction[0][2])
    elif int_dict[2] == 'Classical':
        Classical_list.append(prediction[0][2])
    else:
        raise Exception('Unknown format in int_dict.',int_dict)
    if int_dict[3] == 'Detached':
        Detached_list.append(prediction[0][3])
    elif int_dict[3] == 'Resonant':
        Resonant_list.append(prediction[0][3])
    elif int_dict[3] == 'Scattering':
        Scattering_list.append(prediction[0][3])
    elif int_dict[3] == 'Classical':
        Classical_list.append(prediction[0][3])
    else:
        raise Exception('Unknown format in int_dict.',int_dict)
    return Detached_list,Resonant_list,Scattering_list,Classical_list
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
def get_GMdict():
    GM_sun = 1
    GM_mercury = 1/6023600
    GM_venus = 1/408523.71
    GM_earthmoon = 1/328900.56
    GM_mars = 1/3098708
    GM_jupiter = 1/1047.3486
    GM_saturn = 1/3497.898
    GM_uranus = 1/22902.98
    GM_neptune = 1/19412.24
    GMdict = {'sun':GM_sun,'mercury':GM_mercury,'venus':GM_venus,\
              'earth':GM_earthmoon,'mars':GM_mars,'jupiter':GM_jupiter,\
              'saturn':GM_saturn,'uranus':GM_uranus,'neptune':GM_neptune}
    return GMdict
#%%
def get_sv20classifier():
# Machine learning classifier for Kuiper belt objects. Code is unchanged from
# code distributed with 'Machine learning classification of Kuiper belt populations'.
# Smullen/Volk 2020 (sv20), MNRAS 497:2, September 2020, pg 1391-1403,
# https://doi.org/10.1093/mnras/staa1935
# Code is found at https://github.com/rsmullen/KBO_Classifier
    import pandas as pd
    from sklearn.model_selection import train_test_split
    from sklearn.ensemble import GradientBoostingClassifier
    training_file = '00_KBO_features.csv'
    all_KBOs = pd.read_csv(training_file,skipinitialspace=True)
    secure_KBOs = all_KBOs[all_KBOs['Securely Classified']==True]
    all_types = list(set(secure_KBOs['Class']))
    types_dict = { all_types[i] : i for i in range( len(all_types) ) }
    int_dict = { i : all_types[i] for i in range( len(all_types) ) }
    classes = secure_KBOs['Class'].map(types_dict)
    features_train, features_test, classes_train, classes_test = train_test_split(secure_KBOs, classes, test_size=0.3, random_state=30)
    features_train.drop(['MPC ID', 'Securely Classified', 'Class'], axis=1, inplace=True)
    features_train = features_train.to_numpy()
    features_test.drop(['MPC ID', 'Securely Classified', 'Class'], axis=1, inplace=True)
    features_test = features_test.to_numpy()
    classifier = GradientBoostingClassifier( learning_rate=0.1, loss='deviance', max_depth=3, max_features='log2', n_estimators=130, random_state=30 )
    classifier.fit(features_train, classes_train)
    return classifier, int_dict, types_dict
#%%
def parsedata_sv20classifier(data,classifier,int_dict):
# This function computes the necessary features to classify a KBO.
# Data MUST be a 101 row x 6 column array.
# Columns are t, a, e, i, Omega, omega.
# Rows are different time outputs: MUST be 1000yr outputs, ie [0, 1E3, 2E3....99E3,100E3].
# Returns features for classification.
# This code is unchanged from sv20, ie
# Smullen, Rachel A., and Kathryn Volk.
# 'Machine learning classification of Kuiper belt populations.'
# Monthly Notices of the Royal Astronomical Society 497.2 (2020): 1391-1403.
# The classification simulation and data parsing is copied from sv20 code with
# minimal changes to not have to query Horizons through Rebound when setting up a sim,
# because Rebound's Horizons query is very slow.
    import numpy as np
    # Take stats of simulations.
    initials = data[0,1:] # a, e, i, Omega, omega
    finals = data[-1,1:]
    mins = np.amin(data[:,1:],axis = 0)
    maxes = np.amax(data[:,1:],axis = 0)
    dels = maxes-mins
    means = np.mean(data[:,1:],axis = 0)
    stdev = np.std(data[:,1:],axis = 0)
    # Take time derivatives.
    diffs = data[1:,:]-data[:-1,:]
    dxdt = diffs[:,1:]/diffs[:,0, np.newaxis] # Add on new axis to time to give same dimensionality as the numerator.
    mindxdt = np.amin(dxdt,axis = 0)
    meandxdt = np.mean(dxdt,axis = 0)
    maxdxdt = np.amax(dxdt,axis = 0)
    deldxdt = maxdxdt-mindxdt
    # Rearrange data into the order we want.
    arrs = [initials,finals,mins,means,maxes,stdev,dels,mindxdt,meandxdt,maxdxdt,deldxdt]
    inds = [0,1,2,3,4] # a, e, i, Omega, omega
    features = []
    ## Features contains all x values, then all y, etc: xi, xf, xmin, xmean, xmax, xsigma, Deltax, xdotmin, xdotmean, xdotmax.
    for i in inds:
        for a in arrs:
            features += [a[i]]
    features_out = np.array(features).reshape(1,-1) # Make sure features is a 2d array.
    prediction = classifier.predict_proba(features_out) # Predict the probabilities of class membership for object.
    if np.max(prediction) == prediction[0][0]:
        category = int_dict[0]
    elif np.max(prediction) == prediction[0][1]:
        category = int_dict[1]
    elif np.max(prediction) == prediction[0][2]:
        category = int_dict[2]
    elif np.max(prediction) == prediction[0][3]:
        category = int_dict[3]
    return category,prediction
#%%
# This file classifies all the clones of a single Kuiper belt object.
THIS_INSTANCE = 1
import time
import numpy as np
import pandas as pd
import rebound
t0 = time.time()
# Prime the machine learning classifier.
classifier, int_dict, types_dict = get_sv20classifier()
# Get list of objects to classify.
date = '20221012'
year = 2022
month = 1
day = 1
hour = 0
minute = 0
second = 0
datestr = datestr([year,month,day,hour,minute,second])
elements_file = 'horizons_barycentric_nom_noclones_' + date + '_' + datestr + '.csv'
df = pd.read_csv(elements_file)
des_list = df['packed_designation'].tolist()
des = des_list[THIS_INSTANCE-1]
# Read heliocentric orbital elements of clones.
clones_input_file = 'clones_' + des + '.csv'
df = pd.read_csv(clones_input_file)
Nclones = df.shape[0]
ePh     = df['ePh'].tolist() # eccentricity, Plutino, heliocentric
qPh_au  = df['qPh_au'].tolist() # perihelion distance, Plutino, heliocentric, au
tpPh_jd = df['tpPh_jd'].tolist() # time of perihelion passage, Plutino, heliocentric, Julian date TDB
WPh_deg = df['WPh_deg'].tolist() # longitude of ascending node, Plutino, heliocentric, degrees
wPh_deg = df['wPh_deg'].tolist() # argument of perihelion, Plutino, heliocentric, degrees
iPh_deg = df['iPh_deg'].tolist() # inclination, Plutino, heliocentric, degrees
# Read heliocentric orbital elements of planets.
planets_file = 'planets_for_' + des + '.csv'
df = pd.read_csv(planets_file)
epochP  = df['epochP_jd'][0] # epoch of orbital elements, Julian date
eJh     = df['eJh'][0] # Jupiter
qJh_au  = df['qJh_au'][0]
tpJh_jd = df['tpJh_jd'][0]
WJh_deg = df['WJh_deg'][0]
wJh_deg = df['wJh_deg'][0]
iJh_deg = df['iJh_deg'][0]
eSh     = df['eSh'][0] # Saturn
qSh_au  = df['qSh_au'][0]
tpSh_jd = df['tpSh_jd'][0]
WSh_deg = df['WSh_deg'][0]
wSh_deg = df['wSh_deg'][0]
iSh_deg = df['iSh_deg'][0]
eUh     = df['eUh'][0] # Uranus
qUh_au  = df['qUh_au'][0]
tpUh_jd = df['tpUh_jd'][0]
WUh_deg = df['WUh_deg'][0]
wUh_deg = df['wUh_deg'][0]
iUh_deg = df['iUh_deg'][0]
eNh     = df['eNh'][0] # Neptune
qNh_au  = df['qNh_au'][0]
tpNh_jd = df['tpNh_jd'][0]
WNh_deg = df['WNh_deg'][0]
wNh_deg = df['wNh_deg'][0]
iNh_deg = df['iNh_deg'][0]
# Get masses of outer planets.
GMdict = get_GMdict()
# Create lists to store probabilities of class membership for each clone.
clone_class_list = []
classical_list = []
resonant_list = []
scattering_list = []
detached_list = []
# Run a Rebound simulation for each clone and classify it.
for iclone in range(Nclones):
    sim = rebound.Simulation()
    sim.add(m = 1,hash = '0')
    sim.integrator = 'ias15'
    epochobj = epochP
    # Build simulation, outer planets first.
    eobj = eJh # Jupiter
    qobj = qJh_au
    tpobj = tpJh_jd
    Wobj = np.radians(np.mod(WJh_deg,360))
    wobj = np.radians(np.mod(wJh_deg,360))
    iobj = np.radians(iJh_deg)
    aobj = qobj/(1-eobj)
    dt = epochobj - tpobj # time since pericenter passage in days
    dt = dt * (2*np.pi/365.25) # convert days to yr/(2pi)
    n = np.sqrt(1/aobj**3) # radians / (yr/2pi)
    Mobj = np.mod(n*dt,2*np.pi) # radians
    sim.add(primary=sim.particles[0],m=GMdict['jupiter'],hash='jupiter',\
            a=aobj,e=eobj,inc=iobj,omega=wobj,Omega=Wobj,M=Mobj)
    eobj = eSh # Saturn
    qobj = qSh_au
    tpobj = tpSh_jd
    Wobj = np.radians(np.mod(WSh_deg,360))
    wobj = np.radians(np.mod(wSh_deg,360))
    iobj = np.radians(iSh_deg)
    aobj = qobj/(1-eobj)
    dt = epochobj - tpobj # time since pericenter passage in days
    dt = dt * (2*np.pi/365.25) # convert days to yr/(2pi)
    n = np.sqrt(1/aobj**3) # radians / (yr/2pi)
    Mobj = np.mod(n*dt,2*np.pi) # radians
    sim.add(primary=sim.particles[0],m=GMdict['saturn'],hash='saturn',\
            a=aobj,e=eobj,inc=iobj,omega=wobj,Omega=Wobj,M=Mobj)
    eobj = eUh # Uranus
    qobj = qUh_au
    tpobj = tpUh_jd
    Wobj = np.radians(np.mod(WUh_deg,360))
    wobj = np.radians(np.mod(wUh_deg,360))
    iobj = np.radians(iUh_deg)
    aobj = qobj/(1-eobj)
    dt = epochobj - tpobj # time since pericenter passage in days
    dt = dt * (2*np.pi/365.25) # convert days to yr/(2pi)
    n = np.sqrt(1/aobj**3) # radians / (yr/2pi)
    Mobj = np.mod(n*dt,2*np.pi) # radians
    sim.add(primary=sim.particles[0],m=GMdict['uranus'],hash='uranus',\
            a=aobj,e=eobj,inc=iobj,omega=wobj,Omega=Wobj,M=Mobj)
    eobj = eNh # Neptune
    qobj = qNh_au
    tpobj = tpNh_jd
    Wobj = np.radians(np.mod(WNh_deg,360))
    wobj = np.radians(np.mod(wNh_deg,360))
    iobj = np.radians(iNh_deg)
    aobj = qobj/(1-eobj)
    dt = epochobj - tpobj # time since pericenter passage in days
    dt = dt * (2*np.pi/365.25) # convert days to yr/(2pi)
    n = np.sqrt(1/aobj**3) # radians / (yr/2pi)
    Mobj = np.mod(n*dt,2*np.pi) # radians
    sim.add(primary=sim.particles[0],m=GMdict['neptune'],hash='neptune',\
            a=aobj,e=eobj,inc=iobj,omega=wobj,Omega=Wobj,M=Mobj)
    # Add the Plutino clone.
    eobj = ePh[iclone]
    qobj = qPh_au[iclone]
    tpobj = tpPh_jd[iclone]
    Wobj = np.radians(np.mod(WPh_deg[iclone],360))
    wobj = np.radians(np.mod(wPh_deg[iclone],360))
    iobj = np.radians(iPh_deg[iclone])
    aobj = qobj/(1-eobj)
    dt = epochobj - tpobj # time since pericenter passage in days
    dt = dt * (2*np.pi/365.25) # convert days to yr/(2pi)
    n = np.sqrt(1/aobj**3) # radians / (yr/2pi)
    Mobj = np.mod(n*dt,2*np.pi) # radians
    sim.add(primary=sim.particles[0],m=0,hash='plutino',\
            a=aobj,e=eobj,inc=iobj,omega=wobj,Omega=Wobj,M=Mobj)
    # Prepare to integrate simulation.
    sim.N_active = 5
    sim.move_to_com()
    time_outs = np.linspace(0,100E3,101)*2*np.pi # 100 kyr
    data = []
    for i,t in enumerate(time_outs):
        if t>0:
            sim.move_to_com()
            # Integrate to next output.
            sim.integrate(t, exact_finish_time=True)
        orbits = sim.calculate_orbits(primary=sim.calculate_com())
        o = orbits[-1] # take KBO
        # Save t, a, e, i, Omega, omega. Time in data needs to be in years, so divide by 2pi.
        step = np.array([t/2/np.pi, o.a, o.e, np.degrees(o.inc), np.degrees(o.Omega)%360, np.degrees(o.omega)%360])
        # Add step to data.
        if len(data)==0: data = step
        else: data = np.vstack((data,step))
    # Release memory so we don't accidentally keep integrating the same sim with a different clone.
    sim = None
    category,prediction = parsedata_sv20classifier(data,classifier,int_dict)
    clone_class_list.append(category)
    detached_list,resonant_list,scattering_list,classical_list = \
        add_category_lists_sv20classifer(int_dict,prediction,detached_list,\
                        resonant_list,scattering_list,classical_list)
# Make new, separate file of clone class lists.
df = pd.DataFrame()
df['category'] = clone_class_list
df['classical_probability'] = classical_list
df['resonant_probability'] = resonant_list
df['scattering_probability'] = scattering_list
df['detached_probability'] = detached_list
clones_output_file = 'class_lists_' + des + '.csv'
df.to_csv(clones_output_file,index=False)
t1 = time.time()
elapsed_time_hours = (t1-t0)/3600
print('took',elapsed_time_hours,'hrs',THIS_INSTANCE,des)
