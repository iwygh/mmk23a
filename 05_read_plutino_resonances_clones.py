#%%
def JDfun(year,month,day,hour,minute,second):
    # from Vallado pg 183, valid for yrs 1900 to 2100
    import numpy as np
    JD = 367*year - np.floor(7/4*(year+np.floor(1/12*(month+9)))) + \
        np.floor(275*month/9) + day + 1721013.5 + 1/24*(hour+1/60*(second/60+minute))
    return JD
#%% unpack Minor Planet Center MPC ID
def unpack(designation):
    packed_designation = designation.lstrip()
    packed_designation = packed_designation.rstrip()
    Nt = len(packed_designation)
    onedigit_dict = {'0':'0','1':'1','2':'2','3':'3','4':'4','5':'5','6':'6',\
        '7':'7','8':'8','9':'9','A':'10','B':'11','C':'12','D':'13','E':'14',\
        'F':'15','G':'16','H':'17','I':'18','J':'19','K':'20','L':'21','M':'22',\
        'N':'23','O':'24','P':'25','Q':'26','R':'27','S':'28','T':'29','U':'30',\
        'V':'31','W':'32','X':'33','Y':'34','Z':'35','a':'36','b':'37','c':'38',\
        'd':'39','e':'40','f':'41','g':'42','h':'43','i':'44','j':'45','k':'46',\
        'l':'47','m':'48','n':'49','o':'50','p':'51','q':'52','r':'53','s':'54',\
        't':'55','u':'56','v':'57','w':'58','x':'59','y':'60','z':'61'}
    century_dict = {'I':'18','J':'19','K':'20'}
    if Nt == 4:
        unpacked_designation = packed_designation
    if Nt == 5:
        onedigit = packed_designation[0]
        fourdigits = packed_designation[1:5]
        if onedigit in onedigit_dict:
            onedigit = onedigit_dict[onedigit]
            unpacked_designation = onedigit + fourdigits
        else:
            unpacked_designation = '-99999'
    elif Nt == 7:
        century = packed_designation[0]
        year = packed_designation[1:3]
        letter1 = packed_designation[3]
        code1 = packed_designation[4]
        code2 = packed_designation[5]
        letter2 = packed_designation[6]
        space = ' '
        if century in century_dict:
            century = century_dict[century]
            if code1 == '0' and code2 == '0':
                code1 = ''
                code2 = ''
            else:
                if code1 in onedigit_dict:
                    code1 = onedigit_dict[code1]
                    if code1 == '0':
                        code1 = ''
            unpacked_designation = century + year + space + letter1 + letter2 + \
                code1 + code2
        else:
            unpacked_designation = '-99999'
    unpacked_designation = unpacked_designation.lstrip()
    unpacked_designation = unpacked_designation.rstrip()
    return unpacked_designation
#%%
import pandas as pd
import numpy as np
import os
from astroquery.jplhorizons import Horizons
# First, delete any .out files (High Performance Computing output logs).
thisdir = os.getcwd()
files_in_dir = os.listdir(thisdir)
for this_file in files_in_dir:
    if this_file.endswith('.out'):
        os.remove(os.path.join(thisdir,this_file))
# Get epoch at which to retrieve barycentric orbital elements.
epoch_file = 'epoch_settings_clones.csv'
df = pd.read_csv(epoch_file)
year = df['year'][0]
month = df['month'][0]
day = df['day'][0]
hour = df['hour'][0]
minute = df['minute'][0]
second = df['second'][0]
center = '500@0' # solar system barycenter
path = os.getcwd()
JD = JDfun(year,month,day,hour,minute,second)
# Get list of Resonant objects to process.
res_files = 'resonant_objects_clones.csv'
df = pd.read_csv(res_files)
Nres = df.shape[0]
print('Number of Resonant objects = ',Nres)
# Initialize list of Plutinos and their barycentric orbital elements.
plut_des = []
aPb_au = [] # semimajor axis, Plutino, barycentric, au
ePb = [] # eccentricity, Plutino, barycentric
iPb_deg = [] # inclination, Plutino, barycentric, degrees
WPb_deg = [] # longitude of ascending node, Plutino, barycentric, degrees
wPb_deg = [] # argument of perihelion, Plutino, barycentric, degrees
MPb_deg = [] # mean anomaly, Plutino, barycentric, degrees
false_negative_list = [] # objects that look resonant in 3:2 but are tagged False
false_positive_list = [] # objects that don't look resonant in 3:2 but are tagged True
for ires in range(Nres):
    des = df['resonant_des'][ires]
    if des in false_negative_list:
        plot_file = 'plot_False_3_2_' + des + '_clones.pdf'
        if os.path.isfile(os.path.join(path,plot_file)):
            plut_des.append(des)
            unpacked = unpack(des)
            if des == 'D4340':
                unpacked = '9' # We want the Pluto-Charon barycenter, not the Pluto body center.
            obj = Horizons(id=unpacked,location=center,epochs=str(JD))
            el = obj.elements()
            aPb_au.append(float(el['a'])) # au
            ePb.append(float(el['e']))
            iPb_deg.append(float(el['incl'])) # degrees
            WPb_deg.append(np.mod(float(el['Omega']),360))
            wPb_deg.append(np.mod(float(el['w']),360))
            MPb_deg.append(np.mod(float(el['M']),360))
    elif des not in false_positive_list:
        plot_file = 'plot_True_3_2_' + des + '_clones.pdf'
        if os.path.isfile(os.path.join(path,plot_file)):
            plut_des.append(des)
            unpacked = unpack(des)
            if des == 'D4340':
                unpacked = '9' # We want the Pluto-Charon barycenter, not the Pluto body center.
            obj = Horizons(id=unpacked,location=center,epochs=str(JD))
            el = obj.elements()
            aPb_au.append(float(el['a'])) # au
            ePb.append(float(el['e']))
            iPb_deg.append(float(el['incl'])) # degrees
            WPb_deg.append(np.mod(float(el['Omega']),360))
            wPb_deg.append(np.mod(float(el['w']),360))
            MPb_deg.append(np.mod(float(el['M']),360))
# # make new dataframe with reordered columns
df2 = pd.DataFrame()
df2['Packed MPC designation'] = plut_des
df2['Semimajor axis au barycentric'] = aPb_au
df2['Eccentricity barycentric'] = ePb
df2['Inclination ecliptic J2000 barycentric degrees'] = iPb_deg
df2['Longitude of ascending node ecliptic J2000 barycentric degrees'] = WPb_deg
df2['Argument of perihelion ecliptic J2000 barycentric degrees'] = wPb_deg
df2['Mean anomaly ecliptic J2000 barycentric degrees'] = MPb_deg
df2['Epoch JD'] = [JD for i in range(len(aPb_au))]
df2['Barycentric element source'] = ['JPL Horizons via Astroquery' for i in range(len(aPb_au))]
df2.to_csv('plutinos_for_mnras_clones.csv',index=False)
Nplut = df2.shape[0]
print('Number of Plutinos = ',Nplut)
