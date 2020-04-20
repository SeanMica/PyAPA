import sys
from phasepicker import aicdpicker
from associator import assoc3D
from datetime import datetime
from obspy.core import Stream
import associator.pickingSrc_obspy_sean as src
import associator.tables3D as tables3D
import configparser
import os
import time
import pathlib

# get AFile name from input
AFile=sys.argv[1]

# Configure file location
Config_File = f"{pathlib.Path(__file__).parent.absolute()}/PyAP.ini"

# reading Configure file
config = configparser.ConfigParser()
config.read(Config_File)

# some config
zero_phase = False

# Initializing all database will be used
DATABASE = src.Auto_Picking_Initializing(Config_File=Config_File, use_af_db=False)

# Initializing pre-processing module
PRE_PROCESS = src.Pre_processing(A_File = AFile, DB = DATABASE, config = config)

# Pre-processing
st, st_intensity, stmag, FirstStn, fileHeader = PRE_PROCESS.pre_picking_processing()

# copy Stream for small earthquake use
st_small = st.copy()

print(st)

# starting running picking sequence
now=time.time()
print('start picking...')
BB = ['Ch7','Ch8','Ch9']

# Define our picker instance
picker = aicdpicker.AICDPicker(t_ma = 3, nsigma = 6, t_up = 0.78, nr_len = 2, nr_coeff = 2, pol_len = 10, pol_coeff = 10, uncert_coeff = 3)

for tr in st:
    gap_check=tr.copy()
    tr.detrend('linear')
    if tr.stats.channel in BB:
        tr.filter("bandpass",freqmin=float(config['LargeEvt']['BandPassLow']),freqmax=float(config['LargeEvt']['BandPassHigh']),zerophase=zero_phase)
    scnl,picks,polarity,snr,uncert=picker.picks(tr)
    t_create=datetime.utcnow() # Record the time we made the picks
    # Add each pick to the database
    for i in range(len(picks)):
         
        if picks[i]-tr.stats.starttime >= 2.0 and tr.stats.endtime-picks[i] >= 2.0:
            if not src.check_pick_in_gap(gap_check,picks[i]):
                new_pick=tables3D.Pick(scnl,picks[i].datetime,'',snr[i],uncert[i],t_create)
                DATABASE.db_assoc.add(new_pick) # Add pick i to the database
DATABASE.db_assoc.commit() # Commit the pick to the database
print('finish picking... in ', time.time() - now)

# associate picks with phase types
associator = assoc3D.LocalAssociator(DATABASE.db_assoc, \
                                 DATABASE.db_nsta, \
                                 DATABASE.db_tt_curve, \
                                 fileHeader, \
                                 max_Parr = int(config['LargeEvt']['DB_Max_Parr']), \
                                 aggregation = 1, \
                                 aggr_norm = 'FA', \
                                 assoc_ot_uncert = 3, \
                                 nsta_declare=int(config['LargeEvt']['NSTA_Declare']), \
                                 config_par=config)
# candidate events
print('candidate events')
associator.id_candidate_events()
# associate candidates
print('associate candidates')
associator.associate_candidates()

# output to Pfile
print('output to Pfile')
Pfile = src.Pfile(stmag, st_intensity, DATABASE.db_assoc, DATABASE.db_nsta, Config_File=config, Vector_Mode=True)
ids = Pfile.get_assoc_ids()
# show possible earthquake event in database
print(ids)

# with out earthquake location program, remove code after this line
#QA_STATUS = False
#for id in ids:
#    P_header,zerodate,zero_min = Pfile.getNewHeader(id.id, fileHeader)
#    pfilename, QA = Pfile.writePfile(id.id,AFile,P_header,zerodate,zero_min,stmag,hypo3d=True)
#    print(QA)
#    if not QA:
#        os.remove(pfilename)
#    else:
#        QA_STATUS = True
    
