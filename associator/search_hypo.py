#from .tt_stations_3D import *
from .tables3D import *
from .table_nsta24 import *
from subprocess import DEVNULL, STDOUT, call
import numpy as np
import configparser
import os

def HYPO3D_Searching(session_tt_3D, session_nsta, matches, zerotime, config = None):
    # associator db
    assoc_db = session_tt_3D
    nsta_db  = session_nsta
    QA = True
    evt_lat = 23.5
    evt_lon = 120.0
    evt_dep = 30.0
    
    # make header 
    #header = make_header(zerotime)

    # first HYPO3D
    tpfilename = '24312359.T00'
    header = write_tmp_pfile(assoc_db, nsta_db, matches, zerotime, tpfilename)
    
    # first hypo
    call([config['Misc']['HypoRun_Path'], "-m", "1", "-P", tpfilename], stdout=DEVNULL, stderr=STDOUT)
    for i in range(5):
        # if hypo3d error, return failed (False)
        if call([config['Misc']['Hypo3D_Path'], tpfilename], stdout=DEVNULL, stderr=STDOUT) != 0:
        #call([config['Misc']['HypoRun_Path'], "-m", "1", "-P", tpfilename], stdout=DEVNULL, stderr=STDOUT)
        #if call([config['Misc']['HypoRun_Path'], "-m", "1", "-P", tpfilename], stdout=DEVNULL, stderr=STDOUT) != 0:
            QA = False
            return matches, evt_lat, evt_lon, evt_dep, QA
    
    # read first hypo result
    pfile = open(tpfilename,'r')
    line_i = True
    stn_used = []
    for line in pfile:
        if line_i:
            line_i = False
        else:
            P_resi = abs(float(line[29:34]))
            S_resi = abs(float(line[45:50]))
            #print("HYPO", line[1:5].rstrip(), P_resi, S_resi)
            if P_resi <= 3.0 and S_resi <= 3.0:
                stn_used.append(line[1:5].rstrip())
    pfile.close()
    
    new_matches = []
    for match in matches:
        if match[4].sta in stn_used:
            new_matches.append(match)
    
    # second hypo
    if len(new_matches) < 3:
        QA = False
        return matches, evt_lat, evt_lon, evt_dep, QA
    
    header = write_tmp_pfile(assoc_db, nsta_db, new_matches, zerotime, tpfilename)
    call([config['Misc']['HypoRun_Path'], "-m", "1", "-P", tpfilename], stdout=DEVNULL, stderr=STDOUT)
    for i in range(5):
        # if hypo3d error, return failed (False)
        if call([config['Misc']['Hypo3D_Path'], tpfilename], stdout=DEVNULL, stderr=STDOUT) != 0:
        #call([config['Misc']['HypoRun_Path'], "-m", "1", "-P", tpfilename], stdout=DEVNULL, stderr=STDOUT)
        #if call([config['Misc']['HypoRun_Path'], "-m", "1", "-P", tpfilename], stdout=DEVNULL, stderr=STDOUT) != 0:
            QA = False
            return matches, evt_lat, evt_lon, evt_dep, QA
    
    # read second hypo result
    pfile = open(tpfilename,'r')
    line_i = True
    stn_used = []
    for line in pfile:
        if line_i:
            if header[0:71] == line[0:71]:
                QA = False
                return matches, evt_lat, evt_lon, evt_dep, QA
            
            #if len(line[0:71]) >= 73 and line[0:71] == '-':
            #    QA = False
            #    return matches, evt_lat, evt_lon, evt_dep, QA
                
            lat_deg = float(line[19:21])
            lat_min = float(line[21:26])
            evt_lat = lat_deg+lat_min/60.0
            
            lon_deg = float(line[26:29])
            lon_min = float(line[29:34])
            evt_lon = lon_deg+lon_min/60.0
            
            evt_dep = float(line[34:40])
            
            # check TPfile ERZ, ERH, GAP, STN, add here
            STN = int(line[44:46])
            GAP = int(line[51:54])
            RMS = float(line[54:58])
            ERH = float(line[58:62])
            ERZ = float(line[62:66])
            
            if RMS > 9 or ERH > 10 or ERZ > 10:
                QA = False
                return matches, evt_lat, evt_lon, evt_dep, QA
            elif STN <= 4 and GAP >= 300:
                QA = False
                return matches, evt_lat, evt_lon, evt_dep, QA
            
            line_i = False
        else:
            P_resi = abs(float(line[29:34]))
            S_resi = abs(float(line[45:50]))
            #print("HYPO", line[1:5].rstrip(), P_resi, S_resi)
            if P_resi <= 3.0 and S_resi <= 3.0:
                stn_used.append(line[1:5].rstrip())
    pfile.close()
    
    matches = []
    for match in new_matches:
        if match[4].sta in stn_used:
            matches.append(match)
    
    if len(matches) < 3:
        QA = False
    
    return matches, evt_lat, evt_lon, evt_dep, QA
    
def make_header(zerotime, LON, LAT):
    # give initial event location
    lon_deg = '%03d'%int(LON)
    lon_min = '%5.2f'%((LON-int(LON))/60.0)
    lat_deg = '%02d'%int(LAT)
    lat_min = '%5.2f'%((LAT-int(LAT))/60.0)
    dep     = '%6.2f'%(11.0)
    
    # initial time is zero min time
    yy ='%04d'%int(zerotime.year)
    mm ='%02d'%int(zerotime.month)
    dd ='%02d'%int(zerotime.day)
    hh ='%02d'%int(zerotime.hour)
    min='%02d'%int(zerotime.minute)
    ss ='%6.2f'%(float(zerotime.second)+float(zerotime.microsecond)/1000000.0)
    
    header = " %s%s%s%s%s%s%s%s%s%s%s0.00                           D"%(yy,mm,dd,hh,min,ss,lat_deg,lat_min,lon_deg,lon_min,dep)
    
    return header

def write_tmp_pfile(session_tt_3D, session_nsta, matches, zerotime, tpfilename):
    assoc_db = session_tt_3D
    nsta_db  = session_nsta
    
    STN_LON = []
    STN_LAT = [] 
    for match in matches:
        ot_candi = match[4]
        stn = ot_candi.sta
        LAT, LON = nsta_db.query(NSTATable.latitude, NSTATable.longitude).filter(NSTATable.sta==stn).first()
        STN_LON.append(LON)
        STN_LAT.append(LAT)
    INIT_LON = np.mean(np.array(STN_LON))
    INIT_LAT = np.mean(np.array(STN_LAT))
    
    header = make_header(zerotime, INIT_LON, INIT_LAT)
    zeromin = zerotime.minute
    
    pfile = open(tpfilename,'w')
    pfile.write('%s\n'%header)
    
    for match in matches:
        # station & phase infomation
        ot_candi = match[4]
        stn = ot_candi.sta
        P_id = ot_candi.p_modified_id
        S_id = ot_candi.s_modified_id
        
        # get pick time and infomation
        P_modified = assoc_db.query(PickModified).filter(PickModified.id==P_id).first()
        S_modified = assoc_db.query(PickModified).filter(PickModified.id==S_id).first()
        
        # get pick time
        P_arrival = (P_modified.time - zerotime.datetime).total_seconds()
        S_arrival = (S_modified.time - zerotime.datetime).total_seconds()
        
        # get phase weighting
        P_weight = find_phase_weighting(assoc_db, modified_id=P_id, phase='P')
        S_weight = find_phase_weighting(assoc_db, modified_id=S_id, phase='S')
        
        # output to Pfile
        #writefile = ' %-4s   0.0   0   0  %2s%6.2f  .00 %.2f%6.2f .00  %.2f  .00  .00  .01  .00\n'%(stn,zeromin,P_arrival,P_weight,S_arrival,S_weight
        pfile.write(' %-4s   0.0   0   0  %2s%6.2f  .00 %.2f%6.2f .00  %.2f  .00  .00  .00  .00\n'%(stn,zeromin,P_arrival,P_weight,S_arrival,S_weight))
    
    pfile.close()
    
    return header

def find_phase_weighting(session_tt_3D, modified_id=1, phase='P'):
    assoc_db = session_tt_3D
    
    weightings = assoc_db.query(Pick.snr).filter(Pick.modified_id==modified_id).all()
    weighting, = max(weightings)
    
    #print(weighting)
    if phase == 'P':
        if weighting >= 100:
            weight = 0.0
        elif weighting >= 30:
            weight = 1.0
        elif weighting >= 20:
            weight = 2.0
        elif weighting >= 5:
            weight = 3.0
        else:
            weight = 4.0
        
    elif phase == 'S':
        if weighting >= 80:
            weight = 0.0
        elif weighting >= 15:
            weight = 1.0
        elif weighting >= 10:
            weight = 2.0
        elif weighting >= 5:
            weight = 3.0
        else:
            weight = 4.0
    
    else:
        print('input ERROR, GRRRRRRRRR GOONS')
    
    return weight
