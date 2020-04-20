#import scipy
from scipy.fftpack import rfft, irfft, rfftfreq
import numpy as np
from numpy import linalg as LA
from obspy import Trace, Stream
from obspy.core import UTCDateTime
from obspy.geodetics import locations2degrees, degrees2kilometers
from sqlalchemy.orm import *
from sqlalchemy import create_engine
from sqlalchemy import or_

# database format import
from . import table_nsta24 as table_nsta24
from . import table_tt_curve as table_tt_curve
from . import tables3D as tables3D

from .table_nsta24 import *
from .table_tt_curve import *
from .tables3D import *

# import 3d travel-time calculator
#from .pbr import pbr

import os
import sys
import struct
import csv
from datetime import *
import math
import configparser
from shutil import copyfile

# put all pre-processing presure in one
class Pre_processing():
    def __init__(self, A_File = None, DB = None, config = None, EEW_Mode = False):
        if not A_File:
            sys.exit("A_File not given ---- class Pre_processing")
        elif not os.path.isfile(A_File):
            sys.exit("A_File not exist ---- class Pre_processing")
        else:
            self.A_File                 = A_File
        
        if not DB:
            sys.exit("Database not given ---- class Pre_processing")
        else:
            self.DB                     = DB
        
        if not config:
            sys.exit("Config not given ---- class Pre_processing")
        else:
            self.config = config
        
        # EEW mode
        self.EEW_Mode = EEW_Mode
        
    def pre_picking_processing(self):
        # read Afile
        self.st, self.FirstStn, self.fileHeader = self.unpackAfile(self.A_File)
        
        # remove station not in nsta24 database
        self.checkNonDBsta(self.st, self.DB.db_nsta)
        
        # rotate OBS stations
        if not self.EEW_Mode:
            OBS_rotate(self.st, Config_File = self.config)
        
        # remove station in DONT_USE_LIST
        self.remove_sta_not_use(self.st, checklist = self.config['Misc']['CheckList_sta'])
        
        # remove station in addon block list
        self.remove_addon_block_list(self.st, blocklist = self.config['Misc']['Addon_BlockList'])
        
        # remove stations is outside velocty model
        self.Remove_Stn_Outside_Grid(self.st, nsta_db = self.DB.db_nsta)
        
        # de-mean without zeros
        self.demean_wo_zero()

        # copy trace for intensity use
        st_intensity_tmp  = Stream()
        st_intensity_tmp += self.st.select(channel = 'Ch1')
        st_intensity_tmp += self.st.select(channel = 'Ch2')
        st_intensity_tmp += self.st.select(channel = 'Ch3')
        self.st_intensity = st_intensity_tmp.copy()
        
        # copy streams for calculate Magnitude
        self.stmag = self.st.copy()
        
        # remove no data Streams
        self.checkZeroGap_new(self.st)
        
        return self.st, self.st_intensity, self.stmag, self.FirstStn, self.fileHeader
        
    """
    Unpack Binary File ->12bit
    Unpack New Binary File -> 24bit
    """
    def unpackAfile(self, infile):
    
    # == opening Afile ==
        b= os.path.getsize(infile)
        FH = open(infile, 'rb')
        line = FH.read(b)
        fileHeader= struct.unpack("<4s3h6bh6s", line[0:24])
        
        fileLength = fileHeader[3]
        port = fileHeader[10]
        FirstStn = fileHeader[11][0:4].decode('ASCII').rstrip()
        print(fileHeader)
    # =================================Header===================================
        
        portHeader = []
        for i in range(24,port*32,32):
            port_data = struct.unpack("<4s4s3sbh2b4s12b",line[i:i+32])
            portHeader.append(port_data)
    
    # =================================Data===================================    
    
        dataStartByte = 24+int(port)*32
        dataPoint = 3*int(port)*int(fileLength)*100
        times = int(port)*3*4
        data=[]
    
        data = struct.unpack("<%di"%dataPoint,line[dataStartByte:dataStartByte+dataPoint*4])
    
        
        portHeader = np.array(portHeader)    
        data = np.array(data)    
        idata =data.reshape((3,port,fileLength*100),order='F')
        
    #== write to obspy Stream --
        #print(fileHeader)
        #print(len(idata[0][0]))
        sttime = UTCDateTime(fileHeader[1],fileHeader[4],fileHeader[5],fileHeader[6],fileHeader[7],fileHeader[8],fileHeader[2])
        npts = fileHeader[3]*fileHeader[9]
        samp = fileHeader[9]
        #print(sttime)
        # afst = Afile's Stream
        afst = Stream()
        
        for stc in range(fileHeader[10]):
            stn = portHeader[stc][0].decode('ASCII').rstrip()
            instrument = portHeader[stc][1].decode('ASCII').rstrip()
            loc = '0'+str(portHeader[stc][6].decode('ASCII'))
            #net = "TW"
            net = str(portHeader[stc][7].decode('ASCII')).rstrip()
            GPS = int(portHeader[stc][3])
            
            # remove GPS unlock or broken station
            if ( GPS == 1 or GPS == 2 ):
                chc = 0
                if instrument == 'FBA':
                    chc = 1
                elif instrument == 'SP':
                    chc = 4
                elif instrument == 'BB':
                    chc = 7
                
                #print(chc,instrument)
                
                # for each channel in port
                for ch in range(3):
                    #print(num,ch,chc)
                    chn = 'Ch'+str(chc+ch)
                    #print(stc,channel)
                    
                    stats = {'network': net, 'station': stn, 'location': loc,
                            'channel': chn, 'npts': npts, 'sampling_rate': samp,
                            'starttime': sttime}
                    
                    data = np.array(idata[ch][stc], dtype=float)
                    sttmp = Stream([Trace(data=data, header=stats)])
                    afst += sttmp
        
        return afst, FirstStn, fileHeader
        
    def checkNonDBsta(self, InputStream, db_nsta):
        #== check station in database
        for tr in InputStream:
            sta = tr.stats['station']
            stadb = db_nsta.query(NSTATable.id).filter(NSTATable.sta==sta).first()
            #print(stadb)
            
            if stadb == None:
                #print(sta,'will remove')
                InputStream.remove(tr)
        
    def remove_sta_not_use(self, InputStream, checklist = None):
        if checklist:
            not_use_sta = []
            csvfile = open(checklist,'r')
            for row in csv.reader(csvfile, delimiter=' '):
                if row[-1] == '0':
                    not_use_sta.append(row[0])
            
            for tr in InputStream:
                sta = tr.stats['station']
                if sta in not_use_sta:
                    InputStream.remove(tr)
        else:
            sys.exit("checklist not given ---- remove_sta_not_use")

    def remove_addon_block_list(self, InputStream, blocklist = None):
        if blocklist:
        
            BL_DAT = open(blocklist, 'r')
            
            for line in BL_DAT:
                # read block channel from file
                STN = line[0:4].rstrip()
                CHN = line[5:9].rstrip()
                NET = line[10:14].rstrip()
                LOC = '0'+line[15:].rstrip()
                
                if CHN == 'FBA':
                    CHN_DEL = ('Ch1', 'Ch2', 'Ch3')
                elif CHN == 'SP':
                    CHN_DEL = ('Ch4', 'Ch5', 'Ch6')
                elif CHN == 'BB':
                    CHN_DEL = ('Ch7', 'Ch8', 'Ch9')
                else:
                    print('input error in remove_addon_block_list')
                    exit()
                
                for CH in CHN_DEL:
                    st = InputStream.select(station=STN,channel=CH,network=NET,location=LOC)
                    
                    for tr in st:
                        InputStream.remove(tr)
        else:
            sys.exit("blocklist not given ---- remove_addon_block_list")

    def Remove_Stn_Outside_Grid(self, InputStream, nsta_db = None, min_lon = 115, max_lon = 130, min_lat = 15, max_lat = 30):
        #check station distance
        for tr in InputStream:
            stn = tr.stats['station']
            #print(stn)
            tr_lon, tr_lat = nsta_db.query(NSTATable.longitude, NSTATable.latitude).filter(NSTATable.sta==stn).first()
            
            if tr_lon > max_lon or tr_lon < min_lon or tr_lat > max_lat or tr_lat < min_lat:
                InputStream.remove(tr)
        
    def checkZeroGap_new(self, InputStream):
        for tr in InputStream:
            zerocount = 0
            idata = tr.data
            
            if np.count_nonzero(idata) == 0:
                InputStream.remove(tr)
                break
            
            #idata_mask = idata.astype(bool)
            #tmp_data = np.ma.masked_array(idata,mask=idata_mask)
            tmp_data = np.ma.masked_where(idata == 0, idata)
            mean = tmp_data.mean()
            tmp_data = tmp_data - mean
            tmp_data = tmp_data.filled(0)
            
            tr.data = tmp_data
    
    def demean_wo_zero(self):
        # demean with out zero value (gap) in trace
        for tr in self.st:
            idata = tr.data
            zerocount = np.count_nonzero(idata)

            if zerocount == 0:
                self.st.remove(tr)
                continue

            #idata_mask = idata.astype(bool)
            #tmp_data = np.ma.masked_array(idata,mask=idata_mask)
            tmp_data = np.ma.masked_where(idata == 0, idata)
            mean = tmp_data.mean()
            tmp_data = tmp_data - mean
            tmp_data = tmp_data.filled(0)
            tr.data = tmp_data

def check_pick_in_gap(inputStream, pick, time_window=0.2):
    
    result = False
    
    tiny_st = inputStream.slice(pick-time_window, pick+time_window)
    
    if np.count_nonzero(tiny_st.data==0) > 4:
    #if np.count_nonzero(tiny_st.data) == 0:
        result = True
    
    return result
#
def RemoveStnFarAway(InputStream,db_nsta,FirstStn,MaxDeg):
    
    fs_lon, fs_lat = db_nsta.query(NSTATable.longitude, NSTATable.latitude).filter(NSTATable.sta==FirstStn).first()
    
    #check station distance
    for tr in InputStream:
        stn = tr.stats['station']
        tr_lon, tr_lat = db_nsta.query(NSTATable.longitude, NSTATable.latitude).filter(NSTATable.sta==stn).first()
        
        EpicDist=locations2degrees(fs_lat, fs_lon, tr_lat, tr_lon)
        
        if EpicDist > MaxDeg:
            InputStream.remove(tr)
    
    return


class Auto_Picking_Initializing():
    def __init__(self, Config_File='PyAP.ini', use_af_db=False):
        
        self.config = configparser.ConfigParser()
        self.config.read(Config_File)
        
        # Database path:
        self.db_assoc_path      = self.config['DataBase']['Assoc_DB']       # associator database
        self.db_nsta_path       = self.config['DataBase']['NSTA_DB']        # nsta24 info database
        #self.db_AF_path         = self.config['DataBase']['AFile_DB']       # A_File database
        self.db_curve_path      = self.config['DataBase']['TT_CURVE_DB']    # Travel-time Curve database
        
        # Our SQLite databases are:
        self.sqlite_assoc       = 'sqlite:///'+self.db_assoc_path           # associator database
        self.sqlite_nsta        = 'sqlite:///'+self.db_nsta_path            # nsta24 info database
        #self.sqlite_AF          = 'sqlite:///'+self.db_AF_path              # A_File database
        self.sqlite_tt_curve    = 'sqlite:///'+self.db_curve_path           # Travel_time Curve database
        
        # Check database exist or not
        if os.path.exists(self.db_assoc_path):
            os.remove(self.db_assoc_path)
        
        if not os.path.exists(self.db_curve_path):
            self.db_tt_curve = self.read_tt_curve()
        else:
            os.remove(self.db_curve_path)
            self.db_tt_curve = self.read_tt_curve()
    
        if not os.path.exists(self.db_nsta_path):
            self.db_nsta = self.read_nsta24_gain()
        else:
            os.remove(self.db_nsta_path)
            self.db_nsta = self.read_nsta24_gain()
        
        # Connect to our databases
        # Associator
        engine_assoc        =   create_engine(self.sqlite_assoc, echo=False)
        tables3D.Base_Assoc.metadata.create_all(engine_assoc)
        Session1=sessionmaker(bind=engine_assoc)
        self.db_assoc       =   Session1()
        
        # AFile database
        if use_af_db:
            engine_AFdb=create_engine(self.sqlite_AF, echo=False)
            Session2=sessionmaker(bind=engine_AFdb)
            self.db_AF          =   Session2()
        
        
    def read_nsta24_gain(self):
        
        engine_nsta=create_engine(self.sqlite_nsta, echo=False)
        table_nsta24.NSTA24.metadata.create_all(engine_nsta)
        Session=sessionmaker(bind=engine_nsta)
        session_nsta=Session()
        
        f = open(self.config['Misc']['NSTA24_dat'],'r')
        for i in f:
            alive = int(i[30])
            if alive == 1:
                sta        = i[0:4].rstrip()
                net        = i[39:44].rstrip()
                loc        = '0'+i[32]
                instrument = i[45:49].rstrip()
                
                longitude  = float(i[5:14])
                latitude   = float(i[13:22])
                elevation  = float(i[22:30])
                
                gain = []
                gain.append(float(i[49:60]))
                gain.append(float(i[60:71]))
                gain.append(float(i[71:82]))
                            
                if instrument == 'FBA':
                    chc = 1
                elif instrument == 'SP':
                    chc = 4
                elif instrument == 'BB':
                    chc = 7
        
                # for each channel in port
                for ch in range(3):
                    #print(num,ch,chc)
                    chn = 'Ch'+str(chc+ch)
                    
                    station = table_nsta24.NSTATable(sta,chn,net,loc,gain[ch],latitude,longitude,elevation)
                    session_nsta.add(station)
        session_nsta.commit()
        f.close()    
        
        return session_nsta
        
    def read_tt_curve(self):
        # create engine
        engine_tt_curve=create_engine(self.sqlite_tt_curve, echo=False)
        # Create the tables for travel-time curve
        table_tt_curve.BASE_TT_CURVE.metadata.create_all(engine_tt_curve)
        Session=sessionmaker(bind=engine_tt_curve)
        session_tt_curve=Session()
        
        
        f = open(self.config['Misc']['TT_CURVE_DAT'],'r')
        for line in f:
            STN     = line[0:4].rstrip()
            A_VALUE = float(line[4:12])
            B_VALUE = float(line[12:20])
            
            stn_table = table_tt_curve.TT_CURVE(STN,A_VALUE,B_VALUE)
            session_tt_curve.add(stn_table)
        f.close()
        session_tt_curve.commit()
        
        return session_tt_curve
        
    def Begin_Small_Events(self):
        # remove old assoc db
        os.remove(self.db_assoc_path)
        
        # re-create new assoc db
        # Associator
        engine_assoc        =   create_engine(self.sqlite_assoc, echo=False)
        tables3D.Base_Assoc.metadata.create_all(engine_assoc)
        Session1=sessionmaker(bind=engine_assoc)
        self.db_assoc       =   Session1()
        
        return 
        
def run_rehypo(tpfilename, method='3D', config=None):
    copyfile(tpfilename, 'hypo.tmp')
    
    tmp_file = open('hypo.tmp', 'r')
    tp = open(tpfilename, 'w')
    
    f_line = True
    for line in tmp_file:
        if f_line:
            # set initial depth for different method
            if method == '1D':
                depth = 20.0
            elif method == '3D':
                depth = float(line[34:40])
            
            tp.write("%s%6.2f%s"%(line[0:34],depth,line[40:]))
            f_line = False
        else:
            tp.write("%s"%(line))
    tp.close()
    tmp_file.close()
    
    if method == '1D':
        os.system(config['Misc']['HypoRun_Path']+' -m 1 -p '+tpfilename)
    elif method == '3D':
        os.system(config['Misc']['Hypo3D_Path']+' '+tpfilename)
    
    
    return
        

class OBS_rotate():
    def __init__(self, st, Config_File = None):
        # get config information
        self.config                         = Config_File
        if not Config_File:
            sys.exit("Config_File not given ---- class OBS_rotate")
        
        # find EOS station in time period
        self.st_time = st[0].stats.starttime
        
        # get station info from file
        self.get_station_info()
        
        # search if EOS in this A file
        EOS_TMP = Stream()
        for STN in self.STNs:
            EOS_TMP += st.select(station=STN)
        
        #print(len(EOS_TMP))
        if len(EOS_TMP) > 0:
            self.rotate_EOS_station(st)
        
        
    def rotate_EOS_station(self, st):
        CHAN_PAIR = (['Ch1','Ch2','Ch3'], ['Ch4','Ch5','Ch6'])
        
        for STN in self.STNs:
            for PAIR in CHAN_PAIR:
                EOS_st = Stream()
                for CHN in PAIR:
                    EOS_st += st.select(station=STN).select(channel=CHN)
                
                if len(EOS_st) > 0:
                    #print(STN, self.STNs_Roll[STN], self.STNs_Pitch[STN], self.STNs_Az[STN])
                    PHI = np.deg2rad(self.STNs_Roll[STN])
                    PSI = np.deg2rad(self.STNs_Pitch[STN])
                    azimuth = self.STNs_Az[STN] - 90
                    self.rotate_seis(EOS_st, PHI, PSI, azimuth)
    
    def make_YAW_PITCH_ROLL(self, x_vector, y_vector, z_vector):
    
        YZ_len  = LA.norm(np.array([y_vector, z_vector]))
        
        PHI = np.arctan2(y_vector, z_vector)
        PSI = np.arctan2(-x_vector, YZ_len)
        
        ROLL  = np.rad2deg(PHI)
        PITCH = np.rad2deg(PSI)
        
        return PHI, PSI
    
    
    def make_rotate_matrix(self, PHI, PSI, azimuth):
    
        COS_PHI = np.cos(-PHI)
        SIN_PHI = np.sin(-PHI)
        
        ROTATEM_PHI = np.array([[1.0, 0.0, 0.0], \
                            [0.0, COS_PHI,  SIN_PHI], \
                            [0.0, -SIN_PHI,  COS_PHI]])
        
        COS_PSI = np.cos(-PSI)
        SIN_PSI = np.sin(-PSI)
        
        ROTATEM_PSI = np.array([[ COS_PSI, 0.0, -SIN_PSI], \
                                [ 0.0, 1.0, 0.0], \
                                [ SIN_PSI, 0.0, COS_PSI]])
        
        
        COS_AZI = np.cos(np.deg2rad(azimuth))
        SIN_AZI = np.sin(np.deg2rad(azimuth))
        
        ROTATEM_azimuth = np.array([[ COS_AZI, SIN_AZI, 0.0], \
                                    [-SIN_AZI, COS_AZI, 0.0], \
                                    [     0.0,     0.0, 1.0]])
        
        ROTATE_M_TMP = np.matmul(ROTATEM_PSI, ROTATEM_PHI)
        ROTATE_M = np.matmul(ROTATEM_azimuth, ROTATE_M_TMP)
        
        return ROTATE_M
    
    
    def rotate_seis(self, st, PHI, PSI, azimuth):
    
        #new_st = st.copy()
        
        tmp = np.append([st[2].data.data, st[1].data.data], [st[0].data.data], axis=0)
        SEIS_3D = tmp.reshape(3, st[2].stats.npts)
        ROTM_3D = self.make_rotate_matrix(PHI, PSI, azimuth)
    
        SEIS_NEW = np.matmul(ROTM_3D, SEIS_3D)
    
        st[2].data = SEIS_NEW[0]
        st[1].data = SEIS_NEW[1]
        st[0].data = SEIS_NEW[2]
        
        return
        
    def get_station_info(self):
        ORI_DAT = open(self.config['Misc']['OBS_ORI_DAT'], 'r')
        TIME = self.st_time.year*10000+self.st_time.month*100+self.st_time.day
        self.STNs = []
        self.STNs_Roll   = {}
        self.STNs_Pitch  = {}
        self.STNs_Az     = {}
        
        for line in ORI_DAT:
            STN = line[0:5].strip()
            period_start = int(line[37:45])
            period_stop  = int(line[46:54])
            
            if STN not in self.STNs:
                if TIME >= period_start and TIME <= period_stop:
                    self.STNs.append(STN)
                    self.STNs_Roll[STN]  = float(line[8:18])
                    self.STNs_Pitch[STN] = float(line[19:29])
                    self.STNs_Az[STN]    = float(line[30:36])
        ORI_DAT.close()
        return

class Pfile():
    def __init__(self, stmag, st_inten, db_assoc, db_nsta, Config_File = None, polarity_mode='False', Vector_Mode=False):
    
        # Define travel time and associator database
        self.config                         = Config_File
        if not Config_File:
            sys.exit("Config_File not given ---- class Pfile")
        
        
        self.Vector_Mode                    = Vector_Mode
        
        self.st                             = stmag
        self.st_inten                       = st_inten
        self.assoc_db                       = db_assoc 
        self.nsta_db                        = db_nsta
        if polarity_mode == 'True':
            self.polarity_mode              = True
        else:
            self.polarity_mode              = False
        
        self.paz_wa_v = {'sensitivity': 2800, 'zeros': [0j], 'gain': 1, 'poles': [-6.2832 - 4.7124j, -6.2832 + 4.7124j]} # for Velocity
        self.paz_wa_a = {'sensitivity': 2800, 'zeros': [], 'gain': 1, 'poles': [-6.2832 - 4.7124j, -6.2832 + 4.7124j]}   # for Accel
        #self.paz_wa_v = {'sensitivity': 2080, 'zeros': [0j], 'gain': 1, 'poles': [-5.49779 - 5.60886j, -5.49779 + 5.60886j]} # IASPEI 2012, for Velocity
        #self.paz_wa_a = {'sensitivity': 2080, 'zeros': [], 'gain': 1, 'poles': [-5.49779 - 5.60886j, -5.49779 + 5.60886j]}   # IASPEI 2012, for Accel
        
    def get_assoc_ids(self):
        ids = self.assoc_db.query(Associated).filter(Associated.id < 99).all()
    
        return ids