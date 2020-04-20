import numpy as np
from sqlalchemy.orm import *  # session
from sqlalchemy import create_engine
from obspy.core import UTCDateTime
from obspy.geodetics import locations2degrees, degrees2kilometers
from .table_tt_curve import *
from .table_nsta24 import *
from .tables3D import *
from .search_hypo import *
#from .pbr import pbr
from datetime import *
from operator import itemgetter
from itertools import combinations
import logging
import time       #"from datetime import *" will import time, name space will be overwritten

class LocalAssociator():
  """
  The 3D Associator associate picks with travel time curve of 3D velocity.
  """
  def __init__(self, assoc_db, nsta_db, tt_curve_db, fileHeader, max_Parr = 80, aggregation = 1, aggr_norm = 'L2', assoc_ot_uncert = 3, nsta_declare = 6, config_par = None, AS_MODE = False):
    """
    Parameters:
    db_assoc: associator database
    db_tt: travel time table database
    max_km: maximum distance of S-P interval in distance
    aggregation: the coefficient multiplied to minimum travel time
    aggr_norm: L2: median; L1: mean
    assoc_ot_uncert: origin time uncertainty window
    nsta_declare: minimum station number to declare a earthquake
    nt, np, nr: node geometry
    """
    
    self.assoc_db                 =         assoc_db
    self.nsta_db                  =         nsta_db
    self.tt_curve_db              =         tt_curve_db
    self.max_Parr                 =         max_Parr
    self.max_s_p                  =         max_Parr*0.75       # need tuning
    self.min_s_p                  =         1.0
    self.aggregation              =         aggregation
    self.aggr_window              =         self.aggregation * self.min_s_p
    self.aggr_norm                =         aggr_norm           # L1 takes the mean; L2 takes the median
    self.assoc_ot_uncert          =         assoc_ot_uncert     # Uncertainty of origin times of candidates
    self.nsta_declare             =         nsta_declare        # number observation to declare an evnet
    self.AS_MODE                  =         AS_MODE             # aftershock mode
    self.AS_STNs                  =         []
    
    self.zerotime                 =         UTCDateTime(year=fileHeader[1], month=fileHeader[4], day=fileHeader[5], hour=fileHeader[6], minute=fileHeader[7])
    self.config                   =         config_par
  
  def id_candidate_events(self):
    """ Create a set of possible candidate events from our picks table.
    Where session is the connection to the sqlalchemy database.
    This method simply takes all picks with time differences less than our maximum S-P
    times for each station and generates a list of candidate events.
    """
    now1 = time.time()
    #############
    # Get all stations with unnassoiated picks
    stations=self.assoc_db.query(Pick.sta).filter(Pick.assoc_id==None).distinct().all()
    #print('stations:',len(stations))
    
    if self.AS_MODE :
      print('AS_MODE initial')
      AS_LAT = self.config['AS_MODE'].getfloat('Latitude')
      AS_LON = self.config['AS_MODE'].getfloat('Longitude')
      AS_RAD = self.config['AS_MODE'].getfloat('Radius')
      self.AS_SP_TIME = timedelta(seconds=(AS_RAD/8.0)*1.5)
      self.AS_STNs = self.search_stns_in_range(AS_LAT, AS_LON, AS_RAD, stations)
      print('AS_MODE initial end')
    else:
      #self.AS_STNs = []
      for STN, in stations:
        self.AS_STNs.append(STN)
    
    for sta, in stations:  # the comma is needed
      picks=self.assoc_db.query(Pick).filter(Pick.sta==sta).filter(Pick.assoc_id==None).order_by(Pick.time).all()
      # Condense picktimes that are within our pick uncertainty value picktimes are python datetime objects
      if stations.index((sta,))==0: #stupid tuple
        counter0=0
        picktimes_new,counter=pick_cluster(self.assoc_db,picks,self.aggr_window,self.aggr_norm,counter0)
      else:
        picktimes_new,counter=pick_cluster(self.assoc_db,picks,self.aggr_window,self.aggr_norm,counter)
      
      nets = self.assoc_db.query(PickModified.net).filter(PickModified.sta==sta).filter(PickModified.assoc_id==None).all()
      locs = self.assoc_db.query(PickModified.loc).filter(PickModified.sta==sta).filter(PickModified.assoc_id==None).all()
        
      #for net in nets:
      #if sta in STNs:
      for net, in set(nets):
        for loc, in set(locs):
          picks_modified=self.assoc_db.query(PickModified).filter(PickModified.sta==sta,PickModified.net==net,PickModified.loc==loc).filter(PickModified.assoc_id==None).order_by(PickModified.time).all()
        #picks_modified=self.assoc_db.query(PickModified).filter(PickModified.sta==sta).filter(PickModified.assoc_id==None).order_by(PickModified.time).all()
      
          # Generate all possible candidate events
          for i in range(0, len(picks_modified) - 1):
            for j in range(i + 1,len(picks_modified)):
              s_p = (picks_modified[j].time - picks_modified[i].time).total_seconds()#; print s_p
              if s_p <= self.max_s_p and s_p >= self.min_s_p:
                ot = self.find_ot_from_tt_curve(sta, picks_modified[i], s_p)
                new_candidate=Candidate(ot, sta, picks_modified[i].time, picks_modified[i].id, picks_modified[j].time, picks_modified[j].id)
            
                self.assoc_db.add(new_candidate)
      self.assoc_db.commit()
    
    print('id_candidate time in seconds: ',time.time()-now1)
      
  def associate_candidates(self):
    """ Associate all possible candidate events by comparing the projected origin-times.  At
    this point we are not dealing with the condition that more picks and candidate events 
    could be arriving while we do our event associations.
    """
    now2 = time.time()
    
    dt_ot=timedelta(seconds=self.assoc_ot_uncert)

    # Query all candidate ots
    #candidate_ots=self.assoc_db.query(Candidate).filter(Candidate.assoc_id==None).order_by(Candidate.ot).all()
    if self.AS_MODE :
      candidate_ots=self.assoc_db.query(Candidate).filter(Candidate.assoc_id==None, Candidate.sta.in_(self.AS_STNs)).\
                        filter((Candidate.ts-Candidate.tp) <= self.AS_SP_TIME).order_by(Candidate.ot).all()
    else:
      candidate_ots=self.assoc_db.query(Candidate).filter(Candidate.assoc_id==None, Candidate.sta.in_(self.AS_STNs)).order_by(Candidate.ot).all()
    L_ots=len(candidate_ots) #; print(L_ots)
    Array=[]
    for i in range(L_ots):
      #cluster=self.assoc_db.query(Candidate).filter(Candidate.assoc_id==None).filter(Candidate.ot>=candidate_ots[i].ot, Candidate.ot<(candidate_ots[i].ot+dt_ot)).order_by(Candidate.ot).all()
      if self.AS_MODE :
        cluster=self.assoc_db.query(Candidate).filter(Candidate.assoc_id==None, Candidate.sta.in_(self.AS_STNs)).\
                          filter((Candidate.ts-Candidate.tp) <= self.AS_SP_TIME).\
                          filter(Candidate.ot>=candidate_ots[i].ot, Candidate.ot<(candidate_ots[i].ot+dt_ot)).order_by(Candidate.ot).all()
      else:
        cluster=self.assoc_db.query(Candidate).filter(Candidate.assoc_id==None, Candidate.sta.in_(self.AS_STNs)).\
                          filter(Candidate.ot>=candidate_ots[i].ot, Candidate.ot<(candidate_ots[i].ot+dt_ot)).order_by(Candidate.ot).all()
      #cluster_sta=self.assoc_db.query(Candidate.sta).filter(Candidate.assoc_id==None).filter(Candidate.ot>=candidate_ots[i].ot).filter(Candidate.ot<(candidate_ots[i].ot+dt_ot)).order_by(Candidate.ot).all()
      cluster_sta = [candi.sta for candi in cluster]
      #print(cluster_sta)
      l_cluster=len(set(cluster_sta))
      Array.append((i,l_cluster,len(cluster)))
    #print Array
    Array.sort(key=itemgetter(1), reverse=True)  #sort Array by l_cluster, notice Array has been changed
    #print Array
    
    print('candidates_ots:', time.time()-now2, ', Array length:', len(Array))
    
    if not self.AS_MODE:
        Array_count = len(Array)
        if Array_count > 1500:
            print('VERY VERY HUGE EARTHQUAKE MODE!')
            dt_ot=timedelta(seconds=self.assoc_ot_uncert*2)
            self.nsta_declare = 25
        elif Array_count > 800:
            print('HUGE EARTHQUAKE MODE!')
            dt_ot=timedelta(seconds=self.assoc_ot_uncert*2)
            self.nsta_declare = 15
        elif Array_count > 500:
            print('BIG EARTHQUAKE MODE!')
            dt_ot=timedelta(seconds=self.assoc_ot_uncert*2)
            self.nsta_declare = 10
        elif Array_count > 400:
            print('Medium EARTHQUAKE MODE!')
            dt_ot=timedelta(seconds=self.assoc_ot_uncert*2)
            self.nsta_declare = 8


    for i in range(len(Array)):
      index=Array[i][0]
      if Array[i][1]>=self.nsta_declare:
        matches=self.assoc_db.query(Candidate).filter(Candidate.assoc_id==None).\
                    filter(Candidate.ot>=candidate_ots[index].ot).\
                    filter(Candidate.ot<(candidate_ots[index].ot+dt_ot)).order_by(Candidate.ot).all() 
        
        
        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # remove the candidates with the modified picks has been associated
        picks_associated_id=list(set(self.assoc_db.query(PickModified.id).filter(PickModified.assoc_id!=None).all()))
        index_matches=[]
        for id, in picks_associated_id:
          for j,match in enumerate(matches):
            if match.p_modified_id==id or match.s_modified_id==id:
              index_matches.append(j)        
        # delete from the end
        if index_matches:
          for j in sorted(set(index_matches),reverse=True):
            del matches[j]
        # remove the candidates with the modified picks has been associated
        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        #print(i)
        
        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++        
        # 3D Associator
        now = time.time()
        tt = []
        #print(len(matches))
        for match in matches:
          #print('sta:',match.sta)
          match_p = match.tp
          match_s = match.ts
          match_ot = match.ot
          match_ttp = (match_p - match_ot).total_seconds()
          match_tts = (match_s - match_ot).total_seconds()
          tt.append([match.sta, match_ttp, match_tts, match.ot, match])
        
        cb, cb_dupl = self.new_remove_comb(tt)
        
        rms_sort = []
        tt_cb = cb[0]
        if len(tt_cb) >= self.nsta_declare: # self.nsta_declare has to be greater than or equal to 3
          ##PyramidSearching##
          #print(tt_cb)
          tt_new, evt_lat, evt_lon, evt_dep, QA = HYPO3D_Searching(self.assoc_db, self.nsta_db, tt_cb, self.zerotime, config=self.config)
          #rms_sort.append((tt_new, sourcegrid, rms, 1))
          print('HYPO3D_Searching done.', QA)
          
          #if QA:
          #  dp_cb = cb_dupl[0]
          #  cb = self.add_dupl_ots_new(tt_new, dp_cb, evt_lat, evt_lon, evt_dep, 2.0)
          #  
          #  tt_cb = cb
          
          if len(tt_cb) >= self.nsta_declare and QA: # self.nsta_declare has to be greater than or equal to 3
          
            rms = 0.0
            sourcegrid=35000 # useless thing, lazy to remove
            lat, lon, dep = evt_lat, evt_lon, evt_dep
            
            nsta = len(tt_new)
            
            all_ots = []
            for j in range(nsta):  
              all_ots.append(tt_new[j][3])
            origintime, ot_unc = datetime_statistics(all_ots)
            # in 3D Associator, use rms of picks instead of loc_uncertainty
            t_create = datetime.utcnow()
            t_update = datetime.utcnow()
            new_event=Associated(origintime,round(ot_unc,3),lat,lon,dep,round(rms,3),nsta,sourcegrid,t_create,t_update)    
            self.assoc_db.add(new_event)
            self.assoc_db.flush()
            self.assoc_db.refresh(new_event)
            self.assoc_db.commit()
            event_id=new_event.id
            logging.info(str(['event_id:',event_id])) 
            logging.info(str(['ot:', origintime, 'ot_uncert:', round(ot_unc,3), 'loc:', lat,lon,dep, 'rms:', round(rms,3), 'nsta:', nsta]))
            
            #print(event_id)
            for tt_tuple in cb[0]:#[index]:
              match = tt_tuple[4]
              #print(match,match.assoc_id)
              match.set_assoc_id(event_id,self.assoc_db,True)
            self.assoc_db.commit()
            
            # remove all picks near this time
            #not_ots=self.assoc_db.query(Candidate).filter(Candidate.assoc_id==None).filter(Candidate.ot>=(origintime-dt_ot)).filter(Candidate.ot<(origintime+dt_ot)).all()
            mask_time = timedelta(seconds=10)
            not_ots=self.assoc_db.query(Candidate).filter(Candidate.assoc_id==None).filter(Candidate.ot>=(origintime-mask_time)).filter(Candidate.ot<(origintime+mask_time)).all()
            for not_tt in not_ots:
              not_tt.set_assoc_id(99,self.assoc_db,True)
            self.assoc_db.commit()

      else:
        break
  
  
  # remove stations with dupl ots from combinations
  def new_remove_comb(self,tt):
    stns = [item[0] for item in tt]
    dupl_stns = list_duplicates(stns)
    
    f_stns = [stn for stn in dupl_stns if stns.count(stn) < 3]
    
    not_dupl_tt = [item for item in tt if item[0] not in dupl_stns]
    dupl_tt = [item for item in tt if item[0] in f_stns]
    
    cb = []
    cb_dupl = []
    cb.append(tuple(not_dupl_tt))
    cb_dupl.append(tuple(dupl_tt))
    
    # only return combinations of different stations
    return cb, cb_dupl

  def new_remove_comb_2(self,tt):
    stns = [item[0] for item in tt]
    dupl_stns = list_duplicates(stns)
    
    not_dupl_tt = [item for item in tt if item[0] not in dupl_stns]
    dupl_tt = [item for item in tt if item[0] in dupl_stns]
    
    cb = []
    cb_tmp = []
    stn_done = []
    
    for dp_tt in dupl_tt:
        if dp_tt[0] not in stn_done:
            cb_tmp.append(dp_tt)
            stn_done.append(dp_tt[0])
    
    cb.append(tuple(cb_tmp)+tuple(not_dupl_tt))
    
    # only return combinations of different stations
    return cb
    
  def find_ot_from_tt_curve(self, sta, p_arr, s_p_time):
    check_db = self.tt_curve_db.query(TT_CURVE.id).filter(TT_CURVE.sta == sta).first()
    if check_db == None:
        a_value = 1.35 # some default value
        b_value = 0.0
    else:
        a_value, b_value = self.tt_curve_db.query(TT_CURVE.a_value, TT_CURVE.b_value).filter(TT_CURVE.sta == sta).first()
    
    ot = p_arr.time - timedelta(seconds=s_p_time*a_value + b_value)
    
    return ot
   
  def add_dupl_ots_new(self, tt, tt_dupl, evt_lat, evt_lon, evt_dep, tt_residual):
    # find stns
    tt_new = tt
    stn_done = []
    for tt_dupls in tt_dupl:
      if tt_dupls[0] not in stn_done:
        sta = tt_dupls[0]
        ttp = tt_dupls[1]
        tts = tt_dupls[2]
        ttsp = tts-ttp
        
        # check sta in TTtable3D
        sta_id, = self.nsta_db.query(NSTATable.id).filter(NSTATable.sta == sta).first()
        sta_lat, sta_lon, sta_dep = self.nsta_db.query(NSTATable.latitude,NSTATable.longitude,NSTATable.elevation).filter(NSTATable.id == sta_id).first()
        sta_dep = sta_dep * (-0.001)
        
        P_ttime = pbr(evt_lat,evt_lon,evt_dep,sta_lat,sta_lon,sta_dep,1)
        S_ttime = pbr(evt_lat,evt_lon,evt_dep,sta_lat,sta_lon,sta_dep,2)
        S_P_ttime = S_ttime-P_ttime
        
        if abs(P_ttime-ttp) <= tt_residual and abs(S_ttime-tts) <= tt_residual and abs(S_P_ttime-ttsp) <= tt_residual:
          tt_new.append(tt_dupls)
          stn_done.append(sta)
    
    return tt_new
    
  def search_stns_in_range(self, lat, lon, rad, stations):
    
    stns = []
    for stn, in stations:
      stn_lon, stn_lat = self.nsta_db.query(NSTATable.longitude, NSTATable.latitude).filter(NSTATable.sta==stn).first()
      
      Dist=degrees2kilometers(locations2degrees(lat, lon, stn_lat, stn_lon))
      
      if Dist <= rad:
        stns.append(stn)
      
    return stns

  
def list_duplicates(seq):
  seen = set()
  seen_add = seen.add
  # adds all elements it doesn't know yet to seen and all other to seen_twice
  seen_twice = set( x for x in seq if x in seen or seen_add(x) )
  # turn the set into a list (as requested)
  return list( seen_twice )



def datetime_statistics(dt_list,norm='L2'):
  """ mean,std=datetime_statistics(datetime_list)
  Calculate the mean and standard deviations in seconds of a list of datetime values
  """
  offsets=[]
  for dt in dt_list:
    offsets.append((dt-dt_list[0]).total_seconds())
  if norm=='L1':
    mean_offsets=np.mean(offsets)
    new_pick = dt_list[0]+timedelta(seconds=mean_offsets)
  elif norm=='L2':
    mean_offsets=np.median(offsets)
    new_pick = dt_list[0]+timedelta(seconds=mean_offsets)
  elif norm=='FA':
    mean_offsets=np.median(offsets)
    tmp_pick = dt_list.copy()
    tmp_pick.sort()
    new_pick = tmp_pick[0]
  std_offsets=np.std(offsets)
  return new_pick,std_offsets




  
def pick_cluster(session,picks,pickwindow,aggr_norm,counter):
  """ cleaning up very closed picks on different channels of same station
  """
  #                     |    |                     /\
  #                     |    |                    /  \          /\
  #                     |    | /\      /\        /    \        /  \      /\
  #        _____________|/\__|/  \    /  \      /      \      /    \    /  \  /\_________
  #                     |    |    \  /    \    /        \    /      \  /    \/
  #                     |    |     \/      \  /          \  /        \/
  #                     |    |              \/            \/

  # pickwindow:          ----                                      better to set pickwindow==t_up, t_up is to clean closed picks
  # STA1 E   -----------|----|--------------------|--------------
  # STA1 N   ------------|-------------------------|-------------
  # STA1 Z   -------------|-------------------------|------------
  # stack    -----------|||--|--------------------|||------------  
  # cluster STA1 --------|---|---------------------|-------------  chen highly recommend to use norm=='L2' to lower the effect of outlier, L2 takes median
  # ARGUE: whether only take the median or mean of the picks from different stations? won't count the followings after first one
  # 
  
  picks_new=[]
  # only one pick in picks
  if len(picks)==1:
    cluster=[];cluster.append(picks[0]);cluster_time=[];cluster_time.append(picks[0].time)
    picks[0].modified_id=1+counter # assign modified id to picks
    counter+=1
    pickave,pickstd=datetime_statistics(cluster_time,aggr_norm)
    # append the row to the picks_new, not only the pick time
    picks_new.append(picks[0])
    pick_modified=PickModified(picks[0].sta,picks[0].chan,picks[0].net,picks[0].loc,picks[0].time,picks[0].phase,round(pickstd,3),picks[0].assoc_id)
    session.add(pick_modified)
    session.commit()
    
  # more than one pick in picks
  else:
    j=0
    counter=1+counter
    while True:
      i=j
      cluster=[];cluster.append(picks[i]);cluster_time=[];cluster_time.append(picks[i].time);channel=[];channel.append(picks[i].chan)
      picks[i].modified_id=counter
      while True:
        # cluster picks of different channels; notice that for the closed picks on the same station, those picks behind the first pick could be separated lonely or separated cluster
        if picks[i+1].chan not in channel and (picks[i+1].time-picks[i].time).total_seconds()<pickwindow:
          cluster.append(picks[i+1])
          cluster_time.append(picks[i+1].time)
          channel.append(picks[i+1].chan)
          picks[i+1].modified_id=counter # assign modified id to picks     
          i=i+1
          # make sure do not go over the range limit because j=i+1 below, jump out inner while loop
          if i==len(picks)-1:
            break
        # elif is dealing with the exactly same picks, probably from processing same stream twice
        elif picks[i+1].sta==picks[i].sta and picks[i+1].chan==picks[i].chan and picks[i+1].time==picks[i].time: # and picks[i+1].snr==picks[i].snr and picks[i+1].phase==picks[i].phase and picks[i+1].uncert==picks[i].uncert:
          cluster.append(picks[i+1])
          cluster_time.append(picks[i+1].time)
          channel.append(picks[i+1].chan)
          picks[i+1].modified_id=counter # assign modified id to picks     
          i=i+1
          # make sure do not go over the range limit because j=i+1 below, jump out inner while loop
          if i==len(picks)-1:
            break
        else:
          break
      pickave,pickstd=datetime_statistics(cluster_time,aggr_norm)
      
      # append whole rows to the picks_new, not only the pick time
      for pick in cluster:
        if aggr_norm == 'FA' and (pick.time-pickave).total_seconds()==0:
          break
        if (pick.time-pickave).total_seconds()>=0:
          break
      picks_new.append(pick)
      pick_modified=PickModified(pick.sta,pick.chan,pick.net,pick.loc,pick.time,pick.phase,round(pickstd,3),pick.assoc_id)
      session.add(pick_modified)
      session.commit()
      # next cluster
      j=i+1
      counter=counter+1
      
      # jump outer while loop and compare last two picks. For the situation that last one is ungrouped, use if statement to add in picks_new
      if j>=len(picks)-1:
        if (picks[-1].time-picks[-2].time).total_seconds()>pickwindow:
          picks_new.append(picks[-1])
          picks[-1].modified_id=counter # assign modified id to picks
          pick_modified=PickModified(picks[-1].sta,picks[-1].chan,picks[-1].net,picks[-1].loc,picks[-1].time,picks[-1].phase,round(pickstd,3),picks[-1].assoc_id)
          session.add(pick_modified)
          session.commit()
        else:
          if picks[-1] in cluster:
            counter-=1
          else:
            picks[-1].modified_id=counter
            pick_modified=PickModified(picks[-1].sta,picks[-1].chan,picks[-1].net,picks[-1].loc,picks[-1].time,picks[-1].phase,round(pickstd,3),picks[-1].assoc_id)
            session.add(pick_modified)
            session.commit()
        break     
  
  return picks_new, counter
  
