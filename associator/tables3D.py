from sqlalchemy import *   # from sqlalchemy import create_engine, Column, Integer, String, float, Boolean
from sqlalchemy.ext.declarative import declarative_base

# Declare mapping
Base_Assoc=declarative_base()

class Pick(Base_Assoc):
  __tablename__="picks"
  id=Column(Integer,primary_key=True)
  sta=Column(String(5))
  chan=Column(String(3))
  net=Column(String(2))
  loc=Column(String(2))
  time=Column(DateTime)
  snr=Column(Float)
  phase=Column(String(1))
  uncert=Column(Float)
  polarity=Column(String(1))
  locate_flag=Column(Boolean)
  assoc_id=Column(Integer)
  modified_id=Column(Integer)
  t_create=Column(DateTime)
  
  def __init__(self,scnl,picktime,polarity,snr,uncert,t_create,phase=None):
    self.sta=scnl.station
    self.chan=scnl.channel
    self.net=scnl.network
    self.loc=scnl.location
    self.time=picktime
    self.snr=snr
    self.uncert=uncert
    self.modified_id=None
    self.phase=phase
    self.polarity=polarity
    self.locate_flag=None
    self.assoc_id=None
    self.t_create=t_create
    
  def __repr__(self):
    return "Pick <%s.%s.%s.%s %s %s %s %s>" % (self.sta,self.chan,self.net,self.loc,self.time.isoformat("T"),self.phase,self.modified_id,self.assoc_id)


class PickModified(Base_Assoc):
  __tablename__="picks_modified"
  id=Column(Integer,primary_key=True)
  sta=Column(String(5))
  chan=Column(String(3))
  net=Column(String(2))
  loc=Column(String(2))
  time=Column(DateTime)
  phase=Column(String(1))
  error=Column(Float)
  locate_flag=Column(Boolean)
  assoc_id=Column(Integer)
  
  def __init__(self,sta,chan,net,loc,picktime,phase,error,assoc_id):
    self.sta=sta
    self.chan=chan
    self.net=net
    self.loc=loc
    self.time=picktime
    self.phase=phase
    self.error=error
    self.locate_flag=None
    self.assoc_id=assoc_id
    
  def __repr__(self):
    return "Pick <%s.%s.%s.%s %s %s %s %s>" % (self.sta,self.chan,self.net,self.loc,self.time.isoformat("T"),self.phase,self.error,self.assoc_id)
  

class Candidate(Base_Assoc):
  __tablename__="candidate"
  id=Column(Integer,primary_key=True)
  ot=Column(DateTime)
  sta=Column(String(5))
  #d_km=Column(Float)   # dont need this anymore
  #delta=Column(Float)  # dont need this anymore
  weight=Column(Float)
  # P and S travel times are not completely necessary because they can be calculated but simpler to save
  tp=Column(DateTime)
  p_modified_id=Column(Integer)  # modified pick ID
  ts=Column(DateTime)
  s_modified_id=Column(Integer)  # modified pick ID
  locate_flag=Column(Boolean)
  assoc_id=Column(Integer)
  
  def __init__(self, ot, sta, tp, p_modified_id, ts, s_modified_id):
    self.ot=ot
    self.sta=sta
    #self.d_km=d_km
    #self.delta=delta
    self.weight=None
    self.tp=tp
    self.ts=ts
    self.p_modified_id=p_modified_id
    self.s_modified_id=s_modified_id
    self.locate_flag=None
    self.assoc_id=None
  
  def __repr__(self):
    return "Candidate Event <%s %s %d %d>" % (self.ot.isoformat("T"),self.sta,self.p_modified_id, self.s_modified_id)

  #def __str__(self):
  #  return "Candidate Event <%s %s %d %d>" % (self.ot.isoformat("T"),self.sta,self.p_modified_id, self.s_modified_id)

  def set_assoc_id(self,assoc_id,session,FT):
    self.assoc_id=assoc_id
    self.locate_flag=FT
    # Assign phases to modified picks
    
    # Actually only one pick_p and pick_s
    pick_p=session.query(PickModified).filter(PickModified.id==self.p_modified_id)
    for pick in pick_p:
      pick.phase='P'
      pick.assoc_id=assoc_id
      pick.locate_flag=FT
    
    pick_s=session.query(PickModified).filter(PickModified.id==self.s_modified_id)
    for pick in pick_s:
      pick.phase='S'
      pick.assoc_id=assoc_id
      pick.locate_flag=FT
    
    # Assign the phases to picks contribute to a modified picks
    picks_p=session.query(Pick).filter(Pick.modified_id==self.p_modified_id).all()
    for pick in picks_p:
      pick.phase='P'
      pick.assoc_id=assoc_id
      pick.locate_flag=FT
    
    picks_s=session.query(Pick).filter(Pick.modified_id==self.s_modified_id).all()
    for pick in picks_s:
      pick.phase='S'
      pick.assoc_id=assoc_id
      pick.locate_flag=FT

    
class Associated(Base_Assoc):
  __tablename__="associated"
  id=Column(Integer,primary_key=True)
  ot=Column(DateTime)
  ot_uncert=Column(Float)
  latitude=Column(Float)
  longitude=Column(Float)
  depth=Column(Float)
  rms=Column(Float)
  nsta=Column(Integer)
  sourcegrid=Column(Integer)
  t_create=Column(DateTime)
  t_update=Column(DateTime)                          
  
  def __init__(self,ot,ot_uncert,latitude,longitude,depth,rms,nsta,sourcegrid,t_create,t_update):
    self.ot=ot
    self.ot_uncert=ot_uncert
    self.latitude=latitude
    self.longitude=longitude
    self.depth=depth
    self.rms=rms
    self.nsta=nsta
    self.sourcegrid=sourcegrid
    self.t_create=t_create
    self.t_update=t_update
    
  def __repr__(self):
    return "Associated Event <%s %s %.3f %.3f %d %d>" % (self.ot.isoformat("T"),self.ot,self.latitude,self.longitude,self.nsta,self.sourcegrid)
