from sqlalchemy import *
from sqlalchemy.ext.declarative import declarative_base

NSTA24=declarative_base()

class NSTATable(NSTA24):
    __tablename__='NSTATable'
    id=Column(Integer, primary_key=True)
    sta=Column(String(5))
    net=Column(String(4))
    loc=Column(String(2))
    cha=Column(String(3))
    gain=Column(Float)
    latitude=Column(Float)
    longitude=Column(Float)
    elevation=Column(Float)
    starttime=Column(DateTime)
    endtime=Column(DateTime)
    
    def __init__(self,sta,cha,net,loc,gain,latitude,longitude,elevation):
        self.sta=sta
        self.net=net
        self.loc=loc
        self.cha=cha
        self.gain=gain
        self.latitude=latitude
        self.longitude=longitude
        self.elevation=elevation
        self.starttime=None
        self.endtime=None
    
    def __repr__(self):
        return "NSTA <%s.%s.%s.%s %s %s %s>" % (self.sta,self.cha,self.net,self.loc,self.gain, self.latitude,self.longitude)