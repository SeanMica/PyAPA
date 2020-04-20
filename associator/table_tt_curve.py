from sqlalchemy import *
from sqlalchemy.ext.declarative import declarative_base

BASE_TT_CURVE=declarative_base()

class TT_CURVE(BASE_TT_CURVE):
    __tablename__='TT_C_Table'
    id      =   Column(Integer, primary_key=True)
    sta     =   Column(String(5))
    # y = ax+b , where x = s_p time, y = p travel time
    a_value =   Column(Float)
    b_value =   Column(Float)
    
    def __init__(self, sta, a_value, b_value):
        self.sta        =   sta
        self.a_value    =   a_value
        self.b_value    =   b_value
    
    def __repr__(self):
        return "TT_CURVE <%s.%s.%s>" % (self.sta, self.a_value, self.b_value)