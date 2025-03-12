import numpy as np
from sgp4.api import Satrec
from tools.date_trans import calculate_jday
from datetime import datetime
from tools.sgp4 import *

def initialize_time_vector(ConjStartDate: datetime, ConjEndDate: datetime, PropTimeStep: float):
    totalelaps = (ConjEndDate - ConjStartDate).total_seconds() / 60
    timeVec = np.arange(0, np.ceil(totalelaps / PropTimeStep) * PropTimeStep + PropTimeStep, PropTimeStep)
    if timeVec[-1] > totalelaps:
        timeVec[-1] = totalelaps
    return timeVec

class ObjSatDetail:
    def __init__(self, objNum):
        self.satobj = [dict() for _ in range(objNum)]
        self.objpnow = np.zeros((objNum, 3))
        self.objvnow = np.zeros((objNum, 3))
        self.CurrentOr = np.zeros((objNum, 3))
        self.CurrentOt = np.zeros((objNum, 3))
        self.CurrentOh = np.zeros((objNum, 3))

    def calculate_objSat_detail(self, ObjSat, objNum, ConjStartJulian):
        for ii in range(objNum):
            # 解析 TLE 数据
            self.satobj[ii]['struc'] = {'satnum': ObjSat[ii]['CatID']}
            self.satobj[ii]['sattle'] = Satrec.twoline2rv(ObjSat[ii]['Line2'], ObjSat[ii]['Line3'])

            # 计算初始时间和偏移
            self.satobj[ii]['initialepoch'] = self.satobj[ii]['sattle'].jdsatepoch - calculate_jday(1950, 1, 0, 0, 0, 0)
            self.satobj[ii]['initialjulian'] = self.satobj[ii]['sattle'].jdsatepoch
            self.satobj[ii]['offset'] = (ConjStartJulian - self.satobj[ii]['initialjulian']) * 1440 

            # 初始化 SGP4 轨道参数
            self.satobj[ii]['struc'] = self.satobj[ii]['sattle']

            # 计算卫星的位置和速度
            _, p, v = self.satobj[ii]['struc'].sgp4(self.satobj[ii]['initialjulian'], self.satobj[ii]['offset'] / 1440.0)
            self.objpnow[ii, :] = p
            self.objvnow[ii, :] = v

            # 计算径向、沿轨、交叉轨向量
            or_ = p / np.linalg.norm(p)
            h = np.cross(p, v)
            oh = h / np.linalg.norm(h)
            ot = np.cross(oh, or_)
            ot = ot / np.linalg.norm(ot)
            self.CurrentOr[ii, :] = or_
            self.CurrentOt[ii, :] = ot
            self.CurrentOh[ii, :] = oh

    def compute_obj_relative_values(self, objNum):
        if objNum > 1:
            objp = self.objpnow.T 
            objv = self.objvnow.T
            objtemprange = []
            objtemprangerate = []
            for obji in range(objNum - 1):
                # 计算相对位置和速度
                objdp = objp[:, obji + 1:] - np.tile(objp[:, obji].reshape(-1, 1), (1, objNum - obji - 1))
                objdv = objv[:, obji + 1:] - np.tile(objv[:, obji].reshape(-1, 1), (1, objNum - obji - 1))
                # 计算相对距离
                objtemprange.extend(np.linalg.norm(objdp, axis=0))
                # 计算相对距离变化率
                objtemprangerate.extend(np.sum(objdp * objdv, axis=0))
            objCurrentRange = np.array(objtemprange)
            objCurrentRangeRate = np.sign(np.array(objtemprangerate) / objCurrentRange)
            return objCurrentRange, objCurrentRangeRate
        else:
            return None, None
        

class TgtSatDetail:
    def __init__(self, tgtNum, objNum):     
        self.sattgt = [dict() for _ in range(tgtNum)]  
        self.tgtpnow = np.zeros((tgtNum, 3))  
        self.tgtvnow = np.zeros((tgtNum, 3))  
        self.RelativePx = np.zeros((tgtNum, objNum))  
        self.RelativePy = np.zeros((tgtNum, objNum)) 
        self.RelativePz = np.zeros((tgtNum, objNum))  
        self.RelativeVx = np.zeros((tgtNum, objNum))
        self.RelativeVy = np.zeros((tgtNum, objNum))  
        self.RelativeVz = np.zeros((tgtNum, objNum))  
        self.CurrentRangeRate = np.zeros((tgtNum, objNum))  
        self.CurrentRange = np.zeros((tgtNum, objNum))

    def calculate_tgtSat_detail(self, TgtSat, objNum, ConjStartJulian, objpnow, objvnow):
        for ii in range(len(TgtSat)):
            # 解析 TLE 数据
            self.sattgt[ii]['struc'] = {'satnum': TgtSat[ii]['CatID']}
            self.sattgt[ii]['Name'] = TgtSat[ii]['Name']
            self.sattgt[ii]['sattle'] = Satrec.twoline2rv(TgtSat[ii]['Line2'], TgtSat[ii]['Line3'])

            # 计算初始时间和偏移
            self.sattgt[ii]['initialepoch'] = self.sattgt[ii]['sattle'].jdsatepoch - calculate_jday(1950, 1, 0, 0, 0, 0)
            self.sattgt[ii]['initialjulian'] = self.sattgt[ii]['sattle'].jdsatepoch
            self.sattgt[ii]['offset'] = (ConjStartJulian - self.sattgt[ii]['initialjulian']) * 1440  # 时间偏移（分钟）

            # 初始化 SGP4 轨道参数
            self.sattgt[ii]['struc'] = self.sattgt[ii]['sattle']

            # 计算目标卫星的位置和速度
            _, p, v = self.sattgt[ii]['struc'].sgp4(self.sattgt[ii]['initialjulian'], self.sattgt[ii]['offset'] / 1440.0)
            self.tgtpnow[ii, :] = p
            self.tgtvnow[ii, :] = v

            # 计算目标卫星与对象卫星的相对位置和速度
            for kk in range(objNum):
                drtemp = self.tgtpnow[ii, :] - objpnow[kk, :]
                dvtemp = self.tgtvnow[ii, :] - objvnow[kk, :]
                rv = np.sum(drtemp * dvtemp)
                rr = np.linalg.norm(drtemp)
                self.RelativePx[ii, kk] = drtemp[0]
                self.RelativePy[ii, kk] = drtemp[1]
                self.RelativePz[ii, kk] = drtemp[2]
                self.RelativeVx[ii, kk] = dvtemp[0]
                self.RelativeVy[ii, kk] = dvtemp[1]
                self.RelativeVz[ii, kk] = dvtemp[2]
                self.CurrentRangeRate[ii, kk] = rv / rr
                self.CurrentRange[ii, kk] = rr

    

if __name__ == "__main__":
    from setup_TLEfiles import generate_objSat_from_temptle, generate_tarSat_from_temptle
    from tools.date_trans import date_to_julian
    import datetime
    ConjStartDate = datetime.datetime(2024, 11, 15, 0, 0, 0)
    ConjStartJulian = date_to_julian(ConjStartDate)
    tempfile = 'temptle.tle'
    ObjSat_list, objsma= generate_objSat_from_temptle(tempfile)
    TgtSat_list = generate_tarSat_from_temptle(tempfile, ObjSat_list, objsma)
    objNum, tgtNum = len(ObjSat_list), len(TgtSat_list)
    objSatDetail = ObjSatDetail(objNum)
    objSatDetail.calculate_objSat_detail(ObjSat_list, objNum, ConjStartJulian)
    objCurrentRange, objCurrentRangeRate = objSatDetail.compute_obj_relative_values(objNum)
    tgtSatDetail = TgtSatDetail(tgtNum, objNum)
    tgtSatDetail.calculate_tgtSat_detail(TgtSat_list, objNum, ConjStartJulian, objSatDetail.objpnow, objSatDetail.objvnow)
    print("=================Object Sat and Target Sat parameters' initialization completed!=========================")
    
