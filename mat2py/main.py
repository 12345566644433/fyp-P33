import numpy as np
import datetime
from scipy.io import loadmat
from download_TLEs_data import download_tle
from setup_TLEfiles import *
from crash_analysis_prepare import *
from readManeuverDataNewTLEFormat import *
from dataestr import *
from tools.date_trans import *
from tools.common_tools import *
#initialize
ConjStartDate = datetime.datetime(2024, 11, 15, 0, 0, 0)
ConjEndDate = datetime.datetime(2024, 11, 17, 0, 0, 0)
ConjStartJulian = date_to_julian(ConjStartDate)
ConjEndJulian = date_to_julian(ConjEndDate)

#获取卫星数据并生成卫星对象
download_tle()
tempfile = 'temptle.tle'
ObjSat_list, objsma = generate_objSat_from_temptle(tempfile)
TgtSat_list = generate_tarSat_from_temptle(tempfile, ObjSat_list, objsma)


# 初始化用于卫星对象间联合碰撞分析变量
objNum, tgtNum = len(ObjSat_list), len(TgtSat_list)
objNumVec = list(range(objNum))
tgtNumVec = list(range(tgtNum))
ObjPropCond = [1] * objNum
TgtPropCond = [1] * tgtNum
NextRange = [[0.0 for _ in range(objNum)] for _ in range(tgtNum)]
NextRangeRate = [[0.0 for _ in range(objNum)] for _ in range(tgtNum)]
CurrentRange = [[0.0 for _ in range(objNum)] for _ in range(tgtNum)]
CurrentRangeRate = [[0.0 for _ in range(objNum)] for _ in range(tgtNum)]
ConjFlag = [[1 for _ in range(objNum)] for _ in range(tgtNum)]
objConjFlag = [1] * objNum
ConjProb = [[0.0 for _ in range(objNum)] for _ in range(tgtNum)]
ConjDistance = [[0.0 for _ in range(objNum)] for _ in range(tgtNum)]

tgtpnow = [[0.0, 0.0, 0.0] for _ in range(tgtNum)]
tgtvnow = [[0.0, 0.0, 0.0] for _ in range(tgtNum)]
objpnow = [[0.0, 0.0, 0.0] for _ in range(objNum)]
objvnow = [[0.0, 0.0, 0.0] for _ in range(objNum)]
tgtpnext = [[0.0, 0.0, 0.0] for _ in range(tgtNum)]
tgtvnext = [[0.0, 0.0, 0.0] for _ in range(tgtNum)]
objpnext = [[0.0, 0.0, 0.0] for _ in range(objNum)]
objvnext = [[0.0, 0.0, 0.0] for _ in range(objNum)]

PropTimeStep = 5 
thres = 0.001 / 60 
maxcount = 50
dxscalar = np.array([1, 1, 1]) / 60  
tmin = np.zeros((tgtNum, objNum))
tmax = np.ones((tgtNum, objNum)) * PropTimeStep
minDisThres = 5.0  
analysisThres = 3000  
ConjRangeThres = 1000

# data preparation
timeVec=initialize_time_vector(ConjStartDate, ConjEndDate, PropTimeStep)
objSatDetail = ObjSatDetail(objNum)
objSatDetail.calculate_objSat_detail(ObjSat_list, objNum, ConjStartJulian)
objCurrentRange, objCurrentRangeRate = objSatDetail.compute_obj_relative_values(objNum)
tgtSatDetail = TgtSatDetail(tgtNum, objNum)
tgtSatDetail.calculate_tgtSat_detail(TgtSat_list, objNum, ConjStartJulian, objSatDetail.objpnow, objSatDetail.objvnow)
save_to_csv(objCurrentRange, objCurrentRangeRate, filename="obj_relative_values.csv")
