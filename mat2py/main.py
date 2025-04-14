import numpy as np
import datetime
from scipy.io import loadmat
from download_TLEs_data import download_tle
from setup_TLEfiles import *
from crash_analysis_prepare import *
from Propgation_analysis import *
from dataestr import *
from tools.date_trans import *
from tools.common_tools import *
#initialize
ConjStartDate = datetime.datetime(2025, 4, 7, 0, 0, 0)
ConjEndDate = datetime.datetime(2025, 4, 14, 0, 0, 0)
ConjStartJulian = date_to_julian(ConjStartDate)
ConjEndJulian = date_to_julian(ConjEndDate)

#获取卫星数据并生成卫星对象
download_tle()
tempfile = 'temptle.tle'
ObjSat_list, objsma = generate_objSat_from_temptle(tempfile)
TgtSat_list = generate_tarSat_from_temptle(tempfile, ObjSat_list, objsma)


# 初始化用于卫星对象间联合碰撞分析变量
objNum, tgtNum = len(ObjSat_list), len(TgtSat_list)
PropTimeStep = 5 

# data preparation
timeVec=initialize_time_vector(ConjStartDate, ConjEndDate, PropTimeStep)
objSatDetail = ObjSatDetail(objNum)
objSatDetail.calculate_objSat_detail(ObjSat_list, objNum, ConjStartJulian)
tgtSatDetail = TgtSatDetail(tgtNum, objNum)
tgtSatDetail.calculate_tgtSat_detail(TgtSat_list, objNum, ConjStartJulian, objSatDetail.objpnow, objSatDetail.objvnow)

analysis_threshold = 300000  # 分析阈值（km）
conj_range_threshold = 50000  # 会合距离阈值（km）
min_dis_threshold = 8000  # 最小距离阈值（km）
duration_days = 30  # 评估持续时间
time_step_minutes = 5  # 时间步长
report_file = "conjunction_report.csv"
print(f"开始会合评估，总计 {len(timeVec)} 个时间步")
conjunction_assessment(objSatDetail, tgtSatDetail, timeVec, ConjStartDate, time_step_minutes, 
                        analysis_threshold, conj_range_threshold, min_dis_threshold, report_file)
print("Conjunction assessment completed.")