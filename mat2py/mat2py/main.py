import numpy as np
import datetime
from scipy.io import loadmat
from download_TLEs_data import download_tle
from setup_TLEfiles import *
from crash_analysis import *
import ConjunctionReportOutput
from readManeuverDataNewTLEFormat import *
from dataestr import *
# 定义全局变量
d2r = np.pi / 180
mu = 398600.8
req = 6378.135
# 其他全局变量
jdate0 = None
isun = None
imoon = None
idrag = None
isrp = None
smu = None
mmu = None
ad76 = None
thrustvar = None
pressurepoly = None
thrustpoly = None
isppoly = None
maxfuel = None
satdrymass = None
pmat = None
perr = None
lgrav = None
mgrav = None
j2 = None
ccoef = None
scoef = None

# ========= Evaluation Target ==========

ConjStartDate = datetime.datetime(2024, 11, 15, 0, 0, 0)
ConjEndDate = datetime.datetime(2024, 11, 17, 0, 0, 0)
DateTrack = ConjStartDate.replace(hour=0, minute=0, second=0)


def date_to_julian(dt):
    year = dt.year
    month = dt.month
    day = dt.day
    hour = dt.hour
    minute = dt.minute
    second = dt.second
    if month <= 2:
        month += 12
        year -= 1
    A = year // 100
    B = 2 - A + A // 4
    jd_int = int(365.25 * (year + 4716)) + int(30.6001 * (month + 1)) + day + B - 1524.5
    jd_frac = (hour + (minute + second / 60) / 60) / 24
    return jd_int + jd_frac


def read_tle(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    return lines
def read_text_file(filename):
    with open(filename, 'r') as file:
        content = file.read()
    return content
# 计算儒略日
def jday(year, month, day, hour=0, minute=0, second=0):
    return (datetime.datetime(year, month, day, hour, minute, second) - datetime.datetime(year, month, 1)).total_seconds() / 86400 + 1

ConjStartJulian = date_to_julian(ConjStartDate)
ConjEndJulian = date_to_julian(ConjEndDate)

ObjCatID = np.array([['16493'], ['49256'], ['37484'], ['04737'], ['43910']])
objNum=len(ObjCatID)

# 文件名
tempfile="temptle.tle"
objfile ="objtle.tle"
fmnvr = read_text_file('satmaneuver.txt')
data = loadmat('PropulsionDataNewSet.mat')
fmaneverdata = 'eop_data.txt'
# reportfile = loadmat('cdm.dat')
smaThresHold = 50

ObjSat, IsObj=process_satellite_tle(tempfile, objfile, ObjCatID)
tgtNum,TgtSat=setup_TLEfiles()

# 初始化参数
PropTimeStep = 5  # in minutes
# ConjunctionReportOutput(0,reportfile)
thres = 0.001 / 60  # 100 milliseconds
maxcount = 50
dxscalar = np.array([1, 1, 1]) / 60  # one second slope
tmin = np.zeros((tgtNum, objNum))
tmax = np.ones((tgtNum, objNum)) * PropTimeStep
minDisThres = 5.0  # Threshold to report the data, in km
analysisThres = 3000  # Upper limit for conjunction analysis
ConjRangeThres = 1000  # For predicted range

timeVec=initialize_time_vector(ConjStartDate, ConjEndDate, PropTimeStep)
sat_objects, objpnow, objvnow, CurrentOr, CurrentOt, CurrentOh=initialize_satellites(ObjSat, ConjStartJulian)
RelativePx, RelativePy, RelativePz, RelativeVx, RelativeVy, RelativeVz, CurrentRangeRate, CurrentRange=compute_relative_positions(TgtSat, objpnow, objvnow, ConjStartJulian)
objCurrentRange, objCurrentRangeRate=compute_obj_relative(objpnow, objvnow)

maneuver_data = load_maneuver_data(fmaneverdata, ObjSat)#dict
TotalManeuverData=maneuver_data["TotalManeuverData"]
if TotalManeuverData>0:
    print(f"共有{TotalManeuverData}个推进数据可用")
else:
    print("无推进数据")
dtnumber=datestr(DateTrack,'dd-mmm-yyyy')
print(f"在{dtnumber}启动联合评估")
