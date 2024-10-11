import numpy as np
import datetime
from scipy.io import loadmat
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
ObjCatID = np.array([['16493'], ['49256'], ['37484'], ['04737'], ['43910']])

# 转换日期格式
ConjStartDate = datetime.datetime(2024, 10, 7)
ConjEndDate = datetime.datetime(2024, 10, 10, 12)
DateTrack = ConjStartDate.replace(hour=0, minute=0, second=0)
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

ConjStartJulian = jday(ConjStartDate.year, ConjStartDate.month, ConjStartDate.day, ConjStartDate.hour, ConjStartDate.minute, ConjStartDate.second)
ConjEndJulian = jday(ConjEndDate.year, ConjEndDate.month, ConjEndDate.day, ConjEndDate.hour, ConjEndDate.minute, ConjEndDate.second)

objNum = ObjCatID.shape[0]

# 文件名
objfile =read_tle("objtle.tle")
tgtfile = read_tle('tgttle.tle')
tempfile = read_tle('temptle.tle')
fmnvr = read_text_file('satmaneuver.txt')
fmaneverdata = read_text_file('satmaneuver20220311.txt')
data = loadmat('PropulsionDataNewSet.mat')
reportfile = loadmat('cdm.dat')
