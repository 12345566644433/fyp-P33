import os
import numpy as np

# 全局变量
global jdate0, isun, imoon, idrag, isrp, smu, mmu
global bcoeff, thrustvar, pressurepoly, thrustpoly, isppoly, maxfuel, satdrymass
global dtr, rtd, atr
global flat, aunit, omega
global j2, lgrav, mgrav, gst0
global rkcoef, ccoef, scoef, ad76, csrp0

# 添加路径
os.sys.path.append('D:/fyp/Conjunction-Assessment/1toNum')

# 加载TLE误差数据
leotlestdpoly = np.load('leotlestdpoly.npy')  # 需确认数据格式
ucarcov = np.load('ucarcov.npy')  # 需确认数据格式

# 对象半径
RcSat = 7.5e-3
RcCube = 1.5e-3
RcRB = 25e-3
RcRBDEB = 15e-3
RcDEB = 10e-3
objNum = 10  # 需定义对象数量
RcObj = 55e-3 * np.ones(objNum)

# 读取CubeSat名称
cubesat_count = 0
cubesat = []

with open('cubesatname_data.txt', 'r') as fr:
    for fline in fr:
        fline = fline.strip()
        if fline:  # 如果行不为空
            cubesat_count += 1
            cubesat.append(fline)

# CubeSat 约定：最多3U用1m；最多12U用2m

# 加载推进数据
# propulsion_data = np.load('PropulsionData.mat', allow_pickle=True)  # 必须确认数据格式
isSwitch = np.zeros(objNum)  # 检查是否需要切换传播器
startidx = -5
endidx = 5
