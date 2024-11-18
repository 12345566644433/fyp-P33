import os
import numpy as np
from datetime import datetime, timedelta
from scipy.integrate import odeint

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


# 定义常量
flat = 1 / 298.257
dsun = 696000
dmoon = 1738
smu = 132712438000
mmu = 4902.793
aunit = 149597870
omega = 7.292115486e-5
ps = 0.0044
dtr = np.pi / 180
rtd = 180 / np.pi
atr = dtr / 3600

# 读取数据函数，需自行实现
def read76():
    # 请根据实际情况返回数据
    pass

_, ad76 = read76()

# 初始化 Runge-Kutta 方法
rkcoef = 1
# 微分方程的数量
neq = 6

# 读取重力模型系数的函数，需自行实现
def readegm(filename):
    # 请根据实际情况读取文件并返回系数
    pass

ccoef, scoef = readegm('egm96.dat')
j2 = -ccoef[2, 0]

# 传播设置和误差容限
tetol = 1e-8
lgrav = 2  # 重力模型的度数，范围 0 到 18
mgrav = 0  # 重力模型的阶数，范围 0 到 18
isun = 0  # 太阳扰动
imoon = 0  # 月球扰动
isrp = 0  # 日光压力扰动
idrag = 0  # 阻力扰动
thrustvar = np.zeros(3)
thrustvar[2] = 9.81

# ae, cd, PropStartjday, PropEndtjday, satmass 需要定义





# 假设这些函数和变量已经定义
# jday, sgp4, rvteme2eci, computeposcov 等函数需要自行定义

for tstep in range(1, len(timeVec)):
    tsince = timeVec[tstep]

    DateNow = datetime.fromordinal(ConjStartDate.toordinal() + int(timeVec[tstep - 1]))
    DateNext = datetime.fromordinal(ConjStartDate.toordinal() + int(timeVec[tstep]))
    jdaynow = jday(DateNow.year, DateNow.month, DateNow.day, DateNow.hour, DateNow.minute, DateNow.second)
    jdaynext = jday(DateNext.year, DateNext.month, DateNext.day, DateNext.hour, DateNext.minute, DateNext.second)

    if (DateNow - DateTrack).days >= 1.0:
        DateTrack = datetime(DateNow.year, DateNow.month, DateNow.day)
        print(f'The Conjunction Assessment Process Now at {DateTrack.strftime("%d-%b-%Y")}')

    # ----------------------------------------------------
    # 目标卫星传播
    for objk in range(objNum):
        # 进入机动积分的两个条件
        TimeLengthOfManeuver[objk] = 0  # 重置

        if ObjPropCond[objk] > 0:
            # 条件 1
            if (PropStartjday[objk] >= jdaynow) and (PropStartjday[objk] <= jdaynext):
                isSwitch[objk] = 1  # 指示传播方法切换到 J2 扰动模型

                # 计算机动前的额外偏移（以秒为单位）
                TimeToManeuver[objk] = (PropStartjday[objk] - jdaynow) * 1440  # 转换为分钟
                _, initialp, initialv = sgp4(satobj[objk].struc, satobj[objk].offset + tsince + TimeToManeuver[objk] - PropTimeStep)
                manvrpteme[:, objk] = initialp
                manvrvteme[:, objk] = initialv

                # teme2eci 转换
                DateConv = DateNow + timedelta(minutes=TimeToManeuver[objk])
                rtemp, vtemp = rvteme2eci(DateConv, fname, manvrpteme[:, objk], manvrvteme[:, objk])
                manvrptemp[:, objk] = rtemp
                manvrvtemp[:, objk] = vtemp

                poscovtemp = computeposcov(manvrpteme[:, objk], manvrvteme[:, objk],
                                            (satobj[objk].offset + tsince + TimeToManeuver[objk] - PropTimeStep) / 1440,
                                            pcovoffset[objk, :3], leotlecov)
                pcovoffset[objk, :3] = np.diag(poscovtemp)
                pcovoffset[objk, 3:6] = pcovoffset[objk, :3] * np.linalg.norm(manvrvtemp) / np.linalg.norm(manvrptemp) * np.linalg.norm(manvrvtemp) / np.linalg.norm(manvrptemp)
                pcovnow = np.block([[poscovtemp, np.zeros((3, 3))], [np.zeros((3, 3)), poscovtemp * np.linalg.norm(manvrvtemp) / np.linalg.norm(manvrptemp) * np.linalg.norm(manvrvtemp) / np.linalg.norm(manvrptemp)]])

                if PropEndjday[objk] > jdaynext:
                    # 需要下一个循环继续机动过程
                    TimeLengthOfManeuver[objk] = ManeuverDuration[objk] - (PropTimeStep - TimeToManeuver[objk]) * 60
                else:
                    # 此循环将结束
                    isSwitch[objk] = 2  # 指示机动结束，但需要处理 SGP4 TLE
                    TimeLengthOfManeuver[objk] = ManeuverDuration[objk]



# 假设这些函数和变量已经定义
# jday, sgp4, rvteme2eci, computeposcov 等函数需要自行定义

for tstep in range(1, len(timeVec)):
    tsince = timeVec[tstep]

    DateNow = datetime.fromordinal(ConjStartDate.toordinal() + int(timeVec[tstep - 1]))
    DateNext = datetime.fromordinal(ConjStartDate.toordinal() + int(timeVec[tstep]))
    
    jdaynow = jday(DateNow.year, DateNow.month, DateNow.day, DateNow.hour, DateNow.minute, DateNow.second)
    jdaynext = jday(DateNext.year, DateNext.month, DateNext.day, DateNext.hour, DateNext.minute, DateNext.second)

    if (DateNow - DateTrack).days >= 1.0:
        DateTrack = datetime(DateNow.year, DateNow.month, DateNow.day)
        print(f'The Conjunction Assessment Process Now at {DateTrack.strftime("%d-%b-%Y")}')
    
    # ----------------------------------------------------
    # 目标卫星传播
    for objk in range(objNum):
        TimeLengthOfManeuver[objk] = 0  # 重置

        if ObjPropCond[objk] > 0:
            # 条件 1
            if (PropStartjday[objk] >= jdaynow) and (PropStartjday[objk] <= jdaynext):
                isSwitch[objk] = 1  # 指示传播方法切换到 J2 扰动模型

                # 计算机动前的额外偏移（以秒为单位）
                TimeToManeuver[objk] = (PropStartjday[objk] - jdaynow) * 1440  # 转换为分钟
                _, initialp, initialv = sgp4(satobj[objk].struc, satobj[objk].offset + tsince + TimeToManeuver[objk] - PropTimeStep)
                manvrpteme[:, objk] = initialp
                manvrvteme[:, objk] = initialv

                # teme2eci 转换
                DateConv = DateNow + timedelta(minutes=TimeToManeuver[objk])
                rtemp, vtemp = rvteme2eci(DateConv, fname, manvrpteme[:, objk], manvrvteme[:, objk])
                manvrptemp[:, objk] = rtemp
                manvrvtemp[:, objk] = vtemp

                poscovtemp = computeposcov(manvrpteme[:, objk], manvrvteme[:, objk],
                                            (satobj[objk].offset + tsince + TimeToManeuver[objk] - PropTimeStep) / 1440,
                                            pcovoffset[objk, :3], leotlecov)
                                            
                pcovoffset[objk, :3] = np.diag(poscovtemp)
                pcovoffset[objk, 3:6] = pcovoffset[objk, :3] * np.linalg.norm(manvrvtemp[:, objk]) / np.linalg.norm(manvrptemp[:, objk]) * np.linalg.norm(manvrvtemp[:, objk]) / np.linalg.norm(manvrptemp[:, objk])
                pcovnow = np.block([[poscovtemp, np.zeros((3, 3))], [np.zeros((3, 3)), poscovtemp * np.linalg.norm(manvrvtemp[:, objk]) / np.linalg.norm(manvrptemp[:, objk]) * np.linalg.norm(manvrvtemp[:, objk]) / np.linalg.norm(manvrptemp[:, objk])]])

                if PropEndjday[objk] > jdaynext:
                    TimeLengthOfManeuver[objk] = ManeuverDuration[objk] - (PropTimeStep - TimeToManeuver[objk]) * 60
                else:
                    isSwitch[objk] = 2  # 指示机动结束，但需要处理 SGP4 TLE
                    TimeLengthOfManeuver[objk] = ManeuverDuration[objk]

                # 生成日历日期字符串
                cdstr = DateConv.strftime("%d-%b-%Y")
                # 生成通用时间字符串
                utstr = DateConv.strftime("%H:%M:%S.%f")[:-3]  # 只保留毫秒
                print(f'Date: {cdstr}\tTime: {utstr}\t Begin maneuver process for Obj ID: {ObjSat[objk].CatID}')
                
                maxfuel = objmaxfuel[objk]
                satdrymass = objdrymass[objk]
                jdate0 = PropStartjday[objk]
                gst0 = gast(jdate0)
                pressurepoly = pressureprofile[objk]
                thrustpoly = thrustprofile[objk]
                isppoly = ispprofile[objk]

                MnvrTimeVec = np.linspace(0, TimeLengthOfManeuver, int(TimeLengthOfManeuver) + 1)
                for tk in range(len(MnvrTimeVec) - 1):
                    bcoeff = 1.0e-6 * ae[objk] * cd[objk] / satmass[objk]

                    if PropMethod[objk] == 1:  # +ve, -ve 速度方向燃烧
                        thrustvar[4] = PropVsign[objk]
                        perr = 0
                        x0 = np.hstack((manvrptemp[:, objk], manvrvtemp[:, objk], satmass[objk]))
                        t_span = (0, MnvrTimeVec[tk + 1] - MnvrTimeVec[tk])
                        yfinal = odeint(ceqm1hohconj, x0, t_span, tfirst=True)
                    elif PropMethod[objk] == 2:  # 平面变换机动
                        vmb = np.cross(manvrptemp[:, objk], manvrvtemp[:, objk])
                        thrustvar[4:6] = vmb / np.linalg.norm(vmb)
                        x0 = np.hstack((manvrptemp[:, objk], manvrvtemp[:, objk], satmass[objk]))
                        t_span = (0, MnvrTimeVec[tk + 1] - MnvrTimeVec[tk])
                        yfinal = odeint(ceqm1vsigndirectconj, x0, t_span, tfirst=True)

                    # 协方差传播
                    dj2dx = j2dfdx(np.hstack((manvrptemp[:, objk], manvrvtemp[:, objk])))
                    dUdx = np.block([[np.zeros((3, 3))], [np.eye(3)], [dj2dx], [np.zeros((3, 3))]])
                    qmat = 1e-12 * np.eye(3) * np.linalg.norm(manvrvtemp[:, objk])**2
                    gmat = np.block([[np.zeros((3, 3))], [np.eye(3)]])
                    fuelused = 100 - (satmass[objk] - satdrymass) / maxfuel * 100
                    estThrust = cal_poly(cal_poly(fuelused, pressurepoly), thrustpoly)
                    aTmat = thrustErrRatio * estThrust**2 / (satmass[objk]**2) * np.eye(3)

                    p0 = pcovnow.flatten()
                    xx = odeint(lambda t, pv: covprop(t, pv, dUdx, gmat, qmat, dmat, aTmat), p0, t_span, tfirst=True)
                    pcovnow = xx[-1].reshape((6, 6))

                    satmass[objk] = yfinal[-1, -1]
                    manvrptemp[:, objk] = yfinal[-1, :3]
                    manvrvtemp[:, objk] = yfinal[-1, 3:6]

                # 存储协方差数据
                pcovstored[objk].p = pcovnow
                pcovoffset[objk, :] = np.diag(pcovnow)

                sigr = np.sqrt(pcovoffset[objk, :3])
                sigv = np.sqrt(pcovoffset[objk, 3:6])
                DateConv = DateNow + timedelta(minutes=TimeToManeuver[objk] + TimeLengthOfManeuver)
                rtrue, vtrue = rveci2teme(DateConv, fname, manvrptemp[:, objk], manvrvtemp[:, objk])
                satobj[objk].offset = -tsince  # 重置偏移

                # 更新对象
                manvrjday = jday(DateConv.year, DateConv.month, DateConv.day, DateConv.hour, DateConv.minute, DateConv.second)
                satobj[objk].struc.julian = manvrjday
                satobj[objk].sattle.julian = manvrjday
                satobj[objk].struc.epoch = manvrjday - jday(1950, 1, 0, 0, 0, 0)
                satobj[objk].initialepoch = satobj[objk].sattle.julian - jday(1950, 1, 0, 0, 0, 0)
                satobj[objk].initialjulian = satobj[objk].sattle.julian

                if isSwitch[objk] != 2:
                    print('Performing intermediate PV2TLE estimation process.....')
                    xe, Pout, _, _ = kfalgoca(100, rtrue.T, vtrue.T, satobj[objk].struc, sigr, sigv, 1)
                    satobj[objk].struc = sgp4init(72, satobj[objk].struc, satobj[objk].sattle.bstar, xe[2], satobj[objk].initialepoch, xe[4], xe[3], xe[6], xe[1], xe[5])
                    _, p, v = sgp4(satobj[objk].struc, 0.0)
                    print('Completed intermediate PV2TLE estimation process')

            elif (isSwitch[objk] == 1) and (PropEndjday[objk] > jdaynext):
                # 长时间机动的处理
                maxfuel = objmaxfuel[objk]
                satdrymass = objdrymass[objk]
                jdate0 = PropStartjday[objk]
                gst0 = gast(jdate0)
                pressurepoly = pressureprofile[objk]
                thrustpoly = thrustprofile[objk]
                isppoly = ispprofile[objk]

                pcovnow = pcovstored[objk].p
                TimeLengthOfManeuver[objk] = PropTimeStep * 60
                MnvrTimeVec = np.linspace(0, TimeLengthOfManeuver, int(np.ceil(TimeLengthOfManeuver)) + 1)
                for tk in range(len(MnvrTimeVec) - 1):
                    bcoeff = 1.0e-6 * ae[objk] * cd[objk] / satmass[objk]

                    if PropMethod[objk] == 1:  # +ve, -ve 速度方向燃烧
                        thrustvar[4] = PropVsign[objk]
                        perr = 0
                        x0 = np.hstack((manvrptemp[:, objk], manvrvtemp[:, objk], satmass[objk]))
                        t_span = (0, MnvrTimeVec[tk + 1] - MnvrTimeVec[tk])
                        yfinal = odeint(ceqm1hohconj, x0, t_span, tfirst=True)
                    elif PropMethod[objk] == 2:  # 平面变换机动
                        vmb = np.cross(manvrptemp[:, objk], manvrvtemp[:, objk])
                        thrustvar[4:6] = vmb / np.linalg.norm(vmb)
                        x0 = np.hstack((manvrptemp[:, objk], manvrvtemp[:, objk], satmass[objk]))
                        t_span = (0, MnvrTimeVec[tk + 1] - MnvrTimeVec[tk])
                        yfinal = odeint(ceqm1vsigndirectconj, x0, t_span, tfirst=True)

                    # 协方差传播
                    dj2dx = j2dfdx(np.hstack((manvrptemp[:, objk], manvrvtemp[:, objk])))
                    dUdx = np.block([[np.zeros((3, 3))], [np.eye(3)], [dj2dx], [np.zeros((3, 3))]])
                    qmat = 1e-12 * np.eye(3) * np.linalg.norm(manvrvtemp[:, objk])**2
                    gmat = np.block([[np.zeros((3, 3))], [np.eye(3)]])
                    fuelused = 100 - (satmass[objk] - satdrymass) / maxfuel * 100
                    estThrust = cal_poly(cal_poly(fuelused, pressurepoly), thrustpoly)
                    aTmat = thrustErrRatio * estThrust**2 / (satmass[objk]**2) * np.eye(3)

                    p0 = pcovnow.flatten()
                    xx = odeint(lambda t, pv: covprop(t, pv, dUdx, gmat, qmat, dmat, aTmat), p0, t_span, tfirst=True)
                    pcovnow = xx[-1].reshape((6, 6))

                    satmass[objk] = yfinal[-1, -1]
                    manvrptemp[:, objk] = yfinal[-1, :3]
                    manvrvtemp[:, objk] = yfinal[-1, 3:6]

                # 存储协方差数据
                pcovstored[objk].p = pcovnow
                pcovoffset[objk, :] = np.diag(pcovnow)

                DateConv = DateNow + timedelta(minutes=TimeLengthOfManeuver)
                rtrue, vtrue = rveci2teme(DateConv, fname, manvrptemp[:, objk], manvrvtemp[:, objk])
                p = rtrue
                v = vtrue
                satobj[objk].offset = -tsince  # 重置偏移

            elif (PropEndjday[objk] >= jdaynow) and (PropEndjday[objk] <= jdaynext):  # 条件 2
                isSwitch[objk] = 2  # 指示机动结束，但需要处理 SGP4 TLE
                TimeToManeuver[objk] = 0
                pcovnow = pcovstored[objk].p  # 还原协方差数据
                TimeLengthOfManeuver[objk] = (PropEndjday[objk] - jdaynow) * 86400  # 转换为秒
                DateConv = DateNow + timedelta(seconds=TimeLengthOfManeuver[objk])
                
                # 生成日历日期字符串
                cdstr = DateConv.strftime("%d-%b-%Y")
                # 生成通用时间字符串
                utstr = DateConv.strftime("%H:%M:%S.%f")[:-3]
                print(f'Date: {cdstr}\tTime: {utstr}\t End maneuver process for Obj ID: {ObjSat[objk].CatID}')

                maxfuel = objmaxfuel[objk]
                satdrymass = objdrymass[objk]
                jdate0 = jdaynow
                gst0 = gast(jdate0)                
                pressurepoly = pressureprofile[objk]
                thrustpoly = thrustprofile[objk]
                isppoly = ispprofile[objk]

                MnvrTimeVec = np.linspace(0, TimeLengthOfManeuver[objk], int(np.ceil(TimeLengthOfManeuver[objk])) + 1)
                for tk in range(len(MnvrTimeVec) - 1):
                    bcoeff = 1.0e-6 * ae[objk] * cd[objk] / satmass[objk]
                    
                    if PropMethod[objk] == 1:  # +ve, -ve 速度方向燃烧
                        thrustvar[4] = PropVsign[objk]
                        perr = 0
                        x0 = np.hstack((manvrptemp[:, objk], manvrvtemp[:, objk], satmass[objk]))
                        t_span = (0, MnvrTimeVec[tk + 1] - MnvrTimeVec[tk])
                        yfinal = odeint(ceqm1hohconj, x0, t_span, tfirst=True)
                    elif PropMethod[objk] == 2:  # 平面变换机动
                        vmb = np.cross(manvrptemp[:, objk], manvrvtemp[:, objk])
                        thrustvar[4:6] = vmb / np.linalg.norm(vmb)
                        x0 = np.hstack((manvrptemp[:, objk], manvrvtemp[:, objk], satmass[objk]))
                        t_span = (0, MnvrTimeVec[tk + 1] - MnvrTimeVec[tk])
                        yfinal = odeint(ceqm1vsigndirectconj, x0, t_span, tfirst=True)

# 继续进行其他处理或者数据存储




# 假设这些函数和数据结构已经定义
# 需要实现或导入的函数: jday, sgp4, rvteme2eci, computeposcov, kfalgoca, sgp4init, cal_poly, covprop

for tstep in range(1, len(timeVec)):
    tsince = timeVec[tstep]
    
    DateNow = datetime.fromordinal(ConjStartDate.toordinal() + int(timeVec[tstep - 1]))
    DateNext = datetime.fromordinal(ConjStartDate.toordinal() + int(timeVec[tstep]))
    
    jdaynow = jday(DateNow.year, DateNow.month, DateNow.day, DateNow.hour, DateNow.minute, DateNow.second)
    jdaynext = jday(DateNext.year, DateNext.month, DateNext.day, DateNext.hour, DateNext.minute, DateNext.second)
    
    if (DateNow - DateTrack).days >= 1.0:
        DateTrack = datetime(DateNow.year, DateNow.month, DateNow.day)
        print(f'The Conjunction Assessment Process Now at {DateTrack.strftime("%d-%b-%Y")}')
    
    # ----------------------------------------------------
    # 目标卫星传播
    for objk in range(objNum):
        TimeLengthOfManeuver[objk] = 0  # 重置

        if ObjPropCond[objk] > 0:
            # 条件 1
            if (PropStartjday[objk] >= jdaynow) and (PropStartjday[objk] <= jdaynext):
                isSwitch[objk] = 1  # 指示传播方法切换到 J2 扰动模型

                # 计算机动前的额外偏移（以秒为单位）
                TimeToManeuver[objk] = (PropStartjday[objk] - jdaynow) * 1440  # 转换为分钟
                _, initialp, initialv = sgp4(satobj[objk].struc, satobj[objk].offset + tsince + TimeToManeuver[objk] - PropTimeStep)
                manvrpteme[:, objk] = initialp
                manvrvteme[:, objk] = initialv

                # teme2eci 转换
                DateConv = DateNow + timedelta(minutes=TimeToManeuver[objk])
                rtemp, vtemp = rvteme2eci(DateConv, fname, manvrpteme[:, objk], manvrvteme[:, objk])
                manvrptemp[:, objk] = rtemp
                manvrvtemp[:, objk] = vtemp

                # 协方差计算
                poscovtemp = computeposcov(manvrpteme[:, objk], manvrvteme[:, objk],
                                            (satobj[objk].offset + tsince + TimeToManeuver[objk] - PropTimeStep) / 1440,
                                            pcovoffset[objk, :3], leotlecov)
                pcovoffset[objk, :3] = np.diag(poscovtemp)
                pcovoffset[objk, 3:6] = pcovoffset[objk, :3] * np.linalg.norm(manvrvtemp[:, objk]) / np.linalg.norm(manvrptemp[:, objk]) * np.linalg.norm(manvrvtemp[:, objk]) / np.linalg.norm(manvrptemp[:, objk])
                pcovnow = np.block([[poscovtemp, np.zeros((3, 3))], [np.zeros((3, 3)), poscovtemp * np.linalg.norm(manvrvtemp[:, objk]) / np.linalg.norm(manvrptemp[:, objk]) * np.linalg.norm(manvrvtemp[:, objk]) / np.linalg.norm(manvrptemp[:, objk])]])

                if PropEndjday[objk] > jdaynext:
                    TimeLengthOfManeuver[objk] = ManeuverDuration[objk] - (PropTimeStep - TimeToManeuver[objk]) * 60
                else:
                    isSwitch[objk] = 2  # 指示机动结束，但需要处理 SGP4 TLE
                    TimeLengthOfManeuver[objk] = ManeuverDuration[objk]

                # 生成日历和时间字符串
                cdstr = DateConv.strftime("%d-%b-%Y")
                utstr = DateConv.strftime("%H:%M:%S.%f")[:-3]
                print(f'Date: {cdstr}\tTime: {utstr}\t Begin maneuver process for Obj ID: {ObjSat[objk].CatID}')
                
                maxfuel = objmaxfuel[objk]
                satdrymass = objdrymass[objk]
                jdate0 = PropStartjday[objk]
                gst0 = gast(jdate0)
                pressurepoly = pressureprofile[objk]
                thrustpoly = thrustprofile[objk]
                isppoly = ispprofile[objk]

                MnvrTimeVec = np.linspace(0, TimeLengthOfManeuver, int(TimeLengthOfManeuver) + 1)
                for tk in range(len(MnvrTimeVec) - 1):
                    bcoeff = 1.0e-6 * ae[objk] * cd[objk] / satmass[objk]

                    if PropMethod[objk] == 1:  # +ve, -ve 速度方向燃烧
                        thrustvar[4] = PropVsign[objk]
                        perr = 0
                        x0 = np.hstack((manvrptemp[:, objk], manvrvtemp[:, objk], satmass[objk]))
                        t_span = (0, MnvrTimeVec[tk + 1] - MnvrTimeVec[tk])
                        yfinal = odeint(ceqm1hohconj, x0, t_span)
                    elif PropMethod[objk] == 2:  # 平面变换机动
                        vmb = np.cross(manvrptemp[:, objk], manvrvtemp[:, objk])
                        thrustvar[4:6] = vmb / np.linalg.norm(vmb)
                        x0 = np.hstack((manvrptemp[:, objk], manvrvtemp[:, objk], satmass[objk]))
                        t_span = (0, MnvrTimeVec[tk + 1] - MnvrTimeVec[tk])
                        yfinal = odeint(ceqm1vsigndirectconj, x0, t_span)

                    # 协方差传播
                    dj2dx = j2dfdx(np.hstack((manvrptemp[:, objk], manvrvtemp[:, objk])))
                    dUdx = np.block([[np.zeros((3, 3))], [np.eye(3)], [dj2dx], [np.zeros((3, 3))]])
                    qmat = 1e-12 * np.eye(3) * np.linalg.norm(manvrvtemp[:, objk])**2
                    gmat = np.block([[np.zeros((3, 3))], [np.eye(3)]])
                    fuelused = 100 - (satmass[objk] - satdrymass) / maxfuel * 100
                    estThrust = cal_poly(cal_poly(fuelused, pressurepoly), thrustpoly)
                    aTmat = thrustErrRatio * estThrust**2 / (satmass[objk]**2) * np.eye(3)

                    p0 = pcovnow.flatten()
                    xx = odeint(lambda t, pv: covprop(t, pv, dUdx, gmat, qmat, dmat, aTmat), p0, t_span)
                    pcovnow = xx[-1].reshape((6, 6))

                    satmass[objk] = yfinal[-1, -1]
                    manvrptemp[:, objk] = yfinal[-1, :3]
                    manvrvtemp[:, objk] = yfinal[-1, 3:6]

                # 存储协方差数据
                pcovstored[objk].p = pcovnow
                pcovoffset[objk, :] = np.diag(pcovnow)

                sigr = np.sqrt(pcovoffset[objk, :3])
                sigv = np.sqrt(pcovoffset[objk, 3:6])
                DateConv = DateNow + timedelta(minutes=TimeToManeuver[objk] + TimeLengthOfManeuver)
                rtrue, vtrue = rveci2teme(DateConv, fname, manvrptemp[:, objk], manvrvtemp[:, objk])
                satobj[objk].offset = -tsince  # 重置偏移

                # 更新对象
                manvrjday = jday(DateConv.year, DateConv.month, DateConv.day, DateConv.hour, DateConv.minute, DateConv.second)
                satobj[objk].struc.julian = manvrjday
                satobj[objk].sattle.julian = manvrjday
                satobj[objk].struc.epoch = manvrjday - jday(1950, 1, 0, 0, 0, 0)
                satobj[objk].initialepoch = satobj[objk].sattle.julian - jday(1950, 1, 0, 0, 0, 0)
                satobj[objk].initialjulian = satobj[objk].sattle.julian

                if isSwitch[objk] != 2:
                    print('Performing intermediate PV2TLE estimation process.....')
                    xe, Pout, _, _ = kfalgoca(100, rtrue.T, vtrue.T, satobj[objk].struc, sigr, sigv, 1)
                    satobj[objk].struc = sgp4init(72, satobj[objk].struc, satobj[objk].sattle.bstar, xe[2], satobj[objk].initialepoch, xe[4], xe[3], xe[6], xe[1], xe[5])
                    _, p, v = sgp4(satobj[objk].struc, 0.0)
                    print('Completed intermediate PV2TLE estimation process')

            elif (isSwitch[objk] == 1) and (PropEndjday[objk] > jdaynext):
                # 长时间机动的处理
                maxfuel = objmaxfuel[objk]
                satdrymass = objdrymass[objk]
                jdate0 = PropStartjday[objk]
                gst0 = gast(jdate0)
                pressurepoly = pressureprofile[objk]
                thrustpoly = thrustprofile[objk]
                isppoly = ispprofile[objk]

                pcovnow = pcovstored[objk].p
                TimeLengthOfManeuver[objk] = PropTimeStep * 60
                MnvrTimeVec = np.linspace(0, TimeLengthOfManeuver, int(np.ceil(TimeLengthOfManeuver)) + 1)
                for tk in range(len(MnvrTimeVec) - 1):
                    bcoeff = 1.0e-6 * ae[objk] * cd[objk] / satmass[objk]

                    if PropMethod[objk] == 1:  # +ve, -ve 速度方向燃烧
                        thrustvar[4] = PropVsign[objk]
                        perr = 0
                        x0 = np.hstack((manvrptemp[:, objk], manvrvtemp[:, objk], satmass[objk]))
                        t_span = (0, MnvrTimeVec[tk + 1] - MnvrTimeVec[tk])
                        yfinal = odeint(ceqm1hohconj, x0, t_span)
                    elif PropMethod[objk] == 2:  # 平面变换机动
                        vmb = np.cross(manvrptemp[:, objk], manvrvtemp[:, objk])
                        thrustvar[4:6] = vmb / np.linalg.norm(vmb)
                        x0 = np.hstack((manvrptemp[:, objk], manvrvtemp[:, objk], satmass[objk]))
                        t_span = (0, MnvrTimeVec[tk + 1] - MnvrTimeVec[tk])
                        yfinal = odeint(ceqm1vsigndirectconj, x0, t_span)

                    # 协方差传播
                    dj2dx = j2dfdx(np.hstack((manvrptemp[:, objk], manvrvtemp[:, objk])))
                    dUdx = np.block([[np.zeros((3, 3))], [np.eye(3)], [dj2dx], [np.zeros((3, 3))]])
                    qmat = 1e-12 * np.eye(3) * np.linalg.norm(manvrvtemp[:, objk])**2
                    gmat = np.block([[np.zeros((3, 3))], [np.eye(3)]])
                    fuelused = 100 - (satmass[objk] - satdrymass) / maxfuel * 100
                    estThrust = cal_poly(cal_poly(fuelused, pressurepoly), thrustpoly)
                    aTmat = thrustErrRatio * estThrust**2 / (satmass[objk]**2) * np.eye(3)

                    p0 = pcovnow.flatten()
                    xx = odeint(lambda t, pv: covprop(t, pv, dUdx, gmat, qmat, dmat, aTmat), p0, t_span)
                    pcovnow = xx[-1].reshape((6, 6))

                    satmass[objk] = yfinal[-1, -1]
                    manvrptemp[:, objk] = yfinal[-1, :3]
                    manvrvtemp[:, objk] = yfinal[-1, 3:6]

                # 存储协方差数据
                pcovstored[objk].p = pcovnow
                pcovoffset[objk, :] = np.diag(pcovnow)

                # 进行其他必要的处理

# 在循环结束后，将进行后续处理和数据更新
