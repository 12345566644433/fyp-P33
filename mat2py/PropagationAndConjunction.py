    import numpy as np
    
    # 假设需要的库和全局变量在这里初始化
    global jdate0, isun, imoon, idrag, isrp, smu, mmu
    global bcoeff, thrustvar, pressurepoly, thrustpoly, isppoly, maxfuel, satdrymass
    global dtr, rtd, atr
    global flat, aunit, omega
    global j2, lgrav, mgrav, gst0
    global rkcoef, ccoef, scoef, ad76, csrp0
    
    # 假设已有相应的路径和必须的导入
    # addpath('D:\\fyp\\Conjunction-Assessment\\1toNum')
    
    # 加载数据
    leotlestdpoly = np.load('leotlestdpoly.npy')
    ucarcov = np.load('ucarcov.npy')
    
    # 对象半径
    RcSat = 7.5 * 1e-3
    RcCube = 1.5e-3
    RcRB = 25 * 1e-3
    RcRBDEB = 15 * 1e-3
    RcDEB = 10 * 1e-3
    RcObj = 5 * 1e-3
    RcObj = 55 * 1e-3 * np.ones(objNum)
    
    # 读取Cubesat名称
    with open('cubesatname_data.txt', 'r') as fr:
        cubesat_count = 0
        for fline in fr:
            fline = fline.strip()
            if fline:
                cubesat_count += 1
                cubesat.append(fline)
    
    # CubeSat 1U到3U使用1m，12U使用2m
    
    # 加载推进数据
    isSwitch = np.zeros(objNum)  # 检查是否需要切换传播器
    startidx = -5
    endidx = 5
    
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
    # 引用的函数应在此处定义或导入
    # [_, ad76] = read76()
    
    # 初始化Runge-Kutta方法
    rkcoef = 1
    neq = 6  # 微分方程数量
    ccoef, scoef = readegm('egm96.dat')
    j2 = -ccoef[3, 0]
    
    # 传播设置和误差容限
    tetol = 1e-8
    lgrav = 2  # 重力模型的度数（只包含纬度项），范围0到18
    mgrav = 0  # 重力模型的阶数（经差项），范围0到18
    isun = 0   # 太阳扰动
    imoon = 0  # 月球扰动
    isrp = 0   # SRP扰动
    idrag = 0  # 拖曳扰动
    thrustvar[2] = 9.81
    
    for tstep in range(2, len(timeVec)):
        tsince = timeVec[tstep]
        
        DateNow = datetime.datetime.strptime(ConjStartDate, "%Y-%m-%d") + datetime.timedelta(minutes=timeVec[tstep-1])
        DateNext = datetime.datetime.strptime(ConjStartDate, "%Y-%m-%d") + datetime.timedelta(minutes=timeVec[tstep])
        
        jdaynow = jday(DateNow.year, DateNow.month, DateNow.day, DateNow.hour, DateNow.minute, DateNow.second)
        jdaynext = jday(DateNext.year, DateNext.month, DateNext.day, DateNext.hour, DateNext.minute, DateNext.second)
        
        if (DateNow - DateTrack).days >= 1:
            DateTrack = datetime.datetime(DateNow.year, DateNow.month, DateNow.day)
            print(f'The Conjunction Assessment Process Now at {DateTrack.strftime("%d-%b-%Y")}')
      
        # ----------------------------------------------------
        # 对象卫星传播
        for objk in range(1, objNum + 1):
            TimeLengthOfManeuver[objk] = 0  # 重置
    
            if ObjPropCond[objk] > 0:
                # 条件1
                if (PropStartjday[objk] >= jdaynow) and (PropStartjday[objk] <= jdaynext):
                    isSwitch[objk] = 1  # 表示在机动期间传播方法将切换到J2扰动模型
                    # 计算机动前的额外偏移时间（秒）
                    TimeToManeuver[objk] = (PropStartjday[objk] - jdaynow) * 1440  # 转换为分钟
    
                    # 在这里调用相应的函数或计算
                    # [_, initialp, initialv] = sgp4(satobj[objk].struc, satobj[objk].offset + tsince + TimeToManeuver[objk] - PropTimeStep)
                    manvrpteme[:, objk] = initialp.T
                    manvrvteme[:, objk] = initialv.T
    
                    # Teme转Eci转换
                    DateConv = datetime.datetime.strptime(DateNow.strftime("%Y-%m-%d %H:%M:%S"), "%Y-%m-%d %H:%M:%S") + datetime.timedelta(minutes=TimeToManeuver[objk])
                    rtemp, vtemp = rvteme2eci(DateConv, fname, manvrpteme[:, objk], manvrvteme[:, objk])
                    manvrptemp[:, objk] = rtemp
                    manvrvtemp[:, objk] = vtemp
    
                    # [_, objconjp, objconjv] = sgp4(satobj[objk].struc, (satobj[objk].offset + tsince + TimeToManeuver[objk] - PropTimeStep))
                    poscovtemp = computeposcov(manvrpteme[:, objk].T, manvrvteme[:, objk].T, (satobj[objk].offset + tsince + TimeToManeuver[objk] - PropTimeStep) / 1440, pcovoffset[objk, 1:3], leotlecov)
                    pcovoffset[objk, 1:3] = np.diag(poscovtemp)
                    pcovoffset[objk, 4:6] = pcovoffset[objk, 1:3] * np.linalg.norm(manvrvtemp) / np.linalg.norm(manvrptemp) * np.linalg.norm(manvrvtemp) / np.linalg.norm(manvrptemp)  # 使用位置和速度比率方法进行估计
                    pcovnow = np.block([[poscovtemp, np.zeros((3, 3))], [np.zeros((3, 3)), poscovtemp * np.linalg.norm(manvrvtemp) / np.linalg.norm(manvrptemp) * np.linalg.norm(manvrvtemp) / np.linalg.norm(manvrptemp)]])
    
                    if PropEndjday[objk] > jdaynext:
                        # 需要下一个循环以继续机动过程
                        TimeLengthOfManeuver[objk] = ManeuverDuration[objk] - (PropTimeStep - TimeToManeuver[objk]) * 60
                    else:
                        # 本循环将结束
                        isSwitch[objk] = 2  # 指示机动结束，但仍需要处理SGP4 TLE
                        TimeLengthOfManeuver[objk] = ManeuverDuration[objk]
    
                    # 创建日历日期字符串
                    cdstr = DateConv.strftime("%Y-%m-%d")
                    # 创建统一时间字符串
                    utstr = DateConv.strftime("%H:%M:%S.%f")
                    print(f'Date: {cdstr}\tTime: {utstr}\t Begin maneuver process for Obj ID: {ObjSat[objk].CatID}')
    
                    # 机动....
                    maxfuel = objmaxfuel[objk]
                    satdrymass = objdrymass[objk]
                    # 计算0小时UTC的儒略日
                    jdate0 = PropStartjday[objk]
                    # 计算UTC的格林威治恒星时间
                    gst0 = gast(jdate0)
                    pressurepoly = pressureprofile[objk]
                    thrustpoly = thrustprofile[objk]
                    isppoly = ispprofile[objk]
    
                    MnvrTimeVec = np.linspace(0, TimeLengthOfManeuver, TimeLengthOfManeuver + 1)
                    for tk in range(len(MnvrTimeVec) - 1):
                        # 对卫星2进行平面变化机动
                        bcoeff = 1.0e-6 * ae[objk] * cd[objk] / satmass[objk]
    
                        if PropMethod[objk] == 1:  # 这是用于正负速度方向的燃烧
                            thrustvar[3] = PropVsign[objk]
                            perr = 0
                            x0 = np.hstack((manvrptemp[:, objk], manvrvtemp[:, objk], satmass[objk]))
                            # [_, yfinal] = ode45(@ceqm1hohconj, [0 MnvrTimeVec[tk+1]-MnvrTimeVec[tk]], x0, options)
                        elif PropMethod[objk] == 2:  # 这是用于平面变化机动
                            vmb = np.cross(manvrptemp[:, objk], manvrvtemp[:, objk])
                            thrustvar[3:6] = vmb / np.linalg.norm(vmb)
                            pmat = np.eye(3)  # 无误差
                            x0 = np.hstack((manvrptemp[:, objk], manvrvtemp[:, objk], satmass[objk]))
                            # [_, yfinal] = ode45(@ceqm1vsigndirectconj, [0 MnvrTimeVec[tk+1]-MnvrTimeVec[tk]], x0, options)
    
                        # 协方差传播
                        dj2dx = j2dfdx(np.hstack((manvrptemp[:, objk], manvrvtemp[:, objk])))
                        dUdx = np.block([[np.zeros((3, 3)), np.eye(3)], [dj2dx, np.zeros((3, 3))]])
                        qmat = 1e-12 * np.eye(3) * np.linalg.norm(manvrvtemp[:, objk]) ** 2
                        gmat = np.block([[np.zeros((3, 3))], [np.eye(3)]])
                        # 模拟推进器误差方差
                        fuelused = 100 - (satmass[objk] - satdrymass) / maxfuel * 100
                        estThrust = cal_poly(cal_poly(fuelused, pressurepoly), thrustpoly)
                        dmat = np.block([[np.zeros((3, 3))], [np.eye(3)]])
                        aTmat = thrustErrRatio * estThrust**2 / (satmass[objk]**2) * np.eye(3)
                        
                        p0 = pcovnow.flatten()
                        # [_, xx] = ode45(lambda t, pv: covprop(t, pv, dUdx, gmat, qmat, dmat, aTmat), [0, MnvrTimeVec[tk+1] - MnvrTimeVec[tk]], p0, options)
                        pcovnow = xx[-1].reshape((6, 6))
    
                        satmass[objk] = yfinal[-1, -1]
                        manvrptemp[:, objk] = yfinal[-1, :3]
                        manvrvtemp[:, objk] = yfinal[-1, 3:6]
    
                    # 存储协方差数据
                    pcovstored[objk].p = pcovnow
                    pcovoffset[objk, :] = np.diag(pcovnow)
                    # 单点生成估计误差以计算临时TLE数据
                    sigr = np.sqrt(pcovoffset[objk, :3])
                    sigv = np.sqrt(pcovoffset[objk, 4:6])
                    # ECI转TEME转换
                    DateConv = DateNow + datetime.timedelta(minutes=TimeToManeuver[objk] + TimeLengthOfManeuver)
                    rtrue, vtrue = rveci2teme(DateConv, fname, manvrptemp[:, objk], manvrvtemp[:, objk])
                    satobj[objk].offset = -tsince  # 重置偏移到当前时间（tsince的负值，因此现在offset + tsince = 0）
                    # 更新对象
                    manvrjday = jday(DateConv.year, DateConv.month, DateConv.day, DateConv.hour, DateConv.minute, DateConv.second)
                    satobj[objk].struc.julian = manvrjday
                    satobj[objk].sattle.julian = manvrjday
                    satobj[objk].struc.epoch = manvrjday - jday(1950, 1, 0, 0, 0, 0)
                    satobj[objk].initialepoch = satobj[objk].sattle.julian - jday(1950, 1, 0, 0, 0, 0)
                    satobj[objk].initialjulian = satobj[objk].sattle.julian
                    # P Q sigr sigv
                    if isSwitch[objk] != 2:
                        print('Performing intermediate PV2TLE estimation process.....')
                        # [xe, Pout, _, _] = kfalgoca(100, rtrue.T, vtrue.T, satobj[objk].struc, sigr, sigv, 1)
                        satobj[objk].struc = sgp4init(72, satobj[objk].struc, satobj[objk].sattle.bstar, xe[2], satobj[objk].initialepoch, xe[4], xe[3], xe[6], xe[1], xe[5])
                        # [_, p, v] = sgp4(satobj[objk].struc, 0.0)
                        print('Completed intermediate PV2TLE estimation process')
                
                elif isSwitch[objk] == 1 and PropEndjday[objk] > jdaynext:
                    # 在机动持续时间非常长的情况下
                    # 机动....
                    maxfuel = objmaxfuel[objk]
                    satdrymass = objdrymass[objk]
                    # 计算0小时UTC的儒略日
                    jdate0 = PropStartjday[objk]
                    # 计算UTC的格林威治恒星时间
                    gst0 = gast(jdate0)
                    pressurepoly = pressureprofile[objk]
                    thrustpoly = thrustprofile[objk]
                    isppoly = ispprofile[objk]
    
                    # 加载协方差数据
                    pcovnow = pcovstored[objk].p
                    TimeLengthOfManeuver[objk] = PropTimeStep * 60
                    MnvrTimeVec = np.linspace(0, TimeLengthOfManeuver, int(np.ceil(TimeLengthOfManeuver)) + 1)
                    for tk in range(len(MnvrTimeVec) - 1):
                        # 对卫星2进行平面变化机动
                        bcoeff = 1.0e-6 * ae[objk] * cd[objk] / satmass[objk]
    
                        if PropMethod[objk] == 1:  # 这是用于正负速度方向的燃烧
                            thrustvar[3] = PropVsign[objk]
                            perr = 0
                            x0 = np.hstack((manvrptemp[:, objk], manvrvtemp[:, objk], satmass[objk]))
                            # [_, yfinal] = ode45(@ceqm1hohconj, [0 MnvrTimeVec[tk+1]-MnvrTimeVec[tk]], x0, options)
                        elif PropMethod[objk] == 2:  # 这是用于平面变化机动
                            vmb = np.cross(manvrptemp[:, objk], manvrvtemp[:, objk])
                            thrustvar[3:6] = vmb / np.linalg.norm(vmb)
                            pmat = np.eye(3)  # 无误差
                            x0 = np.hstack((manvrptemp[:, objk], manvrvtemp[:, objk], satmass[objk]))
                            # [_, yfinal] = ode45(@ceqm1vsigndirectconj, [0 MnvrTimeVec[tk+1]-MnvrTimeVec[tk]], x0, options)
    
                        # 协方差传播
                        dj2dx = j2dfdx(np.hstack((manvrptemp[:, objk], manvrvtemp[:, objk])))
                        dUdx = np.block([[np.zeros((3, 3)), np.eye(3)], [dj2dx, np.zeros((3, 3))]])
                        qmat = 1e-12 * np.eye(3) * np.linalg.norm(manvrvtemp[:, objk]) ** 2
                        gmat = np.block([[np.zeros((3, 3))], [np.eye(3)]])
                        # 模拟推进器误差方差
                        fuelused = 100 - (satmass[objk] - satdrymass) / maxfuel * 100
                        estThrust = cal_poly(cal_poly(fuelused, pressurepoly), thrustpoly)
                        dmat = np.block([[np.zeros((3, 3))], [np.eye(3)]])
                        aTmat = thrustErrRatio * estThrust**2 / (satmass[objk]**2) * np.eye(3)
    
                        p0 = pcovnow.flatten()
                        # [_, xx] = ode45(lambda t, pv: covprop(t, pv, dUdx, gmat, qmat, dmat, aTmat), [0, MnvrTimeVec[tk+1] - MnvrTimeVec[tk]], p0, options)
                        pcovnow = xx[-1].reshape((6, 6))
                        
                        satmass[objk] = yfinal[-1, -1]
                        manvrptemp[:, objk] = yfinal[-1, :3]
                        manvrvtemp[:, objk] = yfinal[-1, 3:6]
    
                    # 存储协方差数据
                    pcovstored[objk].p = pcovnow
                    pcovoffset[objk, :] = np.diag(pcovnow)
                    # ECI转TEME转换以获取交会评估
                    DateConv = DateNow + datetime.timedelta(minutes=TimeLengthOfManeuver)
                    rtrue, vtrue = rveci2teme(DateConv, fname, manvrptemp[:, objk], manvrvtemp[:, objk])
                    p = rtrue.T
                    v = vtrue.T
                    satobj[objk].offset = -tsince  # 重置偏移到当前时间（tsince的负值，因此现在offset + tsince = 0)
    
                elif (PropEndjday[objk] >= jdaynow) and (PropEndjday[objk] <= jdaynext):  # 条件2
                    isSwitch[objk] = 2  # 指示机动结束，但稍后需要处理SGP4 TLE
                    TimeToManeuver[objk] = 0
                    pcovnow = pcovstored[objk].p  # 恢复协方差数据
                    # 计算剩余的机动时间
                    TimeLengthOfManeuver[objk] = (PropEndjday[objk] - jdaynow) * 86400
                    # ...
                    DateConv = DateNow + datetime.timedelta(minutes=TimeLengthOfManeuver[objk])
                    # 创建日历日期字符串
                    cdstr = DateConv.strftime("%Y-%m-%d")
                    # 创建统一时间字符串
                    utstr = DateConv.strftime("%H:%M:%S.%f")
                    print(f'Date: {cdstr}\tTime: {utstr}\t End maneuver process for Obj ID: {ObjSat[objk].CatID}')
                    # 机动....
                    maxfuel = objmaxfuel[objk]
                    satdrymass = objdrymass[objk]
                    # 计算0小时UTC的儒略日
                    jdate0 = jdaynow
                    # 计算UTC的格林威治恒星时间
                    gst0 = gast(jdate0)
                    pressurepoly = pressureprofile[objk]
                    thrustpoly = thrustprofile[objk]
                    isppoly = ispprofile[objk]
    
                    MnvrTimeVec = np.linspace(0, TimeLengthOfManeuver[objk], int(np.ceil(TimeLengthOfManeuver[objk])) + 1)
                    for tk in range(len(MnvrTimeVec) - 1):
                        # 对卫星2进行平面变化机动
                        bcoeff = 1.0e-6 * ae[objk] * cd[objk] / satmass[objk]
                        if PropMethod[objk] == 1:  # 这是用于正负速度方向的燃烧
                            thrustvar[3] = PropVsign[objk]
                            perr = 0
                            x0 = np.hstack((manvrptemp[:, objk], manvrvtemp[:, objk], satmass[objk]))
                            # [_, yfinal] = ode45(@ceqm1hohconj, [0 MnvrTimeVec[tk+1]-MnvrTimeVec[tk]], x0, options)
                        elif PropMethod[objk] == 2:  # 这是用于平面变化机动
                            vmb = np.cross(manvrptemp[:, objk], manvrvtemp[:, objk])
                            thrustvar[3:6] = vmb / np.linalg.norm(vmb)
                            pmat = np.eye(3)  # 无误差
                            x0 = np.hstack((manvrptemp[:, objk], manvrvtemp[:, objk], satmass[objk]))
                            # [_, yfinal] = ode45(@ceqm1vsigndirectconj, [0 MnvrTimeVec[tk+1]-MnvrTimeVec[tk]], x0, options)
    
                        # 协方差传播
                        dj2dx = j2dfdx(yfinal[-1, :6])
                        dUdx = np.block([[np.zeros((3, 3)), np.eye(3)], [dj2dx, np.zeros((3, 3))]])
                        qmat = 1e-12 * np.eye(3) * np.linalg.norm(manvrvtemp[:, objk]) ** 2
                        gmat = np.block([[np.zeros((3, 3))], [np.eye(3)]])
                        # 模拟推进器误差方差
                        fuelused = 100 - (satmass[objk] - satdrymass) / maxfuel * 100
                        estThrust = cal_poly(cal_poly(fuelused, pressurepoly), thrustpoly)
                        dmat = np.block([[np.zeros((3, 3))], [np.eye(3)]])
                        aTmat = thrustErrRatio * estThrust**2 / (satmass[objk]**2) * np.eye(3)
    
                        p0 = pcovnow.flatten()
                        # [_, xx] = ode45(lambda t, pv: covprop(t, pv, dUdx, gmat, qmat, dmat, aTmat), [0, MnvrTimeVec[tk+1] - MnvrTimeVec[tk]], p0, options)
                        pcovnow = xx[-1].reshape((6, 6))
    
                        satmass[objk] = yfinal[-1, 6]
                        manvrptemp[:, objk] = yfinal[-1, :3]
                        manvrvtemp[:, objk] = yfinal[-1, 3:6]
                    # 存储协方差数据
                    pcovstored[objk].p = pcovnow
                    pcovoffset[objk, :] = np.diag(pcovnow)
    
                    # 更新对象日期、时间和偏移量
                    DateConv = DateNow + datetime.timedelta(minutes=TimeLengthOfManeuver[objk])
                    manvrjday = jday(DateConv.year, DateConv.month, DateConv.day, DateConv.hour, DateConv.minute, DateConv.second)
                    satobj[objk].struc.julian = manvrjday
                    satobj[objk].sattle.julian = manvrjday
                    satobj[objk].struc.epoch = manvrjday - jday(1950, 1, 0, 0, 0, 0)
                    satobj[objk].initialepoch = satobj[objk].sattle.julian - jday(1950, 1, 0, 0, 0, 0)
                    satobj[objk].initialjulian = satobj[objk].sattle.julian
    
                else:
                    try:
                        # [_, p, v] = sgp4(satobj[objk].struc, satobj[objk].offset + tsince)
                        pass  # 需要根据实际逻辑实现
                    except:
                        ObjPropCond[objk] = 0
    
            if isSwitch[objk] > 1:
                # 重置为SGP4
                isSwitch[objk] = 0
    
                totalidx = range(startidx, endidx + 1)
                epochidx = 1 - startidx
                pcovnow = pcovstored[objk].p  # 恢复协方差
                # 从ECI转换为TEME并获取误差
                DateConv = DateNow + datetime.timedelta(minutes=TimeToManeuver[objk] + TimeLengthOfManeuver[objk])
                # 更新对象
                manvrjday = jday(DateConv.year, DateConv.month, DateConv.day, DateConv.hour, DateConv.minute, DateConv.second)
                satobj[objk].struc.julian = manvrjday
                satobj[objk].initialepoch = satobj[objk].sattle.julian - jday(1950, 1, 0, 0, 0, 0)
                satobj[objk].initialjulian = satobj[objk].sattle.julian
    
                # 获取TEME位置和速度，以及误差
                rout, vout = rveci2teme(DateConv, fname, manvrptemp[:, objk], manvrvtemp[:, objk])
                rtrue = np.zeros((len(totalidx), 3))
                vtrue = np.zeros((len(totalidx), 3))
                rtrue[epochidx, :] = rout.T
                vtrue[epochidx, :] = vout.T
                sigr[epochidx, :] = np.sqrt(np.diag(pcovnow[:3, :3]))
                sigv[epochidx, :] = np.sqrt(np.diag(pcovnow[4:6, 4:6]))
    
                bcoeff = 1.0e-6 * ae[objk] * cd[objk] / satmass[objk]
                ptemp = rtrue[epochidx, :].T
                vtemp = vtrue[epochidx, :].T
                for tk in range(-1, startidx - 1, -1):
                    x0 = np.hstack((ptemp, vtemp))
                    # [_, xx] = ode45(@j2eqm, [0, -deltastep], x0, options)
                    ptemp = xx[-1, :3].T
                    vtemp = xx[-1, 4:6].T
                    # 从eci转换到teme
                    DateConv = DateNow + datetime.timedelta(minutes=TimeToManeuver[objk] + TimeLengthOfManeuver[objk] + tk * deltastep)
                    rout, vout = rveci2teme(DateConv, fname, ptemp, vtemp)
                    rtrue[epochidx + tk, :] = rout.T
                    vtrue[epochidx + tk, :] = vout.T
                    # 传播协方差以获取估计误差
                    dj2dx = j2dfdx(np.hstack((ptemp, vtemp)))
                    dUdx = np.block([[np.zeros((3, 3)), np.eye(3)], [dj2dx, np.zeros((3, 3))]])
                    qmat = 1e-12 * np.eye(3) * np.linalg.norm(vtemp) ** 12
                    gmat = np.block([[np.zeros((3, 3))], [np.eye(3)]])
                    dmat = np.zeros((6, 3))
                    aTmat = np.zeros((3))
    
                    p0 = pcovnow.flatten()
                    # [_, xx] = ode45(lambda t, pv: covprop(t, pv, dUdx, gmat, qmat, dmat, aTmat), [0, deltastep], p0, options)
                    pcovnow = xx[-1].reshape((6, 6))
                    pcovdiag = np.diag(pcovnow)
                    sigr[epochidx + tk, :] = np.sqrt(pcovdiag[:3])
                    sigv[epochidx + tk, :] = np.sqrt(pcovdiag[4:6])
    
                # 恢复数据以进行正向传播
                ptemp = rtrue[epochidx, :].T
                vtemp = vtrue[epochidx, :].T
                pcovnow = pcovstored[objk].p  # 恢复协方差
                for tk in range(1, endidx + 1):
                    x0 = np.hstack((ptemp, vtemp))
                    # [_, xx] = ode45(@j2eqm, [0, deltastep], x0, options)
                    ptemp = xx[-1, :3].T
                    vtemp = xx[-1, 4:6].T
                    # 从eci转换到teme
                    DateConv = DateNow + datetime.timedelta(minutes=TimeToManeuver[objk] + TimeLengthOfManeuver[objk] + tk * deltastep)
                    rout, vout = rveci2teme(DateConv, fname, ptemp, vtemp)
                    rtrue[epochidx + tk, :] = rout.T
                    vtrue[epochidx + tk, :] = vout.T
                    # 传播协方差以获取估计误差
                    dj2dx = j2dfdx(np.hstack((ptemp, vtemp)))
                    dUdx = np.block([[np.zeros((3, 3)), np.eye(3)], [dj2dx, np.zeros((3, 3))]])
                    qmat = 1e-12 * np.eye(3) * np.linalg.norm(vtemp) ** 12
                    gmat = np.block([[np.zeros((3, 3))], [np.eye(3)]])
                    dmat = np.zeros((6, 3))
                    aTmat = np.zeros((3))
    
                    p0 = pcovnow.flatten()
                    # [_, xx] = ode45(lambda t, pv: covprop(t, pv, dUdx, gmat, qmat, dmat, aTmat), [0, deltastep], p0, options)
                    pcovnow = xx[-1].reshape((6, 6))
                    pcovdiag = np.diag(pcovnow)
                    sigr[epochidx + tk, :] = np.sqrt(pcovdiag[:3])
                    sigv[epochidx + tk, :] = np.sqrt(pcovdiag[4:6])
    
                # -- 注意 --
                # 无需在此处分配协方差偏移量
                # -----
                # 通过卡尔曼滤波器更新对象
                # [xe, Pout, _, _] = kfalgoca(100, rtrue, vtrue, satobj[objk].struc, sigr, sigv, epochidx)
                satobj[objk].struc = sgp4init(72, satobj[objk].struc, satobj[objk].sattle.bstar, xe[2], satobj[objk].initialepoch, xe[4], xe[3], xe[6], xe[1], xe[5])
                # [_, p, v] = sgp4(satobj[objk].struc, 0.0)
                # 因为位置和速度是为了评估jdaynext
                try:
                    # [_, p, v] = sgp4(satobj[objk].struc, satobj[objk].offset + tsince)
                    pass  # 需要根据实际逻辑实现
                except:
                    ObjPropCond[objk] = 0
    
            # 在这里创建UVW旋转矩阵并存储每个对象的所有位置信息，只要它们仍然符合条件
            if ObjPropCond[objk] > 0:
                objpnext[objk, :] = p
                objvnext[objk, :] = v
    
                # 计算uvw（或径向，顺行，横向）矩阵
                or_ = p / np.linalg.norm(p)
                h = np.cross(p, v)
                oh = h / np.linalg.norm(h)
                ot = np.cross(oh, or_)
                ot = ot / np.linalg.norm(ot)
                # satobj[objk].Auvw = np.array([or_, ot, oh]).T
                NextOr[objk, :] = or_
                NextOt[objk, :] = ot
                NextOh[objk, :] = oh
            else:
                objpnow[objk, :] = np.zeros(3)
                objvnow[objk, :] = np.zeros(3)
                objConjFlag[objk] = 3
                ConjFlag[:, objk] = 3 * np.ones(tgtNum)
    
        # ----------------------------------------------------
        # 目标卫星传播
    
        tgtNonComplied = []  # 用于记录任何目标TLE的错误
        for tgtkIdx in range(len(tgtNumVec)):
            tgtk = tgtkIdx
            if TgtPropCond[tgtk] > 0:
                try:
                    # [satout, p, v] = sgp4(sattgt[tgtk].struc, sattgt[tgtk].offset + tsince)
                    if (satout.error == 4) or (satout.error == 6):
                        TgtPropCond[tgtk] = 0
                except:
                    TgtPropCond[tgtk] = 0
                
                if TgtPropCond[tgtk] > 0:
                    tgtpnext[tgtk, :] = p
                    tgtvnext[tgtk, :] = v
                    
                    tgtpmat = np.tile(p, (objNum, 1))
                    tgtvmat = np.tile(v, (objNum, 1))
                    drtemp = tgtpmat - objpnext
                    dvtemp = tgtvmat - objvnext
                    rv = np.sum(drtemp * dvtemp, axis=1)
                    rr = np.sqrt(np.sum(drtemp ** 2, axis=1))
                    NextRelativePx[tgtk, :] = drtemp[:, 0]
                    NextRelativePy[tgtk, :] = drtemp[:, 1]
                    NextRelativePz[tgtk, :] = drtemp[:, 2]
                    NextRelativeVx[tgtk, :] = dvtemp[:, 0]
                    NextRelativeVy[tgtk, :] = dvtemp[:, 1]
                    NextRelativeVz[tgtk, :] = dvtemp[:, 2]
                    NextRange[tgtk, :] = rr
                    NextRangeRate[tgtk, :] = rv / rr
                    
                    objIdx = np.where((CurrentRange[tgtk, :] <= analysisThres) & 
                                      (CurrentRangeRate[tgtk, :] <= 0) & 
                                      (NextRange[tgtk, :] <= analysisThres) & 
                                      (NextRangeRate[tgtk, :] >= 0))[0]
                    if objIdx.size > 0:
                        # ---- 用于预测的公式 --------------
                        a1 = CurrentRange[tgtk, objIdx]
                        b1 = CurrentRangeRate[tgtk, objIdx]
                        a2 = NextRange[tgtk, objIdx]
                        b2 = NextRangeRate[tgtk, objIdx]
                        ProjectedMinTime = (a2 - a1 - b2 * PropTimeStep * 60) / (b1 - b2)  # 计算达到最小距离的预测时间
                        
                        # 计算两个可能的最小距离，理论上应该接近
                        ProjectedDistanceL = a1 + b1 * ProjectedMinTime
                        ProjectedDistanceR = a2 - b2 * (PropTimeStep * 60 - ProjectedMinTime)
                        maxProjectedDistance = np.maximum(ProjectedDistanceL, ProjectedDistanceR)
    
                        conjIdx = np.where(maxProjectedDistance <= ConjRangeThres)[0]
                        
                        if conjIdx.size > 0:
                            for kk in range(len(conjIdx)):
                                objk = objNumVec[objIdx[conjIdx[kk]]]
                                CurrentUVW = np.array([CurrentOr[objk, :], CurrentOt[objk, :], CurrentOh[objk, :]])
                                NextUVW = np.array([NextOr[objk, :], NextOt[objk, :], NextOh[objk, :]])
                                CurrentRelativePos = np.array([RelativePx[tgtk, objk], RelativePy[tgtk, objk], RelativePz[tgtk, objk]])
                                NextRelativePos = np.array([NextRelativePx[tgtk, objk], NextRelativePy[tgtk, objk], NextRelativePz[tgtk, objk]])
                                
                                CurrentRelativeVel = np.array([RelativeVx[tgtk, objk], RelativeVy[tgtk, objk], RelativeVz[tgtk, objk]])
                                NextRelativeVel = np.array([NextRelativeVx[tgtk, objk], NextRelativeVy[tgtk, objk], NextRelativeVz[tgtk, objk]])
                                
                                CurrentUVWPos = CurrentUVW @ CurrentRelativePos
                                NextUVWPos = NextUVW @ NextRelativePos
                                CurrentUVWVel = CurrentUVW @ CurrentRelativeVel
                                NextUVWVel = NextUVW @ NextRelativeVel
                                
                                if ObjPropCond[objk] > 0:
                                    # ----------------------------------------------------
                                    # 检查交会分析条件然后进行分析
                                    tmax = timeVec[tstep]
                                    tmin = timeVec[tstep - 1]
    
                                    x0 = ProjectedMinTime[conjIdx[kk]] / 60 + tmin
    
                                    # compute julian date at 0 hours ut
                                    jdate0 = jdaynow
                                    # compute greenwich sidereal time at 0 hours ut
                                    gst0 = gast(jdate0)
    
                                    tout = myipm(x0, time4min, satobj[objk].struc, sattgt[tgtk].struc, satobj[objk].offset, sattgt[tgtk].offset, tmin, tmax, dxscalar, thres, maxcount)
                                    minDistance, relspeed, ObjSinceEpoch, TgtSinceEpoch, cdstr, utstr, jdateconj = conjunctionOutput(satobj[objk].struc, sattgt[tgtk].struc, satobj[objk].offset, sattgt[tgtk].offset, tout[0], ConjStartJulian, minDisThres)
                                    ConjFlag[tgtk, objk] = 2  # 设置条件为2以避免重复计算
    
                                    if (minDistance <= minDisThres) and (jdateconj >= jdaynow) and (jdateconj <= jdaynext):
                                        # 预测对象的误差协方差
                                        # [_, objconjp, objconjv] = sgp4(satobj[objk].struc, ObjSinceEpoch * 1440)
                                        objsig = computeposcov(objconjp, objconjv, ObjSinceEpoch, pcovoffset[objk, 1:3], leotlecov)
                                        
                                        # 预测目标的误差协方差
                                        # [_, tgtconjp, tgtconjv] = sgp4(sattgt[tgtk].struc, TgtSinceEpoch * 1440)
                                        tgtsig = computeposcov(tgtconjp, tgtconjv, TgtSinceEpoch, np.zeros(3), leotlecov)
                                        
                                        # 误差Cov = objsig + tgtsig; # 合并误差协方差
                                        errorCov = objsig + tgtsig  # 合并误差协方差
                                        
                                        if 'R/B' in sattgt[tgtk].Name:
                                            if 'DEB' in sattgt[tgtk].Name:
                                                RcTgt = RcRBDEB
                                            else:
                                                RcTgt = RcRB
                                        elif 'DEB' in sattgt[tgtk].Name:
                                            RcTgt = RcDEB
                                        else:
                                            RcTgt = RcSat
                                            RcTgt = satnamecheck(RcTgt, sattgt[tgtk].Name, cubesat, RcCube)
    
                                        Rc = RcTgt + RcObj[objk]  # 合并对象半径
                                        # 计算碰撞概率
                                        collisionProb = collisionprobabilitySimpson(objconjp, tgtconjp, objconjv, tgtconjv, Rc, 1.0, errorCov, 20, 20)
                                        print(f'Date: {cdstr}\tTime: {utstr}\tMinimum Distance is {minDistance * 1000:.4f} meters, Obj ID: {ObjSat[objk].CatID}, Tgt ID: {sattgt[tgtk].struc.satnum}, Obj TLE since: {ObjSinceEpoch:.3f} days and Tgt TLE since: {TgtSinceEpoch:.3f} days, Collision Probability (1e-6): {collisionProb * 1e6:.6f}')
    
                                        # Report to csv format *.dat file
                                        ConjunctionReportOutput(1, reportfile, ObjSat[objk].Name, ObjSat[objk].CatID, TgtSat[tgtk].Name, TgtSat[tgtk].CatID, minDistance, relspeed, ObjSinceEpoch, TgtSinceEpoch, cdstr, utstr, collisionProb)
                    
                    else:
                        tgtNonComplied.append(tgtkIdx)
    
        if tgtNonComplied:
            tgtNumVec = np.delete(tgtNumVec, tgtNonComplied)
            NextRelativePx[tgtNonComplied, :] = []
            NextRelativePy[tgtNonComplied, :] = []
            NextRelativePz[tgtNonComplied, :] = []
            NextRelativeVx[tgtNonComplied, :] = []
            NextRelativeVy[tgtNonComplied, :] = []
            NextRelativeVz[tgtNonComplied, :] = []
            ConjFlag[tgtNonComplied, :] = []
            tgtpnext[tgtNonComplied, :] = []
            tgtvnext[tgtNonComplied, :] = []
            NextRangeRate[tgtNonComplied, :] = []
            NextRange[tgtNonComplied, :] = []
    
        # 更新量信息
        objpnow = objpnext
        objvnow = objvnext
    
        tgtpnow = tgtpnext
        tgtvnow = tgtvnext
    
        RelativePx = NextRelativePx
        RelativePy = NextRelativePy
        RelativePz = NextRelativePz
        RelativeVx = NextRelativeVx
        RelativeVy = NextRelativeVy
        RelativeVz = NextRelativeVz
    
        CurrentRangeRate = NextRangeRate
        CurrentRange = NextRange
