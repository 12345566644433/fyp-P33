import numpy as np
from sgp4.earth_gravity import wgs72
from sgp4.io import twoline2rv
from datetime import datetime, timedelta
import os
import csv
import math
from tools.date_trans import *

def initialize_time_vector(ConjStartDate: datetime, ConjEndDate: datetime, PropTimeStep: float):
    totalelaps = (ConjEndDate - ConjStartDate).total_seconds() / 60
    timeVec = np.arange(0, np.ceil(totalelaps / PropTimeStep) * PropTimeStep + PropTimeStep, PropTimeStep)
    if timeVec[-1] > totalelaps:
        timeVec[-1] = totalelaps
    return timeVec

def time4min(x, satobj_sattle, sattgt_sattle, ConjStartJulian):
    """计算最小距离的时间函数"""
    jdate = ConjStartJulian + x / 1440.0
    jd_day = int(jdate)
    jd_fraction = jdate - jd_day
    if jd_fraction >= 1.0:
            jd_day += 1
            jd_fraction -= 1.0
    elif jd_fraction < 0.0:
            jd_day -= 1
            jd_fraction += 1.0

    error_code_obj, p_obj, v_obj = satobj_sattle.sgp4(jd_day, jd_fraction)
    if error_code_obj != 0:
        return float('inf')
    error_code_tgt, p_tgt, v_tgt = sattgt_sattle.sgp4(jd_day, jd_fraction)
    if error_code_tgt != 0:
        return float('inf')
    dr = np.array(p_tgt) - np.array(p_obj)
    dv = np.array(v_tgt) - np.array(v_obj)
    dot_product = np.dot(dr, dv)
    return dot_product


def myipm(x0, func, satobj_sattle, sattgt_sattle, tmin, tmax, ConjStartJulian):
    """迭代点法求解最小距离时间"""
    dxscalar = 1e-6
    thres = 1e-8
    maxcount = 50
    x = x0
    count = 0
    dx = dxscalar
    while count < maxcount:
        f = func(x, satobj_sattle, sattgt_sattle,ConjStartJulian)
        
        if abs(f) < thres:
            return x
        fp = (func(x + dx, satobj_sattle, sattgt_sattle, ConjStartJulian) - f) / dx
        
        if fp == 0:
            dx = dx / 2
            continue
        xnew = x - f / fp
        
        if xnew < tmin:
            xnew = tmin
        elif xnew > tmax:
            xnew = tmax
        
        if abs(xnew - x) < thres:
            return xnew
        
        x = xnew
        count += 1
    
    return x

def conjunction_output(satobj_sattle, sattgt_sattle, tca, ConjStartJulian):
    jdate_conj = ConjStartJulian + tca / 1440.0
    jd_day, jd_fraction = math.modf(jdate_conj)
    error_code_obj, p_obj, v_obj = satobj_sattle.sgp4(jd_day, jd_fraction)
    if error_code_obj:
        print(f"Error in SGP4 propagation for object satellite: {error_code_obj}")
        raise ValueError(f"SGP4 propagation error for object satellite: {error_code_obj}")
    error_code_tgt, p_tgt, v_tgt = sattgt_sattle.sgp4(jd_day, jd_fraction)
    if error_code_tgt:
        print(f"Error in SGP4 propagation for target satellite: {error_code_tgt}")
        raise ValueError(f"SGP4 propagation error for target satellite: {error_code_tgt}")
    p_obj_np = np.array(p_obj)
    v_obj_np = np.array(v_obj)
    p_tgt_np = np.array(p_tgt)
    v_tgt_np = np.array(v_tgt)
    dr = p_tgt_np - p_obj_np
    dv = v_tgt_np - v_obj_np
    min_distance = np.linalg.norm(dr)
    rel_speed = np.linalg.norm(dv)
    obj_epoch_jd = satobj_sattle.jdsatepoch + satobj_sattle.jdsatepochF
    tgt_epoch_jd = sattgt_sattle.jdsatepoch + sattgt_sattle.jdsatepochF
    obj_since_epoch_days = jdate_conj - obj_epoch_jd
    tgt_since_epoch_days = jdate_conj - tgt_epoch_jd
    try:
        year, mon, day, hr, minute, sec = invjday(jdate_conj)
        cdstr = f"{int(year):04d}-{int(mon):02d}-{int(day):02d}"
        sec_whole = int(sec)
        sec_frac = sec - sec_whole
        if sec_frac < 0: sec_frac = 0.0
        utstr = f"{int(hr):02d}:{int(minute):02d}:{sec_whole:02d}.{int(sec_frac*1000):03d}"
    except Exception as e:
        print(f"Error converting Julian Date {jdate_conj} to calendar date/time: {e}")
        raise e
    return min_distance, rel_speed, obj_since_epoch_days, tgt_since_epoch_days, cdstr, utstr, jdate_conj, p_obj, v_obj, p_tgt, v_tgt
    


def compute_collision_probability(p_obj, p_tgt, v_obj, v_tgt, Rc):
    dr = np.array(p_tgt) - np.array(p_obj)
    dv = np.array(v_tgt) - np.array(v_obj)
    
    # 计算碰撞平面
    rel_speed = np.linalg.norm(dv)
    if rel_speed < 1e-10:
        return 0
    
    v_unit = dv / rel_speed
    if abs(v_unit[0]) < abs(v_unit[1]):
        u1 = np.array([0, -v_unit[2], v_unit[1]])
    else:
        u1 = np.array([-v_unit[2], 0, v_unit[0]])
    
    u1 = u1 / np.linalg.norm(u1)
    u2 = np.cross(v_unit, u1)
    
    miss_vector = dr - np.dot(dr, v_unit) * v_unit
    miss_distance = np.linalg.norm(miss_vector)
    
    if miss_distance > 5 * Rc:
        return 0
    
    position_uncertainty = 0.1
    
    collision_prob = math.exp(-miss_distance**2 / (2 * position_uncertainty**2))
    collision_prob = min(1.0, collision_prob) 
    
    return collision_prob


def satnamecheck(RcTgt, name, cubesat_list, RcCube):
    """检查卫星名称是否为立方体卫星"""
    for cube_name in cubesat_list:
        if cube_name.lower() in name.lower():
            return RcCube
    return RcTgt


def conjunction_assessment(objSatDetail, tgtSatDetail, timeVec, ConjStartDate, PropTimeStep, analysisThres, ConjRangeThres, minDisThres, reportfile):
    """主会合评估函数"""
    ConjStartJulian = calculate_jday(ConjStartDate.year, ConjStartDate.month, ConjStartDate.day, 
                                    ConjStartDate.hour, ConjStartDate.minute, ConjStartDate.second)
    
    RcSat = 7.5e-3
    RcCube = 1.5e-3
    RcRB = 25e-3
    RcRBDEB = 15e-3
    RcDEB = 10e-3
    RcObj = np.ones(len(objSatDetail.satobj)) * 55e-3
    
    
    ConjFlag = np.ones((len(tgtSatDetail.sattgt), len(objSatDetail.satobj)))
    objConjFlag = np.ones(len(objSatDetail.satobj))
    DateTrack = datetime(ConjStartDate.year, ConjStartDate.month, ConjStartDate.day)
    
    ObjSatStatus = np.ones(len(objSatDetail.satobj))
    TgtSatStatus = np.ones(len(tgtSatDetail.sattgt))
    
    if not os.path.exists(reportfile):
        with open(reportfile, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['Date', 'Time', 'Object Name', 'Object ID', 'Target Name', 'Target ID', 
                            'Min Distance (km)', 'Relative Speed (km/s)', 'Object Since Epoch (days)', 
                            'Target Since Epoch (days)', 'Collision Probability'])
    
    for tstep in range(1, len(timeVec)):
        tsince = timeVec[tstep]
        
        DateNow = ConjStartDate + timedelta(minutes=int(timeVec[tstep-1]))
        DateNext = ConjStartDate + timedelta(minutes=int(timeVec[tstep]))
        jdaynow = calculate_jday(DateNow.year, DateNow.month, DateNow.day, DateNow.hour, DateNow.minute, DateNow.second)
        jdaynext = calculate_jday(DateNext.year, DateNext.month, DateNext.day, DateNext.hour, DateNext.minute, DateNext.second)
        
        if (DateNow - DateTrack).days >= 1.0:
            DateTrack = datetime(DateNow.year, DateNow.month, DateNow.day)
            print(f'The Conjunction Assessment Process Now at {DateTrack.strftime("%d-%b-%Y")}')
        
        objpnext = np.zeros_like(objSatDetail.objpnow)
        objvnext = np.zeros_like(objSatDetail.objvnow)
        NextOr = np.zeros_like(objSatDetail.CurrentOr)
        NextOt = np.zeros_like(objSatDetail.CurrentOt)
        NextOh = np.zeros_like(objSatDetail.CurrentOh)
        jd_target = ConjStartJulian + tsince / 1440.0
        jd_minnutes, jd_day = math.modf(jd_target)
        for objk in range(len(objSatDetail.satobj)):
            if not ObjSatStatus[objk]:
                continue
            try:
                # 使用SGP4传播目标卫星
                obj_sgp4_status, p, v = objSatDetail.satobj[objk]['sattle'].sgp4(jd_day, jd_minnutes)
                if obj_sgp4_status:
                    print(f"Error in SGP4 propagation for object satellite {objk}")
                    continue
                objpnext[objk, :] = p
                objvnext[objk, :] = v
                
                
                r_norm = np.linalg.norm(p)
                if r_norm > 0:
                    or_ = p / r_norm
                    h = np.cross(p, v)
                    h_norm = np.linalg.norm(h)
                    if h_norm > 0:
                        oh = h / h_norm
                        ot = np.cross(oh, or_)
                        ot_norm = np.linalg.norm(ot)
                        if ot_norm > 0:
                            ot = ot / ot_norm
                            NextOr[objk, :] = or_
                            NextOt[objk, :] = ot
                            NextOh[objk, :] = oh
            except:
                ObjSatStatus[objk] = 0
                objpnext[objk, :] = np.zeros(3)
                objvnext[objk, :] = np.zeros(3)
                objConjFlag[objk] = 3
                ConjFlag[:, objk] = 3
        
       
        tgtpnext = np.zeros_like(tgtSatDetail.tgtpnow)
        tgtvnext = np.zeros_like(tgtSatDetail.tgtvnow)
        
        NextRelativePx = np.zeros_like(tgtSatDetail.RelativePx)
        NextRelativePy = np.zeros_like(tgtSatDetail.RelativePy)
        NextRelativePz = np.zeros_like(tgtSatDetail.RelativePz)
        NextRelativeVx = np.zeros_like(tgtSatDetail.RelativeVx)
        NextRelativeVy = np.zeros_like(tgtSatDetail.RelativeVy)
        NextRelativeVz = np.zeros_like(tgtSatDetail.RelativeVz)
        NextRange = np.zeros_like(tgtSatDetail.CurrentRange)
        NextRangeRate = np.zeros_like(tgtSatDetail.CurrentRangeRate)
        
        tgtNonComplied = []
        
        for tgtk in range(len(tgtSatDetail.sattgt)):
            if not TgtSatStatus[tgtk]:
                continue
            # try:
                # 使用SGP4传播目标卫星
            tgt_sgp4_status, p, v = tgtSatDetail.sattgt[tgtk]['sattle'].sgp4(jd_day, jd_minnutes)
            if tgt_sgp4_status:
                print(f"Error in SGP4 propagation for target satellite {tgtk}")
                continue
            tgtpnext[tgtk, :] = p
            tgtvnext[tgtk, :] = v
            
        
            for objk in range(len(objSatDetail.satobj)):
                if ObjSatStatus[objk] > 0:
                    drtemp = p - objpnext[objk, :]
                    dvtemp = v - objvnext[objk, :]
                    rr = np.linalg.norm(drtemp)
                    if rr > 0:
                        rv = np.sum(drtemp * dvtemp)
                        
                        NextRelativePx[tgtk, objk] = drtemp[0]
                        NextRelativePy[tgtk, objk] = drtemp[1]
                        NextRelativePz[tgtk, objk] = drtemp[2]
                        NextRelativeVx[tgtk, objk] = dvtemp[0]
                        NextRelativeVy[tgtk, objk] = dvtemp[1]
                        NextRelativeVz[tgtk, objk] = dvtemp[2]
                        NextRange[tgtk, objk] = rr
                        NextRangeRate[tgtk, objk] = rv / rr
            
            # 检查会合条件
            objIdx = np.where((tgtSatDetail.CurrentRange[tgtk, :] <= analysisThres) & 
                                (tgtSatDetail.CurrentRangeRate[tgtk, :] <= 0) & 
                                (NextRange[tgtk, :] <= analysisThres) & 
                                (NextRangeRate[tgtk, :] >= 0))[0]
            
            if len(objIdx) == 0:
                continue
                # 计算预测的最小距离时间
            a1 = tgtSatDetail.CurrentRange[tgtk, objIdx]
            b1 = tgtSatDetail.CurrentRangeRate[tgtk, objIdx]
            a2 = NextRange[tgtk, objIdx]
            b2 = NextRangeRate[tgtk, objIdx]
            
            
            valid_idx = np.where(b1 != b2)[0]
            if len(valid_idx) == 0:
                continue
            ProjectedMinTime = np.zeros_like(a1)
            ProjectedMinTime[valid_idx] = (a2[valid_idx] - a1[valid_idx] - b2[valid_idx] * PropTimeStep * 60) / (b1[valid_idx] - b2[valid_idx])
            
            
            ProjectedDistanceL = a1 + b1 * ProjectedMinTime
            ProjectedDistanceR = a2 - b2 * (PropTimeStep * 60 - ProjectedMinTime)
            maxProjectedDistance = np.maximum(ProjectedDistanceL, ProjectedDistanceR)
            
            
            conjIdx = np.where(maxProjectedDistance <= ConjRangeThres)[0]
            
            if len(conjIdx) == 0:
                continue
            for kk in range(len(conjIdx)):
                objk = objIdx[conjIdx[kk]]
                
                
                if not ObjSatStatus[objk]:
                    continue
                tmax = timeVec[tstep]
                tmin = timeVec[tstep-1]
                
                x0 = ProjectedMinTime[conjIdx[kk]] / 60 + tmin
                
                tout = myipm(x0, time4min,
                            objSatDetail.satobj[objk]['sattle'], 
                            tgtSatDetail.sattgt[tgtk]['sattle'],
                            tmin, tmax, ConjStartJulian)
                
                # 计算会合输出
                min_distance, rel_speed, obj_since_epoch_days, tgt_since_epoch_days, cdstr, utstr, jdate_conj, p_obj, v_obj, p_tgt, v_tgt = conjunction_output(
                    objSatDetail.satobj[objk]['sattle'], 
                    tgtSatDetail.sattgt[tgtk]['sattle'],
                    tout, ConjStartJulian)
                    
                
                ConjFlag[tgtk, objk] = 2 
                
                
                if (min_distance <= minDisThres and jdate_conj >= jdaynow and jdate_conj <= jdaynext):
                    if "R/B" in tgtSatDetail.sattgt[tgtk]['Name']:
                        if "DEB" in tgtSatDetail.sattgt[tgtk]['Name']:
                            RcTgt = RcRBDEB
                        else:
                            RcTgt = RcRB
                    elif "DEB" in tgtSatDetail.sattgt[tgtk]['Name']:
                        RcTgt = RcDEB
                    else:
                        RcTgt = RcSat
                    
                    # 合并对象半径
                    Rc = RcTgt + RcObj[objk]
                    
                    collision_prob = compute_collision_probability(
                        p_obj, p_tgt, v_obj, v_tgt, Rc)
                    if collision_prob == 0:
                        continue
                    print(f'Date: {cdstr}\tTime: {utstr}\tMinimum Distance is {min_distance*1000:.4f} meters, '
                            f'Obj ID: {objSatDetail.satobj[objk]["struc"]["satnum"]}, '
                            f'Tgt ID: {tgtSatDetail.sattgt[tgtk]["struc"]["satnum"]}, '
                            f'Obj TLE since: {obj_since_epoch_days:.3f} days and '
                            f'Tgt TLE since: {tgt_since_epoch_days:.3f} days, '
                            f'Collision Probability: {collision_prob:.6e}')
                    
                    with open(reportfile, 'a', newline='') as f:
                        writer = csv.writer(f)
                        writer.writerow([cdstr, utstr, 
                                        objSatDetail.satobj[objk]['Name'], 
                                        objSatDetail.satobj[objk]['struc']['satnum'], 
                                        tgtSatDetail.sattgt[tgtk]['Name'], 
                                        tgtSatDetail.sattgt[tgtk]['struc']['satnum'],
                                        min_distance, rel_speed, obj_since_epoch_days, 
                                        tgt_since_epoch_days, collision_prob])
            # except:
            #     TgtSatStatus[tgtk] = 0
            #     tgtNonComplied.append(tgtk)
        
        objSatDetail.objpnow = objpnext
        objSatDetail.objvnow = objvnext
        objSatDetail.CurrentOr = NextOr
        objSatDetail.CurrentOt = NextOt
        objSatDetail.CurrentOh = NextOh
        
        tgtSatDetail.tgtpnow = tgtpnext
        tgtSatDetail.tgtvnow = tgtvnext
        tgtSatDetail.RelativePx = NextRelativePx
        tgtSatDetail.RelativePy = NextRelativePy
        tgtSatDetail.RelativePz = NextRelativePz
        tgtSatDetail.RelativeVx = NextRelativeVx
        tgtSatDetail.RelativeVy = NextRelativeVy
        tgtSatDetail.RelativeVz = NextRelativeVz
        tgtSatDetail.CurrentRangeRate = NextRangeRate
        tgtSatDetail.CurrentRange = NextRange

    