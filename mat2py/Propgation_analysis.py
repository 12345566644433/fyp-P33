import numpy as np
from sgp4.earth_gravity import wgs72
from sgp4.io import twoline2rv
from datetime import datetime, timedelta
import os
import csv
import math
from tools.date_trans import *



def time4min(x, satobj_sattle, sattgt_sattle, satobj_offset, sattgt_offset):
    """计算最小距离的时间函数"""
    t_since = x * 60  
    
    try:
        p_obj, v_obj = satobj_sattle.propagate(0, 0, (satobj_offset + t_since) / (60 * 24)) 
        p_tgt, v_tgt = sattgt_sattle.propagate(0, 0, (sattgt_offset + t_since) / (60 * 24))  
        
        dr = np.array(p_tgt) - np.array(p_obj)
        dv = np.array(v_tgt) - np.array(v_obj)
        
        r = np.linalg.norm(dr)
        if r > 0:
            rdot = np.dot(dr, dv) / r
            return rdot
        else:
            return float('inf')
    except:
        return float('inf')


def myipm(x0, func, satobj_sattle, sattgt_sattle, satobj_offset, sattgt_offset, tmin, tmax):
    """迭代点法求解最小距离时间"""
    dxscalar = 1e-6
    thres = 1e-8
    maxcount = 50
    
    x = x0
    count = 0
    dx = dxscalar
    
    while count < maxcount:
        f = func(x, satobj_sattle, sattgt_sattle, satobj_offset, sattgt_offset)
        
        if abs(f) < thres:
            return x
        fp = (func(x + dx, satobj_sattle, sattgt_sattle, satobj_offset, sattgt_offset) - f) / dx
        
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


def conjunction_output(satobj_sattle, sattgt_sattle, satobj_offset, sattgt_offset, tca, ConjStartJulian):
    """输出会合信息"""
    t_since_days = (satobj_offset + tca * 60) / (60 * 24)  
    t_since_days_tgt = (sattgt_offset + tca * 60) / (60 * 24)  
    
    p_obj, v_obj = satobj_sattle.propagate(0, 0, t_since_days)
    p_tgt, v_tgt = sattgt_sattle.propagate(0, 0, t_since_days_tgt)
    
    dr = np.array(p_tgt) - np.array(p_obj)
    dv = np.array(v_tgt) - np.array(v_obj)
    min_distance = np.linalg.norm(dr)
    rel_speed = np.linalg.norm(dv)
    
    obj_since_epoch = t_since_days
    tgt_since_epoch = t_since_days_tgt
    
    jdate_conj = ConjStartJulian + tca / 1440.0

    days_since_j2000 = jdate_conj - 2451545.0
    dt_conj = datetime(2000, 1, 1) + timedelta(days=days_since_j2000)
    
    cdstr = dt_conj.strftime('%Y-%m-%d')
    utstr = dt_conj.strftime('%H:%M:%S')
    
    return min_distance, rel_speed, obj_since_epoch, tgt_since_epoch, cdstr, utstr, jdate_conj


def compute_collision_probability(p_obj, p_tgt, v_obj, v_tgt, Rc, error_cov):
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
    
    cubesat = []
    try:
        with open('cubesatname_data.txt', 'r') as fr:
            for line in fr:
                name = line.strip().split()[0]
                cubesat.append(name)
    except:
        print("无法读取立方体卫星数据文件")
    
    ConjFlag = np.ones((len(tgtSatDetail.sattgt), len(objSatDetail.satobj)))
    objConjFlag = np.ones(len(objSatDetail.satobj))
    DateTrack = datetime(ConjStartDate.year, ConjStartDate.month, ConjStartDate.day)
    
    ObjPropCond = np.ones(len(objSatDetail.satobj))
    TgtPropCond = np.ones(len(tgtSatDetail.sattgt))
    
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
        
        for objk in range(len(objSatDetail.satobj)):
            if ObjPropCond[objk] > 0:
                try:
                    # 使用SGP4传播对象卫星
                    t_since_days = (objSatDetail.satobj[objk]['offset'] + tsince) / (60 * 24)  
                    p, v = objSatDetail.satobj[objk]['sattle'].propagate(0, 0, t_since_days)
                    
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
                    ObjPropCond[objk] = 0
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
            if TgtPropCond[tgtk] > 0:
                try:
                    # 使用SGP4传播目标卫星
                    t_since_days = (tgtSatDetail.sattgt[tgtk]['offset'] + tsince) / (60 * 24)  
                    p, v = tgtSatDetail.sattgt[tgtk]['sattle'].propagate(0, 0, t_since_days)
                    
                    tgtpnext[tgtk, :] = p
                    tgtvnext[tgtk, :] = v
                    
                
                    for objk in range(len(objSatDetail.satobj)):
                        if ObjPropCond[objk] > 0:
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
                    
                    if len(objIdx) > 0:
                        # 计算预测的最小距离时间
                        a1 = tgtSatDetail.CurrentRange[tgtk, objIdx]
                        b1 = tgtSatDetail.CurrentRangeRate[tgtk, objIdx]
                        a2 = NextRange[tgtk, objIdx]
                        b2 = NextRangeRate[tgtk, objIdx]
                        
                        
                        valid_idx = np.where(b1 != b2)[0]
                        if len(valid_idx) > 0:
                            ProjectedMinTime = np.zeros_like(a1)
                            ProjectedMinTime[valid_idx] = (a2[valid_idx] - a1[valid_idx] - b2[valid_idx] * PropTimeStep * 60) / (b1[valid_idx] - b2[valid_idx])
                            
                            
                            ProjectedDistanceL = a1 + b1 * ProjectedMinTime
                            ProjectedDistanceR = a2 - b2 * (PropTimeStep * 60 - ProjectedMinTime)
                            maxProjectedDistance = np.maximum(ProjectedDistanceL, ProjectedDistanceR)
                            
                            
                            conjIdx = np.where(maxProjectedDistance <= ConjRangeThres)[0]
                            
                            if len(conjIdx) > 0:
                                for kk in range(len(conjIdx)):
                                    objk = objIdx[conjIdx[kk]]
                                    
                                    
                                    if ObjPropCond[objk] > 0:
                                        tmax = timeVec[tstep]
                                        tmin = timeVec[tstep-1]
                                        
                                        x0 = ProjectedMinTime[conjIdx[kk]] / 60 + tmin
                                        
                                        tout = myipm(x0, time4min, 
                                                   objSatDetail.satobj[objk]['sattle'], 
                                                   tgtSatDetail.sattgt[tgtk]['sattle'],
                                                   objSatDetail.satobj[objk]['offset'], 
                                                   tgtSatDetail.sattgt[tgtk]['offset'],
                                                   tmin, tmax)
                                        
                                        # 计算会合输出
                                        min_distance, rel_speed, obj_since_epoch, tgt_since_epoch, cdstr, utstr, jdate_conj = conjunction_output(
                                            objSatDetail.satobj[objk]['sattle'], 
                                            tgtSatDetail.sattgt[tgtk]['sattle'],
                                            objSatDetail.satobj[objk]['offset'], 
                                            tgtSatDetail.sattgt[tgtk]['offset'],
                                            tout, ConjStartJulian)
                                        
                                        ConjFlag[tgtk, objk] = 2 
                                        
                                        
                                        if (min_distance <= minDisThres and jdate_conj >= jdaynow and jdate_conj <= jdaynext):
                                           
                                            t_since_days_obj = obj_since_epoch
                                            t_since_days_tgt = tgt_since_epoch
                                            
                                            p_obj, v_obj = objSatDetail.satobj[objk]['sattle'].propagate(0, 0, t_since_days_obj)
                                            p_tgt, v_tgt = tgtSatDetail.sattgt[tgtk]['sattle'].propagate(0, 0, t_since_days_tgt)
                                            
                                            
                                            if "R/B" in tgtSatDetail.sattgt[tgtk]['Name']:
                                                if "DEB" in tgtSatDetail.sattgt[tgtk]['Name']:
                                                    RcTgt = RcRBDEB
                                                else:
                                                    RcTgt = RcRB
                                            elif "DEB" in tgtSatDetail.sattgt[tgtk]['Name']:
                                                RcTgt = RcDEB
                                            else:
                                                RcTgt = RcSat
                                                RcTgt = satnamecheck(RcTgt, tgtSatDetail.sattgt[tgtk]['Name'], cubesat, RcCube)
                                            
                                            # 合并对象半径
                                            Rc = RcTgt + RcObj[objk]
                                            
                                            error_cov = np.eye(6) 
                                            collision_prob = compute_collision_probability(
                                                p_obj, p_tgt, v_obj, v_tgt, Rc, error_cov)
                                            
                                            print(f'Date: {cdstr}\tTime: {utstr}\tMinimum Distance is {min_distance*1000:.4f} meters, '
                                                  f'Obj ID: {objSatDetail.satobj[objk]["struc"]["satnum"]}, '
                                                  f'Tgt ID: {tgtSatDetail.sattgt[tgtk]["struc"]["satnum"]}, '
                                                  f'Obj TLE since: {obj_since_epoch:.3f} days and '
                                                  f'Tgt TLE since: {tgt_since_epoch:.3f} days, '
                                                  f'Collision Probability: {collision_prob:.6e}')
                                            
                                            with open(reportfile, 'a', newline='') as f:
                                                writer = csv.writer(f)
                                                writer.writerow([cdstr, utstr, 
                                                                objSatDetail.satobj[objk]['Name'], 
                                                                objSatDetail.satobj[objk]['struc']['satnum'], 
                                                                tgtSatDetail.sattgt[tgtk]['Name'], 
                                                                tgtSatDetail.sattgt[tgtk]['struc']['satnum'],
                                                                min_distance, rel_speed, obj_since_epoch, 
                                                                tgt_since_epoch, collision_prob])
                except:
                    TgtPropCond[tgtk] = 0
                    tgtNonComplied.append(tgtk)
        
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


def run_conjunction_assessment(objSatDetail, tgtSatDetail, start_date, duration_days, time_step_minutes, 
                              analysis_threshold, conj_range_threshold, min_dis_threshold, report_file):

    total_minutes = duration_days * 24 * 60
    timeVec = np.arange(0, total_minutes + time_step_minutes, time_step_minutes)
    
    print(f"开始会合评估，总计 {len(timeVec)} 个时间步")
    conjunction_assessment(objSatDetail, tgtSatDetail, timeVec, start_date, time_step_minutes, 
                          analysis_threshold, conj_range_threshold, min_dis_threshold, report_file)

    