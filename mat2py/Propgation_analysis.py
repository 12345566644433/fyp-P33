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

def satnamecheck(RcTgt, name, RcCube, cubesat_filename='cubesatname_data.txt'):
    """检查卫星名称是否为立方体卫星"""
    if not os.path.exists(cubesat_filename):
        return RcTgt
    with open(cubesat_filename, 'r') as f:
        cubesat_list = f.read().splitlines()
    for cube_name in cubesat_list:
        if cube_name.lower() in name.lower():
            return RcCube
    return RcTgt


def time4min(x, satobj_sattle, sattgt_sattle, ConjStartJulian):
    """计算最小距离的时间函数"""
    jdate = ConjStartJulian + x / 1440.0
    try:
        jd_day = int(jdate)
    except:
        print(f"Invalid Julian date: {jdate}")
        return float('inf')
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
    pobj_km = np.array(p_obj)
    ptgt_km = np.array(p_tgt)
    relative_dis = np.sum((pobj_km - ptgt_km)**2)
    return relative_dis


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
        return None, None, None, None, None, None, None, None, None, None, None, None
    return min_distance, rel_speed, obj_since_epoch_days, tgt_since_epoch_days, cdstr, utstr, jdate_conj, p_obj, v_obj, p_tgt, v_tgt
    
def compute_pos_cov(pos_km, vel_km_s, time_since_epoch_days, pcov_offset_km2, leotlecov_coeffs):
    pos_km = np.asarray(pos_km).reshape(3)
    vel_km_s = np.asarray(vel_km_s).reshape(3)
    pcov_offset_km2 = np.asarray(pcov_offset_km2).reshape(3)

    if leotlecov_coeffs.shape[0] != 3:
        raise ValueError("leotlecov_coeffs must have 3 rows")

    sig_v = np.polyval(leotlecov_coeffs[0, :], time_since_epoch_days)
    sig_n = np.polyval(leotlecov_coeffs[1, :], time_since_epoch_days)
    sig_r = np.polyval(leotlecov_coeffs[2, :], time_since_epoch_days)

    var_v = sig_v**2
    var_n = sig_n**2
    var_r = sig_r**2

    norm_vel = np.linalg.norm(vel_km_s)
    norm_pos = np.linalg.norm(pos_km)

    if norm_vel < np.finfo(float).eps or norm_pos < np.finfo(float).eps:
        raise ValueError("position and velocity vectors cannot be zero")

    vnc_i = vel_km_s / norm_vel
    h_vec = np.cross(pos_km, vel_km_s)
    norm_h = np.linalg.norm(h_vec)
    if norm_h < np.finfo(float).eps:
         raise ValueError("position and velocity vectors are parallel")
    vnc_n = h_vec / norm_h
    vnc_c = np.cross(vnc_i, vnc_n)

    rot_mat_vnc_to_cartesian = np.vstack((vnc_i, vnc_n, vnc_c)).T
    cov_vnc = np.diag([var_v, var_n, var_r])
    cov_cartesian = rot_mat_vnc_to_cartesian @ cov_vnc @ rot_mat_vnc_to_cartesian.T
    cov_cartesian += np.diag(pcov_offset_km2)

    return cov_cartesian

def calculate_combined_error_covariance(
        pos_obj_km, vel_obj_km_s, time_obj_days, pcov_offset_obj_km2,
        pos_tgt_km, vel_tgt_km_s, time_tgt_days):

    leotlecov_coeffs = np.array([
    [0.000130506662256288, -0.00660267348262973, 0.675201587247618, 0.319582814409902, 4.28641693034372],
    [0.0, 0.0, 0.0, 0.0, 2.82947383703629],
    [0.000289322335966444, -0.00528484813631999, 0.0346960233226241, -0.0102440383785667, 0.240391752770729]])

    cov_obj = compute_pos_cov(pos_obj_km, vel_obj_km_s, time_obj_days,
                              pcov_offset_obj_km2, leotlecov_coeffs)

    pcov_offset_tgt_km2 = np.zeros(3)
    cov_tgt = compute_pos_cov(pos_tgt_km, vel_tgt_km_s, time_tgt_days,
                              pcov_offset_tgt_km2, leotlecov_coeffs)

    error_cov = cov_obj + cov_tgt
    return error_cov



def collision_probability_simpson(pobj, ptgt, vobj, vtgt, Rc, AR, errorCov, hx, hy):
    if errorCov.shape != (3, 3):
         raise ValueError("errorCov must be 3*3")
    pobj = np.asarray(pobj).reshape(3)
    ptgt = np.asarray(ptgt).reshape(3)
    vobj = np.asarray(vobj).reshape(3)
    vtgt = np.asarray(vtgt).reshape(3)

    vr = vobj - vtgt
    vr_norm = np.linalg.norm(vr)

    if vr_norm < np.finfo(float).eps:
         print("warn: relative velocity is approching zero, collision probability is undefined")
         return 0.0

    jk = vr / vr_norm 
    cross_prod_vtgt_vobj = np.cross(vtgt, vobj)
    kk_norm = np.linalg.norm(cross_prod_vtgt_vobj)
    if kk_norm < np.finfo(float).eps:
        print("warning: relative velocity is parallel to the relative position, another method is needed")
        relative_pos_vec = pobj - ptgt
        kk_alt = np.cross(relative_pos_vec, vr)
        kk_alt_norm = np.linalg.norm(kk_alt)
        if kk_alt_norm < np.finfo(float).eps:
             raise ValueError("unable to define the encounter plane")
        kk = kk_alt / kk_alt_norm
    else:
        kk = cross_prod_vtgt_vobj / kk_norm

    ik = np.cross(jk, kk) 
    Mer = np.vstack((ik, jk, kk))

    relative_pos = pobj - ptgt
    rRd = Mer @ relative_pos.reshape(3, 1)

    xm = rRd[0, 0] 
    ym = rRd[2, 0]
    Pcov = Mer @ errorCov @ Mer.T
    epsilon = 1e-15
    sigx_sq = np.abs(Pcov[0, 0])
    sigy_sq = np.abs(Pcov[2, 2])
    sigx = np.sqrt(sigx_sq) + epsilon
    sigy = np.sqrt(sigy_sq) + epsilon

    Lx = Rc   
    Ly = AR * Rc 
    xd, xu = -Lx, Lx  
    yd, yu = -Ly, Ly  
    dx = (xu - xd) / hx 
    dy = (yu - yd) / hy
    nx_points = hx + 1
    ny_points = hy + 1
    nx_vec = np.linspace(xd, xu, nx_points)
    ny_vec = np.linspace(yd, yu, ny_points)
    x_grid, y_grid = np.meshgrid(nx_vec, ny_vec, indexing='ij') 
    exponent = -0.5 * ( ((x_grid - xm) / sigx)**2 + ((y_grid - ym) / sigy)**2 )
    pdf_values = np.exp(exponent)
    wy = np.ones(ny_points)
    wy[1:-1:2] = 4 
    wy[2:-1:2] = 2 
    wx = np.ones(nx_points)
    wx[1:-1:2] = 4 
    wx[2:-1:2] = 2

    inner_integral = np.sum(pdf_values * wy.reshape(1, ny_points), axis=1)
    total_integral_sum = np.sum(wx * inner_integral)
    simpson_coeff = (dx / 3.0) * (dy / 3.0)
    integral_result = simpson_coeff * total_integral_sum
    normalization = 1.0 / (2 * np.pi * sigx * sigy)
    P = normalization * integral_result
    P = max(0.0, min(1.0, P))

    return P

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
        
        
        for tgtk in range(len(tgtSatDetail.sattgt)):
            if not TgtSatStatus[tgtk]:
                continue
            # 使用SGP4传播目标卫星
            tgt_sgp4_status, p, v = tgtSatDetail.sattgt[tgtk]['sattle'].sgp4(jd_day, jd_minnutes)
            if tgt_sgp4_status:
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
                if min_distance is None:
                    continue
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
                        RcTgt = satnamecheck(RcSat, tgtSatDetail.sattgt[tgtk]['Name'], RcCube)

                    
                    # 合并对象半径
                    Rc = RcTgt + RcObj[objk]
                    pcov_offset_km2 = np.array([0.5, 0.2, 0.1])
                    err_cov = calculate_combined_error_covariance(
                        p_obj, v_obj, obj_since_epoch_days, pcov_offset_km2,
                        p_tgt, v_tgt, tgt_since_epoch_days)
                    collision_prob = collision_probability_simpson(p_obj, p_tgt, v_obj, v_tgt, Rc, 1.0, err_cov, 100, 100)
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

    