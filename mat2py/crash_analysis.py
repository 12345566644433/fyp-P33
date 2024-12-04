import numpy as np
from sgp4.api import Satrec, jday, SatrecArray
from datetime import datetime

# 初始化时间向量
def initialize_time_vector(ConjStartDate, ConjEndDate, PropTimeStep):
    totalelaps = (ConjEndDate - ConjStartDate).total_seconds() / 60  # 转换为分钟
    timeVec = np.arange(0, totalelaps + PropTimeStep, PropTimeStep)
    if timeVec[-1] > totalelaps:
        timeVec[-1] = totalelaps
    return timeVec

# 初始化卫星对象
def initialize_satellites(ObjSat, ConjStartJulian):
    """ 
    参数:
    - ObjSat: 包含卫星 TLE 数据的列表
    - ConjStartJulian: 碰撞分析的开始儒略日期
    返回:
    - sat_objects: 包含所有初始化卫星数据的列表
    - objpnow: 卫星当前位置向量 (km)
    - objvnow: 卫星当前速度向量 (km/s)
    - CurrentOr, CurrentOt, CurrentOh: 分别为径向、轨道和横向单位向量
    """
    objNum = len(ObjSat)
    
    # 初始化存储变量
    objpnow = np.zeros((objNum, 3))
    objvnow = np.zeros((objNum, 3))
    CurrentOr = np.zeros((objNum, 3))
    CurrentOt = np.zeros((objNum, 3))
    CurrentOh = np.zeros((objNum, 3))
    
    sat_objects = []

    # 遍历所有卫星对象，初始化其轨道数据
    for ii in range(objNum):
        # 使用 TLE 数据创建 SGP4 卫星对象
        line2 = ObjSat[ii]['line2']
        line3 = ObjSat[ii]['line3']
        satrec = Satrec.twoline2rv(line2, line3)
        
        # 计算初始儒略日
        initial_julian = satrec.jdsatepoch
        jd, fr = jday(1950, 1, 1, 0, 0, 0)
        initial_epoch = initial_julian - (jd+fr)
        
        # 计算时间偏移量（单位：分钟）
        offset = (ConjStartJulian - initial_julian) * 1440  # 转换为分钟
        # print(f"Initial Julian: {initial_julian}")
        # print(f"Target Julian (jd): {jd}")
        # print(f"Time offset (days): {offset / 1440}")
        # print(f"Conjunction Start Julian: {ConjStartJulian}")
        # 使用 SGP4 传播器计算当前位置和速度
        e, p, v = satrec.sgp4(jd, offset / 1440)  # 将偏移时间转换为天
        if e != 0:
            raise ValueError(f"SGP4 传播失败，错误代码: {e}")
        
        objpnow[ii, :] = p
        objvnow[ii, :] = v

        # 计算 uvw 矩阵（径向、轨道和横向单位向量）
        or_vec = p / np.linalg.norm(p)
        h = np.cross(p, v)  # 角动量向量
        oh_vec = h / np.linalg.norm(h)
        ot_vec = np.cross(oh_vec, or_vec)
        ot_vec = ot_vec / np.linalg.norm(ot_vec)

        CurrentOr[ii, :] = or_vec
        CurrentOt[ii, :] = ot_vec
        CurrentOh[ii, :] = oh_vec

        # 将卫星对象添加到列表中
        sat_objects.append({
            "CatID": ObjSat[ii]['CatID'],
            "satrec": satrec,
            "initial_epoch": initial_epoch,
            "initial_julian": initial_julian,
            "offset": offset,
            "position": p,
            "velocity": v,
            "or_vec": or_vec,
            "ot_vec": ot_vec,
            "oh_vec": oh_vec
        })
    
    return sat_objects, objpnow, objvnow, CurrentOr, CurrentOt, CurrentOh


def compute_relative_positions(TgtSat, objpnow, objvnow, ConjStartJulian):
    """
    计算目标卫星与物体之间的相对位置、速度、径向速度和距离
    """
    tgtNum = len(TgtSat)
    objNum = len(objpnow)
    
    # 初始化存储变量
    tgtpnow = np.zeros((tgtNum, 3))
    tgtvnow = np.zeros((tgtNum, 3))
    RelativePx = np.zeros((tgtNum, objNum))
    RelativePy = np.zeros((tgtNum, objNum))
    RelativePz = np.zeros((tgtNum, objNum))
    RelativeVx = np.zeros((tgtNum, objNum))
    RelativeVy = np.zeros((tgtNum, objNum))
    RelativeVz = np.zeros((tgtNum, objNum))
    CurrentRangeRate = np.zeros((tgtNum, objNum))
    CurrentRange = np.zeros((tgtNum, objNum))
    
    # 遍历每一个目标卫星
    for ii in range(tgtNum):
        # 初始化卫星数据
        sattgt = {}
        sattgt['struc'] = {}
        sattgt['struc']['satnum'] = TgtSat[ii]['CatID']
        sattgt['Name'] = TgtSat[ii]['Name']
        sattgt['sattle'] =  Satrec.twoline2rv(TgtSat[ii]['line2'], TgtSat[ii]['line3'])
        
        # 计算初始纪元时间
        sattgt['initialjulian'] = sattgt['sattle'].jdsatepoch
        jd, fr = jday(1950, 1, 1, 0, 0, 0)
        sattgt['initialepoch'] = sattgt['initialjulian'] - (jd+fr)
        
        # 时间偏移（分钟）
        sattgt['offset'] = (ConjStartJulian - sattgt['initialjulian']) * 1440
             
        # 使用 SGP4 计算当前位置和速度
        tsince = sattgt['offset'] / 1440.0  # 计算时间偏移量
        e, p, v = sattgt['sattle'].sgp4(jd,tsince) 
        tgtpnow[ii, :] = p
        tgtvnow[ii, :] = v
        
        # 计算相对位置和速度
        for kk in range(objNum):
            drtemp = tgtpnow[ii, :] - objpnow[kk, :]
            dvtemp = tgtvnow[ii, :] - objvnow[kk, :]
            rv = np.sum(drtemp * dvtemp)
            rr = np.linalg.norm(drtemp)
            
            # 存储相对位置
            RelativePx[ii, kk] = drtemp[0]
            RelativePy[ii, kk] = drtemp[1]
            RelativePz[ii, kk] = drtemp[2]
            
            # 存储相对速度
            RelativeVx[ii, kk] = dvtemp[0]
            RelativeVy[ii, kk] = dvtemp[1]
            RelativeVz[ii, kk] = dvtemp[2]
            
            # 计算当前的径向速度和距离
            CurrentRangeRate[ii, kk] = rv / rr if rr != 0 else 0
            CurrentRange[ii, kk] = rr

    return RelativePx, RelativePy, RelativePz, RelativeVx, RelativeVy, RelativeVz, CurrentRangeRate, CurrentRange


def compute_obj_relative(objpnow, objvnow):
    """
    计算物体之间的相对距离和相对径向速度
    """
    objNum = len(objpnow)
    
    # 如果物体数量大于 1，才计算相对距离和径向速度
    if objNum > 1:
        objp = objpnow.T  
        objv = objvnow.T  

        objtemprange = []
        objtemprangerate = []

        # 遍历所有物体对
        for obji in range(objNum - 1):
            # 计算相对位置向量
            objdp = objp[:, obji + 1:] - objp[:, [obji]]
            # 计算相对速度向量
            objdv = objv[:, obji + 1:] - objv[:, [obji]]

            # 计算相对距离（欧几里得范数）
            ranges = np.linalg.norm(objdp, axis=0)
            objtemprange.extend(ranges)

            # 计算相对径向速度
            rangerates = np.sum(objdp * objdv, axis=0)
            objtemprangerate.extend(rangerates)

        # 转换为 NumPy 数组
        objtemprange = np.array(objtemprange)
        objtemprangerate = np.array(objtemprangerate)
        
        # 计算当前的距离和径向速度
        objCurrentRange = objtemprange
        objCurrentRangeRate = np.sign(objtemprangerate / objtemprange)
        
        return objCurrentRange, objCurrentRangeRate

    
