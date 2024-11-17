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
        initial_epoch = initial_julian - jday(1950, 1, 1, 0, 0, 0)
        
        # 计算时间偏移量（单位：分钟）
        offset = (ConjStartJulian - initial_julian) * 1440  # 转换为分钟

        # 使用 SGP4 传播器计算当前位置和速度
        e, r, v = satrec.sgp4(offset / 1440)  # 将偏移时间转换为天
        if e != 0:
            raise ValueError(f"SGP4 传播失败，错误代码: {e}")
        
        objpnow[ii, :] = r
        objvnow[ii, :] = v

        # 计算 uvw 矩阵（径向、轨道和横向单位向量）
        or_vec = r / np.linalg.norm(r)
        h = np.cross(r, v)  # 角动量向量
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
            "position": r,
            "velocity": v,
            "or_vec": or_vec,
            "ot_vec": ot_vec,
            "oh_vec": oh_vec
        })
    
    return sat_objects, objpnow, objvnow, CurrentOr, CurrentOt, CurrentOh

# 计算卫星位置和速度
def compute_position_velocity(sat_objects):
    positions = []
    velocities = []
    for sat in sat_objects:
        satrec = sat['satrec']
        e, r, v = satrec.sgp4(sat['offset'] / 1440)  # offset 转换为天
        if e == 0:
            positions.append(np.array(r))
            velocities.append(np.array(v))
        else:
            positions.append(np.zeros(3))
            velocities.append(np.zeros(3))
    return np.array(positions), np.array(velocities)

# 计算相对位置和速度
def compute_relative_metrics(obj_positions, obj_velocities, tgt_positions, tgt_velocities):
    tgtNum, objNum = len(tgt_positions), len(obj_positions)
    RelativePx, RelativePy, RelativePz = np.zeros((tgtNum, objNum)), np.zeros((tgtNum, objNum)), np.zeros((tgtNum, objNum))
    RelativeVx, RelativeVy, RelativeVz = np.zeros((tgtNum, objNum)), np.zeros((tgtNum, objNum)), np.zeros((tgtNum, objNum))
    CurrentRange = np.zeros((tgtNum, objNum))
    CurrentRangeRate = np.zeros((tgtNum, objNum))
    
    for i in range(tgtNum):
        for k in range(objNum):
            dr = tgt_positions[i] - obj_positions[k]
            dv = tgt_velocities[i] - obj_velocities[k]
            rv = np.dot(dr, dv)
            rr = np.linalg.norm(dr)
            
            RelativePx[i, k], RelativePy[i, k], RelativePz[i, k] = dr
            RelativeVx[i, k], RelativeVy[i, k], RelativeVz[i, k] = dv
            CurrentRange[i, k] = rr
            CurrentRangeRate[i, k] = rv / rr if rr != 0 else 0
    return CurrentRange, CurrentRangeRate

# 主函数
def crash_analysis(ObjSat,TgtSat,ConjStartJulian,PropTimeStep):
    # 设置时间范围
    ConjStartDate = datetime(2024, 11, 15)
    ConjEndDate = datetime(2024, 11, 16)
    timeVec = initialize_time_vector(ConjStartDate, ConjEndDate, PropTimeStep)
    
    # 初始化对象和目标卫星
    objSatObjects = initialize_satellites(ObjSat,ConjStartJulian)
    tgtSatObjects = initialize_satellites(TgtSat,ConjStartJulian)
    
    # 计算对象卫星的位置和速度
    obj_positions, obj_velocities = compute_position_velocity(objSatObjects)
    
    # 计算目标卫星的位置和速度
    tgt_positions, tgt_velocities = compute_position_velocity(tgtSatObjects)
    
    # 计算相对位置和速度
    CurrentRange, CurrentRangeRate = compute_relative_metrics(obj_positions, obj_velocities, tgt_positions, tgt_velocities)
    
    print("相对距离 (km):")
    print(CurrentRange)
    print("相对速度变化率 (km/s):")
    print(CurrentRangeRate)

# 初始化向量和矩阵
def init_vector_matrix(objNUm,tgtNum):
    objNumVec = np.arange(1, objNUm + 1)  # 对象卫星的索引向量
    tgtNumVec = np.arange(1, tgtNum + 1)  # 目标卫星的索引向量

    # 传播条件（如果为 0，则不传播，主要用于 Try-Catch 方法）
    ObjPropCond = np.ones(objNUm, dtype=int)
    TgtPropCond = np.ones(tgtNum, dtype=int)

    # 初始化相对距离和速度的矩阵
    NextRange = np.zeros((tgtNum, objNUm))
    NextRangeRate = np.zeros((tgtNum, objNUm))
    CurrentRange = np.zeros((tgtNum, objNUm))
    CurrentRangeRate = np.zeros((tgtNum, objNUm))

    # 碰撞分析标志矩阵
    ConjFlag = np.ones((tgtNum, objNUm), dtype=int)
    objConjFlag = np.ones(objNUm, dtype=int)

    # 初始化碰撞概率和碰撞距离矩阵
    ConjProb = np.zeros((tgtNum, objNUm))
    ConjDistance = np.zeros((tgtNum, objNUm))

    # 初始化当前和下一个时间步长的卫星位置和速度
    tgtpnow = np.zeros((tgtNum, 3))
    tgtvnow = np.zeros((tgtNum, 3))
    objpnow = np.zeros((objNUm, 3))
    objvnow = np.zeros((objNUm, 3))
    tgtpnext = np.zeros((tgtNum, 3))
    tgtvnext = np.zeros((tgtNum, 3))
    objpnext = np.zeros((objNUm, 3))
    objvnext = np.zeros((objNUm, 3))

    # ============================================================
    # 设置对象之间的碰撞分析
    if objNUm > 1:
        icount = 0
        objObjIdxVector = []
        objTgtIdxVector = []

        # 生成对象卫星对之间的索引
        for obji in range(objNUm - 1):
            for objj in range(obji + 1, objNUm):
                objObjIdxVector.append(obji)
                objTgtIdxVector.append(objj)
                icount += 1

        # 将索引转换为 numpy 数组
        objObjIdxVector = np.array(objObjIdxVector)
        objTgtIdxVector = np.array(objTgtIdxVector)

        # 初始化对象之间的相对距离和速度变化率矩阵
        objNextRange = np.zeros(len(objObjIdxVector))
        objCurrentRange = np.zeros(len(objObjIdxVector))
        objNextRangeRate = np.zeros(len(objObjIdxVector))
        objCurrentRangeRate = np.zeros(len(objObjIdxVector))
    return ObjPropCond,TgtPropCond,CurrentRange,ConjFlag,objCurrentRange
if __name__ == "__main__":
    crash_analysis()
