import numpy as np
from download_TLEs_data import download_tle
import time 
from sgp4.api import jday
import math

tempfile = './temptle.tle'
objfile = './objects.tle'
tgtfile = './targets.tle'


mu = 398600.8 # 地球标准重力参数 (km^3/s^2)
smaThresHold = 50  # 半长轴阈值 (可调整)
ConjStartJulian = 2451545.0  # 2000年1月1日的儒略日



def remove_blank_lines(file_path):#移除原文件空行
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()
        non_empty_lines = [line for line in lines if line.strip()]
        with open(file_path, 'w') as file:
            file.writelines(non_empty_lines)
        print(f"空行已删除，文件 '{file_path}' 已更新。")
    except FileNotFoundError:
        print(f"错误: 文件 '{file_path}' 未找到。")
    except Exception as e:
        print(f"发生错误: {e}")


def process_satellite_tle(tempfile, objfile, ObjCatID):
    objNum = len(ObjCatID)  
    ObjSat = [{} for _ in range(objNum)] 
    IsObj = 0 
    with open(tempfile, 'r') as fid, open(objfile, 'w') as fobj:
        l1 = fid.readline() 
        while l1:
            l2 = fid.readline() 
            l3 = fid.readline() 
            tempCatID = l2[2:7]  
            tempSatName = l1[2:].strip()
            for objIdx in range(objNum):
                ObjSat[objIdx]['CatID'] = ObjCatID[objIdx, 0] 
                if ObjCatID[objIdx, 0] == tempCatID:                  
                    ObjSat[objIdx]['Name'] = tempSatName
                    ObjSat[objIdx]['line1'] = l1.strip()
                    ObjSat[objIdx]['line2'] = l2.strip()
                    ObjSat[objIdx]['line3'] = l3.strip()
                    tempMeanMotion = float(l3[52:63].strip())  
                    tempn = tempMeanMotion * 2 * np.pi / 86400  
                    tempSMA = (mu / tempn**2) ** (1/3)  
                    ObjSat[objIdx]['SMA'] = tempSMA
                    ObjSat[objIdx]['ecc']=float('0.' + l3[26:33])
                    IsObj += 1
                    fobj.write(f"{l1[2:].strip()}\n")
                    fobj.write(f"{l2.strip()}\n")
                    fobj.write(f"{l3.strip()}\n")

            l1 = fid.readline() 
    ObjSat = [sat for sat in ObjSat if 'Name' in sat and sat['Name']]
    return ObjSat, IsObj

# 解析 TLE 数据
def parse_tle(filepath):
    satellites = []
    with open(filepath, 'r') as f:
        while True:
            l1 = f.readline().strip()
            if not l1:
                break
            l2 = f.readline().strip()
            l3 = f.readline().strip()
            if not (l1.startswith('0') and l2.startswith('1') and l3.startswith('2')):
                continue
            cat_id = l2[2:7].strip()
            name = l1[2:].strip()
            ecc = float(f"0.{l3[26:33]}")
            mean_motion = float(l3[52:63].strip())
            # 计算半长轴 (SMA)
            temp_n = mean_motion * 2 * np.pi / 86400 
            temp_sma = (mu / temp_n**2)**(1/3)  
            
            satellites.append({
                "CatID": cat_id,
                "Name": name,
                "line1": l1,
                "line2": l2,
                "line3": l3,
                "eccentricity": ecc,
                "SMA": temp_sma,
            })
    return satellites

# 判断卫星是否符合碰撞分析条件
def qualify_satellite(sat, objsma, smaThresHold):
    temp_sma = sat["SMA"]
    ecc = sat["eccentricity"]
    if (temp_sma <= (max(objsma) + smaThresHold) and 
        temp_sma >= (min(objsma) - smaThresHold)):
        return True
    if (((1 + ecc) * temp_sma >= max(objsma)) and 
        ((1 - ecc) * temp_sma <= min(objsma))):
        return True
    return False

# 生成目标卫星文件targets.tle
def generate_target_file(satellites, objsma):
    targets = []
    with open(tgtfile, 'w') as ftgt:
        for sat in satellites:
            if qualify_satellite(sat, objsma, smaThresHold):
                ftgt.write(f"{sat['line1'][2:]}\n")
                ftgt.write(f"{sat['line2']}\n")
                ftgt.write(f"{sat['line3']}\n")
                targets.append(sat)
    print(f"目标卫星已保存到 {tgtfile}")
    return targets,len(targets)

def setup_TLEfiles():
    if not tempfile:
        download_tle()
    # 解析已下载的 TLE 数据
    remove_blank_lines(tempfile)
    TmpSat = parse_tle(tempfile)
    # 提取对象卫星的数据
    objsma = [sat['SMA'] for sat in TmpSat]
    # 生成目标卫星文件
    starttime=time.time()
    TgtSat,tgtNum = generate_target_file(TmpSat, objsma)
    endtime=time.time()
    cost_time=endtime-starttime
    print(f"处理完成，总共找到 {len(TgtSat)} 个目标卫星用于碰撞分析")
    print(f"共消耗{cost_time}s处理数据")
    return tgtNum,TgtSat





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


