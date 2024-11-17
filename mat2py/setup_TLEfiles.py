import numpy as np
from download_TLEs_data import download_tle
import time 

# 文件路径
tempfile = 'temptle.tle'
objfile = 'objects.tle'
tgtfile = 'targets.tle'

# 常量
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
            temp_n = mean_motion * 2 * np.pi / 86400  # 平均运动 (弧度/秒)
            temp_sma = (mu / temp_n**2)**(1/3)  # 半长轴 (km)
            
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
    ObjSat = parse_tle(tempfile)
    # 提取对象卫星的数据
    objNum = len(ObjSat)
    objsma = [sat["SMA"] for sat in ObjSat]
    # 生成目标卫星文件
    starttime=time.time()
    TgtSat,tgtNum = generate_target_file(ObjSat, objsma)
    endtime=time.time()
    cost_time=endtime-starttime
    print(f"处理完成，总共找到 {len(TgtSat)} 个目标卫星用于碰撞分析")
    print(f"共消耗{cost_time}s处理数据")
    return objNum,tgtNum,ObjSat,TgtSat
if __name__ == "__main__":
    setup_TLEfiles()