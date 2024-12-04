def ConjunctionReportOutput(flag, reportfile, ObjSatName=None, ObjSatCatID=None,
                            TgtSatName=None, TgtSatCatID=None, minDistance=None,
                            relspeed=None, ObjSinceEpoch=None, TgtSinceEpoch=None,
                            cdstr=None, utstr=None, eprob=None):
    #flag: 0 表示创建新文件, 1 表示追加内容
    if flag == 0:
        # 创建文件并写入标题
        with open(reportfile, 'w') as freport:
            freport.write("Date,Time,Object Satellite,Object CatID,Object Days Since Epoch,"
                          "Target Satellite,Target CatID,Target Days Since Epoch,"
                          "Minimum Distance (meter),Relative Speed (km/s),Probability of Collision\n")
        print(f"报告文件 '{reportfile}' 已创建")
    elif flag == 1:
        # 追加分析数据到已有文件
        with open(reportfile, 'a') as freport:
            # 确保所有参数不为空
            if None in [ObjSatName, ObjSatCatID, TgtSatName, TgtSatCatID, 
                        minDistance, relspeed, ObjSinceEpoch, TgtSinceEpoch, 
                        cdstr, utstr, eprob]:
                raise ValueError("追加数据时所有字段均不能为空")
            # 格式化输出内容并写入文件
            freport.write(f"{cdstr},{utstr},{ObjSatName},{ObjSatCatID},{ObjSinceEpoch:.4f},"
                          f"{TgtSatName},{TgtSatCatID},{TgtSinceEpoch:.4f},"
                          f"{minDistance*1000:.4f},{relspeed:.5f},{eprob:.8f}\n")
        print(f"数据已追加到报告文件 '{reportfile}'")
