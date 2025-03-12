import requests

def remove_blank_lines(file_path):
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

def download_tle():#下载实时卫星数据
    username = 'elegst@nus.edu.sg'
    password = 'thisiswonderfulworld'
    login_url = "https://www.space-track.org/ajaxauth/login"
    api_url = "https://www.space-track.org/basicspacedata/query/class/gp/decay_date/null-val/epoch/%3Enow-30/orderby/norad_cat_id/format/3le"
    payload = {
        'identity': username,
        'password': password
    }

    with requests.Session() as session:
        login_response = session.post(login_url, data=payload)
        if login_response.status_code != 200:
            print("登录失败！请检查用户名和密码。")
            return
        
        print("登录成功，正在下载 TLE 数据...")
        tle_response = session.get(api_url)
        
        if tle_response.status_code == 200:
            with open('temptle.tle', 'w') as f:
                f.write(tle_response.text)
            print("数据下载成功！")
        else:
            print("数据下载失败！请检查 API URL。")
        session.cookies.clear()
    remove_blank_lines('temptle.tle')


if __name__=="__main__":
    download_tle()
    
