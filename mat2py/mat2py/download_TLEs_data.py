import requests

def download_tle():#下载实时卫星数据
    username = 'elegst@nus.edu.sg'
    password = 'thisiswonderfulworld'
    login_url = "https://www.space-track.org/ajaxauth/login"
    api_url = "https://www.space-track.org/basicspacedata/query/class/gp/decay_date/null-val/epoch/%3Enow-30/orderby/norad_cat_id/format/3le"
    tempfile = 'temptle.tle'
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
            with open(tempfile, 'w') as file:
                file.write(tle_response.text)
            print(f"TLE 数据已保存到 {tempfile}")
        else:
            print("数据下载失败！请检查 API URL。")
        session.cookies.clear()

if __name__=="__main__":
    download_tle()
