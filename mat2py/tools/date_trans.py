import math


def calculate_jday(yr, mon, day, hr, min, sec):
    jd = (367.0 * yr
          - math.floor((7 * (yr + math.floor((mon + 9) / 12.0))) * 0.25)
          + math.floor(275 * mon / 9.0)
          + day + 1721013.5
          + ((sec / 60.0 + min) / 60.0 + hr) / 24.0)
    return jd

def date_to_julian(dt):
    yr = dt.year
    mon = dt.month
    day = dt.day
    hr = dt.hour
    minute = dt.minute
    sec = dt.second
    if mon <= 2:
        mon += 12
        yr -= 1
    A = yr // 100
    B = 2 - A + A // 4
    jd = (367.0 * yr
          - math.floor((7 * (yr + math.floor((mon + 9) / 12.0))) * 0.25)
          + math.floor(275 * mon / 9.0)
          + day + 1721013.5
          + ((sec / 60.0 + minute) / 60.0 + hr) / 24.0)
    return jd