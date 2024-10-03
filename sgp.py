import spacetrack.operators as op
from datetime import datetime, date, timedelta
from spacetrack import SpaceTrackClient
from pyorbital.orbital import Orbital
from pyorbital import astronomy
from sgp4.api import Satrec, jday
import json
import math

USERNAME = 'ee3068316@gmail.com'
PASSWORD = 'soskub-0mekme-sEgboj'

satellite_id = [53433, 46, 55, 59, 44713, 44714 , 44715, 44716, 44717, 44718, 44719, 44720, 44721, 44723, 44724]


def get_spacetrack_tle (sat_id, start_date, end_date, username, password, latest=False):
    st = SpaceTrackClient(identity=username, password=password)
    if not latest:
        daterange = op.inclusive_range(start_date, end_date)
        data = st.tle(norad_cat_id=sat_id, orderby='epoch desc', limit=1, format='tle', epoch = daterange)
    else:
        data = st.tle_latest(norad_cat_id=sat_id, orderby='epoch desc', limit=1, format='tle')

    if not data:
        return 0, 0
    
    tle_1 = data[0:69]
    tle_2 = data[70:139]
    incl = data[79:86]
    period = data[122:133]
    print(period)
    return tle_1, tle_2


def create_orbital_track_shapefile_for_day (sat_id, track_day, step_minutes):
    # Для начала получаем TLE    
    # Если запрошенная дата наступит в будущем, то запрашиваем самый последний набор TLE 
    if track_day == date.today():
        tle_1, tle_2 = get_spacetrack_tle(sat_id, None, None, USERNAME, PASSWORD, True)
    # Иначе на конкретный период, формируя запрос для указанной даты и дня после неё
    else:
        tle_1, tle_2 = get_spacetrack_tle(sat_id, track_day, track_day + timedelta(days = 1), USERNAME, PASSWORD, False)

    # Если не получилось добыть    
    if not tle_1 or not tle_2:
        print('Impossible to retrieve TLE')      
        return
    
     # Создаём экземляр класса Orbital
    orb = Orbital("N", line1=tle_1, line2=tle_2)
    # print(f'{tle_1}\n{tle_2}')
    
    # Объявляем счётчики, i для идентификаторов, minutes для времени
    i = 0
    minutes = 0
    output = []
    main_dict = {'id': sat_id, 'loc':[]}
    
    while minutes < 2880:
        minn = timedelta(minutes = minutes)
        

        # utc_hour = int(minutes // 60)
        # utc_days = int(utc_hour // 24)
        # utc_minutes = int((minutes - (utc_hour*60)) // 1)
        # utc_seconds = int(round((minutes - (utc_hour*60) - utc_minutes)*60))

        # utc_time = datetime(track_day.year,track_day.month, track_day.day, utc_hour,utc_minutes ,utc_seconds)
        utc_time = datetime(track_day.year,track_day.month, track_day.day, 0 ,0 ,0) + minn

        jd, fr = jday(utc_time.year,utc_time.month, utc_time.day, utc_time.hour,utc_time.minute,utc_time.second)

        satellite = Satrec.twoline2rv(tle_1, tle_2)
        e, r, v = satellite.sgp4(jd, fr)


        res_dict = {}
        res_tuple = ()
        Cx = r[1]*v[2] - r[2]*v[1]
        Cy = r[2]*v[0] - r[0]*v[2]
        Cz = r[0]*v[1] - r[1]*v[0]
        c = math.sqrt(Cx**2 + Cy**2 + Cz**2)
        r0 = math.sqrt(r[0]**2 + r[1]**2 + r[2]**2)
        mu = 398600.4415
        Lx = -(mu*r[0])/r0 + Cz*v[1] - Cy*v[2]
        Ly = -(mu*r[1])/r0 + Cx*v[2] - Cz*v[0]
        Lz = -(mu*r[2])/r0 + Cy*v[0] - Cx*v[1]
        l = math.sqrt(Lx**2 + Ly**2 + Lz**2)
        e = l/mu
        p = c**2/mu
        a = p/(1-e**2)
        b = a*(math.sqrt(1-(e**2)))
        znach = Cz/c 
        i = math.acos(znach)
        J2 = 1.08263e-3
        n = math.sqrt(mu)/a**(3/2)
        Re = 6378.1
        omega = 3/2*J2*n*((Re/p)**2)*math.cos(i)
        velocity = math.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
        
        
        res_tuple = utc_time.strftime('%Y-%m-%dT%H:%M:%S.%f'), r[0], r[1], r[2], b, a, omega, i, velocity
        print(res_tuple)
        print(utc_time)
        main_dict['loc'].append(res_tuple)
        

        i += 1
        minutes += step_minutes
    return main_dict

final_list = []
for i in satellite_id:
    final_list.append(create_orbital_track_shapefile_for_day(i, date.today(), 5))

with open('res.json', 'w') as f:
    json.dump(final_list, f, indent=4, default=str)
