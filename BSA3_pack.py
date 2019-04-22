import numpy as np

from astropy.coordinates import Angle
from astropy.time import Time
from astropy import units as u

def read_pnt(filename):
    header = {}
    with open(filename, 'rb') as F: #работает
        for i in range(16):
            line = F.readline()
            a, *b = line.decode("utf-8").strip('\n').split()
            header[a] = b  
        
        data = np.fromfile(F, dtype=np.float32)
        data = data.reshape((int(header['npoints'][0]), 6, 8, len(header['fbands']) + 1))
    
    return header, data

def get_time_begin_and_end(header):
    if 'UTC' in header['date_begin']:
        # begin
        day, month, year = header['date_begin'][2].split('.')
        hour, minute, second = header['time_begin'][2].split(':')
        isot_time = (year + '-' + month + '-' + day + 'T' +
                     hour + ':' + minute + ':' + second 
        )
        time_start = Time(isot_time,
                          format='isot', scale='utc', location=('37.36d', '54.50d'),precision=7)
        # end
        day, month, year = header['date_end'][0].split('.')
        hour, minute, second = header['time_end'][0].split(':')
        isot_time = (year + '-' + month + '-' + day + 'T' +
                     hour + ':' + minute + ':' + second 
        )
        time_end = Time(isot_time,
                        format='isot', scale='utc', location=('37.36d', '54.50d'),precision=7)
    else:
        # begin
        day, month, year = header['date_begin'][0].split('.')
        hour, minute, second = header['time_begin'][0].split(':')
        isot_time = (year + '-' + month + '-' + day + 'T' +
                     hour + ':' + minute + ':' + second 
        )
        time_start = Time(isot_time,
                          format='isot', scale='utc', location=('37.36d', '54.50d'),precision=7)
        time_start -= 4*u.hour
        # end
        day, month, year = header['date_end'][0].split('.')
        hour, minute, second = header['time_end'][0].split(':')
        isot_time = (year + '-' + month + '-' + day + 'T' +
                     hour + ':' + minute + ':' + second 
        )
        time_end = Time(isot_time,
                        format='isot', scale='utc', location=('37.36d', '54.50d'),precision=7)
        time_end -= 4*u.hour
        
    return time_start, time_end

def my_sidereal_time(time):
    t2000 = Time('2000-01-01T00:00:00', format='isot', scale='utc', precision=7)
    t = time.jd - t2000.jd - 1
    t /= 36525
    s0 = 6 + 41 / 60.0 + 50.55 / 3600.0 + 8640184 / 3600.0 * t + 0.093104 / 3600.0 * t * t - 6.27 / 3600.0 * (1e-6) * t * t * t
    t_culm = (time.datetime.hour*u.hour 
          + time.datetime.minute*u.minute
          + time.datetime.second*u.second
          + time.datetime.microsecond*u.microsecond
             )
    alambda = 2 + 30/60.0 + 34/3600.0
    cnst = 2.7379093e-3
    s_culm = s0 + (cnst + 1) * t_culm.value + alambda
    delta_lucha = 0.89
    fi = 0.956829
    be = 0.008436
    delt = delta_lucha
    culm = s_culm
    aa = ((np.sin(fi))**2) * ((np.cos(be))**2) + ((np.cos(fi))**2)
    bb = 2 * np.sin(fi) * np.cos(be) * np.sin(delt)
    cc = ((np.sin(delt))**2) - ((np.cos(fi))**2)

    x = (bb + np.abs(bb) - 4 * aa * cc) / (2 * aa)
    y = x * np.sin(be) / np.cos(delt)
    z = np.sqrt(1 - y * y)

    t = y / z
    alf = culm * np.pi / 12 + np.arctan(t)
    
    rs = 4.8481368e-6
    am = 46.1 * rs
    an = 20.4 * rs

    alfa1 = copy.copy(alf)
    delta1 = copy.copy(delt)
    alfa = copy.copy(alf)

    for i in range(2):
        alf_res = alfa - (am + an * np.sin(alfa1) * np.tan(delta1)) * t
        delt_res = delt - an * np.cos(alfa1) * t
        alfa1 = (alf_res + alfa) / 2
        delta1 = (delt_res + delt) / 2

    alf_res /= (15 * rs);
    delt_res /= rs;

    alf_res /= 3600;
    while alf_res >= 24:
        alf_res -= 24

    while alf_res < 0:
        alf_res += 24

    return Angle(alf_res*u.hourangle)