#collection of code for processing KuKa data Willatt, UCL 2020-2023

def get_magna_mosaic(date,  #date as ascii string e.g.'20200119'
                  ):
    import fnmatch
    import os
    import pandas as pd
    import numpy as np
    
    #find magnaprobe files
    mfile_path = '/Users/rosie/Documents/mosaic/mosaic_data/magnaprobe_old/magna_xy_new_sd/'
    
    if date == '20191220':
        mf = fnmatch.filter(os.listdir(mfile_path),'*'+'20191219'+'*raw*')
    elif date == '20200124':
        mf = fnmatch.filter(os.listdir(mfile_path),'*'+'20200123'+'*raw*')
    else:
        mf = fnmatch.filter(os.listdir(mfile_path),'*'+date+'*raw*')
    print('date: ', date, 'mf: ', mf)
    
    #make dictionary to put values in 
    mg = {'x': [],
          'y': [],
          'SnowDepth_m': []}

    #loop through magnaprobe files and put in values
    for m in mf:
        mag = pd.read_csv(mfile_path + m, delimiter = ',')
    #                     print(mag)
    #                     ax0.scatter(mag['xc'], mag['yc'], c = mag['DepthCm']/100, vmin = 0, vmax = .8)
        mg['x'].extend(mag['xc'])
        mg['y'].extend(mag['yc'])
        mg['SnowDepth_m'].extend(mag['DepthCm']/100)

    #remove any points where the snow depth is ~ 1.2 m as these are for calibration
    mg_calib = np.where(abs(np.array(mg['SnowDepth_m']) - 1.2) <= .01)[0]
#     print('mg_calib', mg_calib)
    if len(mg_calib > 0):
        for mgc in mg_calib:
            mg['SnowDepth_m'][mgc] = np.nan
            
    return mg



#==================================================================================   



#find the closest data points from a comparison instrument
def find_closest(x_ins, y_ins, #arrays of x and y coords of main instrument
                 comp_x, comp_y, plot = 0): #arrays of x and y coords of instrument to find comparison indices of
    
    import matplotlib.pyplot as plt
    import numpy as np

    
    if plot != 0:
        plt.figure(figsize = (10,10))
        plt.plot(x_ins, y_ins, 'bo')
        plt.plot(comp_x, comp_y, 'm.')
        plt.show()

    if len(x_ins) != len(y_ins):
        print('x_ins and y_ins must have same dimensions!')
        return

    else:
        closest = np.full(len(x_ins),np.nan)
        closest_dist = np.full(len(x_ins),np.nan)

        if plot != 0:
            plt.figure(figsize = (15,15))
        for i in range(len(x_ins)):
#             if x_ins[i] == x_ins[i] and y_ins[i] == y_ins[i] and \
#                         comp_x[i] == comp_x[i] and comp_y[i] == comp_y[i]: #filter out nans
                distance = (((comp_x-x_ins[i])**2 + (comp_y-y_ins[i])**2)**.5)
                closest[i] = np.where(distance == min(distance))[0][0]
                closest_dist[i] = distance[closest[i].astype(int)]

                if plot != 0:
                        plt.plot(x_ins[i],y_ins[i],'c.')
                        plt.plot(comp_x[closest[i].astype(int)],comp_y[closest[i].astype(int)], 'g.')
                        plt.plot((x_ins[i],comp_x[closest[i].astype(int)]),(y_ins[i],comp_y[closest[i].astype(int)]))   
        if plot != 0:
            plt.axis('equal')
            plt.show()
    return closest.astype(int), closest_dist



#==================================================================================   


def remove_static(x, y, cutoff = 0.2): #returns echo indices where the instrument has moved move than cutoff, default 0.2
    import numpy as np
    dist = (x - np.roll(x,1))**2 + (y - np.roll(y,1))**2
    dist[0] = 0 #set to same due to effect of np.roll 
    moving = np.where(dist > cutoff)[0].astype(int)
    return dist, moving




#================================================================================== 

def move_processed(dir_processed, 
                   dir_raw, 
                   mode, #'scan' or 'stare'
                   string_match = '*.nc'):
    
    import fnmatch
    import os
    import numpy as np
    
    print('dir_processed: ',dir_processed)
    print('dir_raw: ',dir_raw)

    processed = np.array(fnmatch.filter(os.listdir(dir_processed), string_match))

#     print('processed:', processed)
    raw = fnmatch.filter(os.listdir(dir_raw),'*.dat')
    raw.sort()
    
    print('processed.shape', processed.shape)
    print('len(raw)', len(raw))
    
    for i in range(len(raw)):
#         print(i)
        r = raw[i]
#         print('r',r)
        to_match = 'kuka_' + mode + '_decon_' + r[0:1] + r[1:2].lower() + r[2:4]+ r[4:7].lower() + r[8:23] + '.nc'
        print(to_match)
        done = np.where(processed == to_match)[0]
        print(done)
        if len(done) > 0:
            print(r, processed[done])
            os.rename(dir_raw + r, dir_raw + '/done/' + r)
        print('---')
    return


#================================================================================== 

#code from Stefan Hendricks as part of floenavi: https://gitlab.awi.de/floenavi-crs/floenavi

def dt64todt(dt64):
    import numpy as np
    from datetime import datetime
    """
    Convert datetime 64 objects to datetime.datetime
    :param dt64:
    :return:
    """
    timestamps = (dt64 - np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's')
    dt = np.ndarray(timestamps.shape, dtype=object)
    for i, ts in enumerate(timestamps):
        dt[i] = datetime.utcfromtimestamp(ts)
    return dt


#================================================================================== 



def areas():
    areas = {'North' :  {'x': [380,830],   'y': [-530,-30]},
             'South' :  {'x': [-660,-170], 'y': [220,620]},
             'Runway':  {'x':[-320,-170],  'y': [-440,70]},
             'Mini':    {'x':[100,290],    'y': [0,150]},
#              'bettrans':{'x':[-200,380],   'y': [-500,240]},
             'Lead':    {'x':[-1200,-750], 'y': [400,620]}}
#             'Lead':    {'x':[-1400,-660], 'y': [400,620]}} #now (above) excluding loops where kuka and magnaprobe not coincident 


    return areas
             

#================================================================================== 


    
def index_areas(x, y):
    
    import numpy as np
    
    a = areas()

    index = {}

    for key in a:
        index[key] = np.where((np.array(x) <= max(a[key]['x'])) & 
                              (np.array(x) >= min(a[key]['x'])) &
                              (np.array(y) <= max(a[key]['y'])) & 
                              (np.array(y) >= min(a[key]['y'])))[0]

        
    return index
             
             
#================================================================================== 