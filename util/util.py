import os
import sys
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf

print('\n\n util path: ~/magnetosphere/util/util.py \n\n')


def tpad(time, length=7):

    # TODO: Check that time is valid
    time = list(time)

    assert(len(time) > 2)
    
    if len(time) > length:
        time = time[0:length]
    else:
        pad = length - len(time)
        time = time + pad*[0]

    return tuple(time)


def tstr(time, length=7):
    """Create date/time string of the convention to tag files with given array of integers
    
    tstr((2000, 1, 1, 2)) # 2000:01:01T02:00:00
    tstr((2000, 1, 1, 2, 3)) # 2000:01:01T02:03:00
    tstr((2000, 1, 1, 2, 3, 4)) # 2000:01:01T02:03:04
    tstr((2000, 1, 1, 2, 3, 4, 567)) # 2000:01:01T02:03:04.567
    """

    # ISO 8601
    assert(len(time) > 2)

    if length == 7:
        return '%04d-%02d-%02dT%02d:%02d:%02d.%03d' % tpad(time, length=length)
    elif length == 6:
        return '%04d-%02d-%02dT%02d:%02d:%02d' % tpad(time, length=length)        
    elif length == 5:
        return '%04d-%02d-%02dT%02d:%02d' % tpad(time, length=length)        
    elif length == 4:
        return '%04d-%02d-%02dT%02d' % tpad(time, length=length)        
    elif length == 3:
        return '%04d-%02d-%02d' % tpad(time, length=length)        

# TODO: maketag(time) -> '_' + tstr(time)
def maketag(time):
    
    return '_%04d:%02d:%02dT%02d:%02d:%02d.%03d' % tpad(time, length=7)


def time2datetime(t):
    import datetime as dt
    
    for i in range(len(t)):
        if int(t[i]) != t[i]:
            raise ValueError("int(t[{0:d}] != t[{0:d}] = {1:f}".format(i, t[i]))\
            
    if len(t) < 3:
        raise ValueError('Time list/tuple must have 3 or more elements')
    if len(t) == 3:
        return dt.datetime(int(t[0]), int(t[1]), int(t[2]))    
    if len(t) == 4:
        return dt.datetime(int(t[0]), int(t[1]), int(t[2]), int(t[3]))    
    if len(t) == 5:
        return dt.datetime(int(t[0]), int(t[1]), int(t[2]), int(t[3]), int(t[4]))    
    if len(t) == 6:
        return dt.datetime(int(t[0]), int(t[1]), int(t[2]), int(t[3]), int(t[4]), int(t[5]))    
    if len(t) == 7:
        return dt.datetime(int(t[0]), int(t[1]), int(t[2]), int(t[3]), int(t[4]), int(t[5]), int(t[6]))    

'''
def filename2time(filename): #TODO: finish
    """Extract time stamp from file name"""
    tstr = filename[11:] 
    y, m, d = int(tstr[0:4]), int(tstr[4:6]), int(tstr[6:8])
    h, M, s = int(tstr[9:11]), int(tstr[11:13]), int(tstr[13:15])
    f = int(tstr[16:19])
    return [y, m, d, h, M, s, f]


def time2filename(time, extension='.out.cdf', split=False):

    filename = '3d__var_3_e' \
        + '%04d%02d%02d-%02d%02d%02d-%03d' % tpad(time, length=7) + extension
    if split:
        return filename
    return conf["run_path"] + filename


def time2SWPCfile(time):
    import numpy as np
    t = np.array(time)
    if len(t.shape) != 1:
        #filenames = np.empty((t.shape[0],), dtype=str)
        filenames = []
        for i in range(t.shape[0]):
            ret = time2SWPCfile(t[i,:])
            if ret != None:
                filenames.append(ret)
        return filenames
    time = tpad(time)
    listnames = conf['SWPC_cdf']+'SWPC_SWMF_052811_2_GM_cdf_list'
    a = np.loadtxt(listnames, dtype=str, skiprows=1) #(N,5)
    Tr = np.logical_and(a[:,2] == '{0:04d}/{1:02d}/{2:02d}'.format(*time[0:3]),
                        a[:,4] == '{0:02d}:{1:02d}:{2:02d}'.format(*time[3:6]))
    if a[Tr, 0].size != 0:
        return conf['SWPC_cdf'] + a[Tr, 0][0]
'''

def dirlist(rootdir, **kwargs):
    """Recursive file list constrained by regular expression
    
    dirlist(rootdir)
    dirlist(rootdir, regex=regex)

    Example:
    -------
    
    from util import dirlist
    dl = dirlist('.', regex='\.py$') # Find files ending in .py
    print(dl)

    import os
    import util

    direct = '/home/gary/magnetosphere/data/DIPTSUR2-derived/cutplanes'

    files = util.dirlist(direct, regex='            ', fullpath=True)

    for f in files:
        os.system('mv %s %s'%(f.replace('            ', '\\ \\ \\ \\ \\ \\ \\ \\ \\ \\ \\ \\ '), f.replace('            ', '')))
    
    """
    import re

    file_keep = []
    for subdir, dirs, files in os.walk(rootdir):
        for file in files:
            if 'regex' in kwargs:
                if re.search(kwargs['regex'], file):
                    if 'fullpath' in kwargs:
                        if kwargs['fullpath']:
                            file_keep.append(subdir+file)
                        else:
                            file_keep.append(file)
                    else:
                        file_keep.append(file)
            else:
                if 'fullpath' in kwargs:
                    if kwargs['fullpath']:
                        file_keep.append(subdir+file)
                    else:
                        file_keep.append(file)
                else:
                    file_keep.append(file)
    return file_keep

'''
def timelist(listtxt='ls-1.txt'):
    times = []
    for file_name in filelist(listtxt = 'ls-1.txt'):
        times.append(filename2time(file_name))
    return times
'''

def safeprep_fileout(outname):
    if not outname[-4] == '.':
        raise ValueError('invalid extension')
    extension = outname[-4:]
    fname = outname
    while os.path.exists(fname):
        fname = fname + '-old' + extension
    if outname != fname:
        os.system('mv %s %s'%(outname,fname))

def safenumpy_tofile(array, outname):
    safeprep_fileout(outname)
    array.tofile(outname)


def get_available_station_times(run, station):
    import numpy as np
    import read_ccmc_datafiles as r_ccmc
    import read_mag_grid_files as rmg

    if run == 'SWPC':
        if type(station) != str:
            raise ValueError('to get data from SWMF output for SWPC run, \
            need station to be a string for valid magnetometer station name')
        data, headers = r_ccmc.getdata([2006, station])
        times = np.array(data[:, 0:6], dtype=int)

    else:
        magfilenames = np.loadtxt(conf[run + '_cdf'] + 'magfilelist.txt', dtype=str)[:,0]
        times = []
        for file_name in list(magfilenames):
            times.append(rmg.mag_grid_file2time(file_name))
        times = np.array(times1, dtype=int)

    return times


def generate_TESTANALYTIC_cdflist():
    f = open(conf['TESTANALYTIC_derived'] + 'cdflist.txt', 'w')
    times = [(2000, 1, 1, h, mi, 0, 0) for h in range(4) for mi in range(60) ]
    for time in times:
        tup = (tstr(time, length=7), ) + time
        f.write('%s.cdf %d %d %d %d %d %d %d\n'%tup)


def get_available_slices(run):
    import numpy as np
    a = np.loadtxt(conf[run+'_derived'] + 'cdflist.txt', dtype=str)
    files = a[:, 0]
    times = np.array(a[:, 1:], dtype=int)
    return files, times


def time2CDFfilename(run, time, split=False, debug=True, ignore_miliseconds=True):
    """
    >>> u.time2CDFfilename('SCARR5',[2003,11,20,7,7,0])
    '/home/gary/magnetosphere/data/SCARR5_GM_IO2/IO2/3d__var_3_e20031120-070700-000.out.cdf'
    >>> u.time2CDFfilename('SCARR5',[2003,11,20,7,7,0,99])
    '/home/gary/magnetosphere/data/SCARR5_GM_IO2/IO2/3d__var_3_e20031120-070700-099.out.cdf'
    >>> u.time2CDFfilename('SWPC',[2006,12,15,7,7,0])
    '/home/gary/magnetosphere/data/SWPC_SWMF_052811_2/GM_CDF/3d__var_1_t00240700_n0290130.out.cdf'
    >>> u.time2CDFfilename('SWPC',[2006,12,15,7,7,0,99])
    '/home/gary/magnetosphere/data/SWPC_SWMF_052811_2/GM_CDF/3d__var_1_t00240700_n0290130.out.cdf'
    >>> u.time2CDFfilename('SWPC',[2003,11,20,7,7,0])
    >>> u.time2CDFfilename('SWPC',[2003,11,20,7,7,0]) == None
    True
    """
    import numpy as np
    t_arr = np.array(time)
    if len(t_arr.shape) != 1:
        filenames = []
        for i in range(t_arr.shape[0]):
            ret = time2CDFfilename(run, t_arr[i,:])
            if ret != None:
                filenames.append(ret)
        return filenames

    time = tpad(time, length=7)

    files, times = get_available_slices(run)
    ts = times[:,5] + 60*times[:,4] + 60*60*times[:,3] + 60*60*24*times[:,2]
    t = time[5] + 60*time[4] + 60*60*time[3] + 60*60*24*time[2]
    if not ignore_miliseconds:
        ts = ts + 0.001*times[:,6]
        t = t + 0.001*time[6]

    ret = np.argwhere(ts == t)
    if np.size(ret) != 0:
        if split:
            return files[ret[0][0]]
        else:
            return conf[run+'_cdf'] + files[ret[0][0]]


    '''

    if run == 'SCARR5':
        filename = '3d__var_3_e' \
            + '%04d%02d%02d-%02d%02d%02d-%03d' % tpad(time, length=7) + '.out.cdf'

    if run == 'SWPC':
        time = tpad(time)
        listnames = conf['SWPC_cdf']+'SWPC_SWMF_052811_2_GM_cdf_list'
        a = np.loadtxt(listnames, dtype=str, skiprows=1) #(N,5)
        Tr = np.logical_and(a[:,2] == '{0:04d}/{1:02d}/{2:02d}'.format(*time[0:3]),
                            a[:,4] == '{0:02d}:{1:02d}:{2:02d}'.format(*time[3:6]))
        if a[Tr, 0].size == 0:
            print('\n WARNING invalid time \n')
            return None

        filename = a[Tr, 0][0]

    if run == 'CARR_IMPULSE':
        if tpad(time, length=6) != (2019,9,2,6,30,0):
            print('\n WARNING invalid time \n')
            return None
        filename = '3d__var_2_e20190902-063000-009.out.cdf'

    dlfile(conf[run + '_cdf'] + filename, debug=debug)

    if split:
        return filename
    return conf[run + '_cdf'] + filename
    '''

def CDFfilename2time(run, filename): #TODO: finish
    """Extract time stamp from file name"""

    import numpy as np

    if not isinstance(filename, str):
        times = []
        for onefilename in filename:
            ret = CDFfilename2time(run, onefilename)
            if ret is not None:
                times.append(ret)
        return times

    files, times = get_available_slices(run)

    ret = np.argwhere(filename == files)

    if np.size(ret) != 0:
        return times[ret[0][0], :]

    '''
    if run == 'SCARR5':
        tstr = filename[11:]
        y, m, d = int(tstr[0:4]), int(tstr[4:6]), int(tstr[6:8])
        h, M, s = int(tstr[9:11]), int(tstr[11:13]), int(tstr[13:15])
        f = int(tstr[16:19])
        return [y, m, d, h, M, s, f]
    if run == 'CARR_IMPULSE':
        if filename == '3d__var_2_e20190902-063000-009.out.cdf':
            return [2019, 9, 2, 6, 30, 0, 9]
        else:
            return None
    '''

'''
def filelist(run):

    if run == 'CARR_IMPULSE':
        return ['3d__var_2_e20190902-063000-009.out.cdf']

    if run == 'SWPC':
        listnames = conf['SWPC_cdf'] + 'SWPC_SWMF_052811_2_GM_cdf_list'
        a = np.loadtxt(listnames, dtype=str, skiprows=1) #(_,5)
        return a[:, 0]

    if run == 'SCARR5':
        listtxt = 'ls-2.txt'
        # Get list of run files
        #urlretrieve(conf['run_url'] + listtxt, conf['run_path'] + listtxt)
        # Read list of run files
        ls = conf[run + '_cdf'] + listtxt
        print('Reading ' + ls)
        with open(ls,'r') as f:
            files = f.readlines()
        print('Read ' + ls)
        
        # Keep only certain files
        i = 0
        while i < len(files):
            files[i] = files[i].rstrip()
            if not files[i][0:9] == '3d__var_3':
                files.pop(i)
            else:
                i = i + 1

        return files
'''


def urlretrieve(url, fname):
    """Python 2/3 urlretrieve compatability function.
    """
    #https://stackoverflow.com/questions/17285464/whats-the-best-way-to-download-file-using-urllib3
    import urllib3
    import shutil
    http = urllib3.PoolManager()

    with http.request('GET', url, preload_content=False) as req:
        if req.status == 200:
            with open(fname,'wb') as fl:
                shutil.copyfileobj(req, fl)
        else:
            raise ConnectionError ('error in dowloading from '+url)
        #req.release_conn() # with statement should do this automattically??


def dlfile(filename, debug=False):

    assert(filename[0] == '/')
    fdir, fname_split = os.path.split(filename)
    fdir = fdir + '/'
    assert(fdir in conf.values())

    if not os.path.exists(filename):
        if fdir == conf['TESTANALYTIC_cdf']:
            if debug: print('writting dummy ' + filename)
            f = open(filename, 'w')
            f.write('this is a dummy with full filename '+ filename)
            f.close()
            return None

        fileurl =  filename.replace(conf['storage'], conf['mag_server_url'], 1)

        if debug:
            print('Downloading ' + fileurl)
            print('to')
        fname_tmp = filename + ".tmp"
        if debug:
            print(fname_tmp)
        # TODO: Catch download error
        ret = urlretrieve(fileurl, fname_tmp)
        if debug:
            print('Downloaded ' + fileurl)
        os.rename(fname_tmp, filename)
        if debug:
            print('Renamed *.tmp')
        return ret
    return None

def dlfile_test(filename, debug=True):

    assert(filename[0] == '/')
    fdir, fname_split = os.path.split(filename)
    fdir = fdir + '/'
    assert(fdir in conf.values())

    if not os.path.exists(filename):
        if fdir == conf['TESTANALYTIC_cdf']:
            if debug: print('writting dummy ' + filename)
            print("f = open(%s, 'w')"%filename)
            print("f.write('this is a dummy with full filename '+ filename")
            print("f.close()")
            return None

        fileurl =  filename.replace(conf['storage'], conf['mag_server_url'], 1)

        if debug:
            print('Downloading ' + fileurl)
            print('to')
        fname_tmp = filename + ".tmp"
        if debug:
            print(fname_tmp)
        print("ret = urlretrieve(fileurl, fname_tmp)")
        ret = 1
        if debug:
            print('Downloaded ' + fileurl)
        print("os.rename(%s, %s)"%(fname_tmp,filename))
        if debug:
            print('Renamed *.tmp')
        return ret
    return None

    assert(filename[0] == '/')
    fdir, fname_split = os.path.split(filename)
    fdir = fdir + '/'
    assert(fdir in conf.values())
    if not os.path.exists(filename):
        print('file doesnt exits')
    return conf['mag_server_url'] + fdir.split('/data/')[-1] + fname_split


def filemeta(filename):
    sys.path.append(conf['interpolator'] + 'kameleon/lib/python2.7/site-packages/ccmc')
    import _CCMC as ccmc

    if type(filename) != str:
        filename = time2filename(filename) #!!!!!! need to update

    #filename = os.path.split(filename)[1]
    #filename = conf['run_path'] + filename

    # Look-up tables for units and labels. Must be extended for each
    # new simulation model.
    units = {            
            'muA/m^2': '$\\mu$A/m$^2$',
            'amu/cm^3': 'amu/cm$^3$',
            'J/m^3': 'J/m$^3$',
            'R': 'R_E'
         }

    names = {
                'x': '$x$',
                'y': '$y$',
                'z': '$z$',
                'bx': '$B_x$',
                'by': '$B_y$',
                'bz': '$B_z$',
                'b1x': '$B_{1x}$',
                'b1y': '$B_{1y}$',
                'b1z': '$B_{1z}$',
                'ux': '$U_x$',
                'uy': '$U_y$',
                'uz': '$U_z$',
                'jx': '$J_x$',
                'jy': '$J_y$',
                'jz': '$J_z$',
                'rho': '$\\rho$',
                'p': '$p$',
                'e': '$e$'
             }
    
    kameleon = ccmc.Kameleon()
    
    if not os.path.exists(filename):
        raise ValueError('Not found: ' + filename)
        return
    
    kameleon.open(filename)

    meta = {}
    meta["filename"] = filename
    meta["parameters"] = {}
    for i in range(kameleon.getNumberOfVariables()):
        varname  = kameleon.getVariableName(i)
        if varname in names:
            meta["parameters"][varname] = {}
            actual_min = kameleon.getVariableAttribute(varname, 'actual_min').getAttributeFloat()
            actual_max = kameleon.getVariableAttribute(varname, 'actual_max').getAttributeFloat()
            native_unit = kameleon.getNativeUnit(varname)
            vis_unit = kameleon.getVisUnit(varname)

            plot_unit = native_unit
            if native_unit in units:
                plot_unit = units[native_unit]

            plot_name = varname
            if plot_name in names:
                plot_name = names[varname]
    
            meta["parameters"][varname] = {
                "actual_min": actual_min,
                "actual_max": actual_max,
                "vis_unit": vis_unit,
                "native_unit": native_unit,
                "plot_unit": plot_unit,
                "plot_name": plot_name
            }


    kameleon.close()
    
    return meta


def printmetadata(filename):
    """Print metadata in a Kameleon .cdf or SWMF .out file
    
    Example
    -------
    import os, sys
    sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
    from config import conf
    from util import printmetadata
    
    printmetadata(conf['run_path'] + '3d__var_3_e20031120-070000-000.out.cdf')
    printmetadata(conf['run_path'] + '3d__var_3_e20031120-070000-000.out')
    """    
    
    if type(filename) != str:
        filename = time2filename(filename)

    if not os.path.exists(filename):
        raise Exception("File not found: " + filename)


    if filename[-3:] == "cdf":
        import _CCMC as ccmc

        print('----------------------------------------------------------------------')
        
        kameleon = ccmc.Kameleon()
        print("Opening " + filename)
        kameleon.open(filename)
        
        for i in range(kameleon.getNumberOfGlobalAttributes()):
            gname = kameleon.getGlobalAttributeName(i)
            gattr = kameleon.getGlobalAttribute(gname)
            if gname != 'README':
                print(gname, gattr.toString())
                
        for i in range(kameleon.getNumberOfVariables()):
            varname  = kameleon.getVariableName(i)
            # Does not seem to be a way to get keys for all variable attributes
            min_attr = kameleon.getVariableAttribute(varname, 'actual_min').getAttributeFloat()
            max_attr = kameleon.getVariableAttribute(varname, 'actual_max').getAttributeFloat()
            units = kameleon.getVisUnit(varname)
            units2 = kameleon.getNativeUnit(varname)
            print(varname, '\t', min_attr,'\t', max_attr, units, units2)
        
        kameleon.close()    
        print('----------------------------------------------------------------------')

    elif filename[-3:] == "out":
        import spacepy.pybats.bats as bats

        print('----------------------------------------------------------------------')
        
        print("Opening " + filename)
        
        data3d = bats.Bats2d(filename)
        
        # look at keys:
        print(data3d.keys())
        
        print(data3d.attrs)
        print(data3d.meta)
        
        print('----------------------------------------------------------------------')
    else:
        raise ValueError("File name does not end in cdf or out. " + filename)

def isclose_32(arr1, arr2):
    if arr1.shape != arr2.shape: return False

    ret = np.empty(arr1.shape ,dtype=bool)
    ret[arr1 == arr2] = True
    tr = arr1 > arr2
    ret[tr] = (1./eps)*(arr1-arr2)/(arr1)
    ret[np.logical_not(tr)] = (1./eps)*(arr1-arr2)/(arr2)

def get_rCurrents(run):  # todo:  get from PARAM file
    if run == 'DIPTSUR2':
        return 1.8
    if run == 'DIPTSUR2_4HOUR':
        return 1.8
    if run == 'UNALT_DIPOLE':
        return 1.8
    if run == 'UNALT_DIPOLE4':
        return 1.8
    if run == 'IMP10_RUN_SAMPLE':
        return 1.7
    if run == 'TESTANALYTIC':
        return 1.5
    if run == 'LUHMANN1979':
        return 1.1
    assert(False)

def get_rBody(run):  # todo:  get from PARAM file
    if run == 'DIPTSUR2':
        return 1.5
    if run == 'DIPTSUR2_4HOUR':
        return 1.5
    if run == 'UNALT_DIPOLE':
        return 1.5
    if run == 'UNALT_DIPOLE4':
        return 1.5
    if run == 'IMP10_RUN_SAMPLE':
        assert(False)
    if run == 'TESTANALYTIC':
        assert(False)
    if run == 'LUHMANN1979':
        return 1.
    assert(False)


def make_data_file(run, time, data_name, data_arr, note=''):
    direct = conf[run+'_derived'] \
              + time2CDFfilename(run, time, split=True)[:-4]+'/'
    if not os.path.exists(direct): os.makedirs(direct)
    arrname = direct + data_name +'.bin'

    if np.isfortran(data_arr):
        order = 'F'
    else:
        order = 'C'
    with open(arrname+'-meta.txt','w') as headerfile:
        headerfile.write('shape = '+str(data_arr.shape)+'\n')
        headerfile.write('arrayorder = '+order+'\n')
        headerfile.write('dtype = '+str(data_arr.dtype)+'\n')
        headerfile.write('byteorder = '+str(data_arr.dtype.byteorder)+'\n')
        headerfile.write(note+'\n')

    data_arr.tofile(arrname)



