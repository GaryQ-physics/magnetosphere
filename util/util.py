import os
import sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf

    
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

            
def filename2time(filename):
    """Extract time stamp from file name"""

    tstr = filename[11:] 
    y, m, d = int(tstr[0:4]), int(tstr[4:6]), int(tstr[6:8])
    h, M, s = int(tstr[9:11]), int(tstr[11:13]), int(tstr[13:15])
    f = int(tstr[16:19])
    return [y, m, d, h, M, s, f]


def time2filename(time):

    filename = conf["run_path"] + '3d__var_3_e' \
        + '%04d%02d%02d-%02d%02d%02d-%03d' % tpad(time, length=7) + '.out.cdf'
        
    return filename


def dirlist(rootdir, **kwargs):
    """Recursive file list constrained by regular expression
    
    dirlist(rootdir)
    dirlist(rootdir, regex=regex)

    Example:
    -------
    
    from util import dirlist
    dl = dirlist('.', regex='\.py$') # Find files ending in .py
    print(dl)
    
    """
    import re

    file_keep = []
    for subdir, dirs, files in os.walk(rootdir):
        for file in files:
            if 'regex' in kwargs:
                if re.search(kwargs['regex'], file):
                    file_keep.append(file)
            else:
                file_keep.append(file)
    return file_keep


def filelist(listtxt = 'ls-1.txt'):

    # Get list of run files
    urlretrieve(conf['run_url'] + listtxt, conf['run_path'] + listtxt)
    
    # Read list of run files
    ls = conf['run_path'] + listtxt
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


def urlretrieve(url, fname):
    """Python 2/3 urlretrieve compatability function.

    If Python 3, returns
    urllib.request.urlretrieve(url, fname)

    If Python 2, returns
    urllib.urlretrieve(url, fname)
    """

    import sys

    if sys.version_info[0] > 2:
        import urllib.request, urllib.error
        try:
            res = urllib.request.urlretrieve(url, fname)
            return res
        except urllib.error.URLError as e:
            print(e)
        except ValueError as e:
            print("'" + url + "' is not a valid URL")
    else:
        import urllib, urllib2, ssl
        try:
            context = ssl._create_unverified_context()
            res = urllib.urlretrieve(url, fname, context=context)
            return res
        except urllib2.URLError as e:
            print(e)
        except ValueError:
            print("'" + url + "' is not a valid URL")

            
def dlfile(filename, debug=False):

    if type(filename) != str:
        filename = time2filename(filename)
    
    filename = os.path.split(filename)[1]
    
    fname_full = conf['run_path'] + filename
    if not os.path.exists(fname_full):
        fileurl = conf['run_url'] + filename
        if debug:
            print('Downloading ' + fileurl)
            print('to')
        fname_tmp = fname_full + ".tmp"
        if debug:
            print(fname_tmp)
        # TODO: Catch download error
        urlretrieve(fileurl, fname_tmp)
        if debug:
            print('Downloaded ' + fileurl)
        os.rename(fname_tmp, fname_full)
        if debug:
            print('Renamed *.cdf.tmp to *.cdf')

    return True


def filemeta(filename):
    
    import _CCMC as ccmc

    if type(filename) != str:
        filename = time2filename(filename)

    filename = os.path.split(filename)[1]
    
    filename = conf['run_path'] + filename

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
                print gname, gattr.toString()
                
        for i in range(kameleon.getNumberOfVariables()):
            varname  = kameleon.getVariableName(i)
            # Does not seem to be a way to get keys for all variable attributes
            min_attr = kameleon.getVariableAttribute(varname, 'actual_min').getAttributeFloat()
            max_attr = kameleon.getVariableAttribute(varname, 'actual_max').getAttributeFloat()
            units = kameleon.getVisUnit(varname)
            units2 = kameleon.getNativeUnit(varname)
            print varname, '\t', min_attr,'\t', max_attr, units, units2
        
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
