
def niceticks(ymin, ymax, n, debug=False):
    
    import numpy as np
    
    ymin = np.double(ymin)
    ymax = np.double(ymax)
    n = np.double(n)
    
    ystep  = (ymax-ymin)/n
    
    if debug:
        print('1.')
        print(ymin,ymax,ystep)
    
    Mi = np.power(10.,-np.floor(np.log10(ystep)))
    
    ymin = np.floor(ymin*Mi)/Mi
    ymax = np.ceil(ymax*Mi)/Mi
    ystep  = (ymax-ymin)/n
    
    if debug:
        print('2.')
        print(ymin,ymax,ystep)
    
    if ystep == 0:
        Mi = np.power(10.,-np.floor(np.log10(ymin)))
        ystep = np.round(ymin*Mi)/Mi
    else:
        Mi = np.power(10.,-np.floor(np.log10(ystep)))
        ystep = np.round(ystep*Mi)/Mi
    
    if debug:
        print('3.')
        print(ymin,ymax,ystep)
    
    if (ystep*Mi > 1) and (ystep*Mi <= 3):
        ystep = 2./Mi
    elif (ystep*Mi > 3) and (ystep*Mi < 4):
        ystep = 4./Mi
    elif (ystep*Mi > 4) and (ystep*Mi <= 7):
        ystep = 5./Mi
    elif (ystep*Mi > 7) and (ystep*Mi <= 10):
        ystep = 10./Mi
    
    if debug:
        print('4.')
        print(ymin,ymax,ystep)
    
    ymin = ystep*np.floor(ymin/ystep);
    ymax = ystep*np.ceil(ymax/ystep);
    
    if debug:
        print('5.')
        print(ymin,ymax,ystep)
    
    if ymin < 0:
        ymino = ymin
        ymin = np.arange(0, ymin-ystep, -ystep)
        if ymin[-1] < ymino:
            ymin = ymin[-1]
        else:
            ymin = ymin[-1] - ystep
    
    if debug:
        print('6.')
        print(ymin,ymax,ystep)
    
    ticks = np.arange(ymin, ymax + ystep, ystep)
    
    if len(ticks) > np.round(1.3*n):
      ystep = ystep*2.
    
    if debug:
        print('7.')
        print(ymin,ymax,ystep)
    
    ticks = np.arange(ymin, ymax + ystep, ystep)
    
    if debug:
        print('8.')
        print(ticks)
        
    return ticks
