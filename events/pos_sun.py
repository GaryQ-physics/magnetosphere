# pos_sun

import numpy as np

# units
deg=np.pi/180
amin=deg/60.
hr=1.
minn=hr/60.
s=minn/60.

def GSMtoMAG_matrix(month,day,year,UT):
    a = int((14-month)/12)
    y = year+4800-a
    m = month + 12*a - 3
    t = (UT-12.*hr)/(24.*hr)
    JD = day + int((153*m+2)/5) + y*365 + int(y/4) - int(y/100) + int(y/400) - 32045 + t
    MJD = JD - 2400000.5 #modified julian day
    
    # https://www.spenvis.oma.be/help/background/coortran/coortran.html#Transformations
    T0 = (MJD-51544.5)/36525.0       
    thetaGST = 100.461 + 36000.770*T0+15.04107*UT

    n = JD-2451545.0   #wikipedia
    L = 280.460*deg+0.9856474*deg*n
    g = 357.528*deg+0.9856003*deg*n
    lamb = L+1.915*deg*np.sin(g) + 0.020*deg*np.sin(2*g)
    epsilon = 23.439*deg-0.0000004*deg*n

    X = np.cos(lamb) #pos of sun in GEI
    Y = np.sin(lamb)*np.cos(epsilon)
    Z = np.sin(lamb)*np.sin(epsilon)

    # fgh is pos sun in GEO    [a_dot_x](t)
    f = np.cos(thetaGST)*X + np.sin(thetaGST)*Y    
    g = -np.sin(thetaGST)*X + np.cos(thetaGST)*Y
    h = Z
    # coordinates of magnetic dipole in GEO coordinates (in 2009; assume constant)
    AA = 83*deg+6*amin   # North 
    BB = 117*deg+48*amin   # West 
    m1 = np.cos(BB)*np.cos(AA)   
    m2 = -np.sin(BB)*np.cos(AA)
    m3 = np.sin(AA)

    # a,b,c  x,y,z  u,v,w   are coordinate vector in GEO for the
    # unit vectors for GEO , GSM , MAG respectively
    a = np.array([1,0,0])
    b = np.array([0,1,0])
    c = np.array([0,0,1])
    m = np.array([m1,m2,m3])
    x = np.array([f,g,h]) #to be normalized
    x = (1./np.linalg.norm(x))*x
    w = (1./np.linalg.norm(m))*m
    v = (1./np.linalg.norm(np.cross(c,m)))*np.cross(c,m)
    u = np.cross(v,w)
    y = (1./np.linalg.norm(np.cross(m,x)))*np.cross(m,x)
    z = np.cross(x,y)
    
    # Coordinate transformaton matrix
    A = np.row_stack([u,v,w])
    B = np.row_stack([x,y,z])

    T = np.matmul(A, np.linalg.inv(B))
    return(T)

#function to convert coordinate vector in MAG coordinates to one in GSM coordinates
def MAGtoGSM(v_MAG,month,day,year,UT):
    v_MAG = np.array(v_MAG)

    T = GSMtoMAG_matrix(month,day,year,UT)
    T_inv = np.linalg.inv(T)

    v_GSM = np.matmul(T_inv,v_MAG)
    return(v_GSM)

#function to convert coordinate vector in MAG coordinates to one in GSM coordinates
def GSMtoMAG(v_GSM,month,day,year,UT):
    v_GSM = np.array(v_GSM)

    T = GSMtoMAG_matrix(month,day,year,UT)

    v_MAG = np.matmul(T,v_GSM)
    return(v_MAG)
