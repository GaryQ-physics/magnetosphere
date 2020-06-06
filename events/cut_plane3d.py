from mpl_toolkits import mplot3d
from matplotlib.patches import Circle, PathPatch, Rectangle
from matplotlib.text import TextPath
from matplotlib.transforms import Affine2D
# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import mpl_toolkits.mplot3d.art3d as art3d

# functions from https://stackoverflow.com/questions/18228966/how-can-matplotlib-2d-patches-be-transformed-to-3d-with-arbitrary-normals
def rotation_matrix(d):
    """
    Calculates a rotation matrix given a vector d. The direction of d
    corresponds to the rotation axis. The length of d corresponds to 
    the sin of the angle of rotation.

    Variant of: http://mail.scipy.org/pipermail/numpy-discussion/2009-March/040806.html
    """
    sin_angle = np.linalg.norm(d)

    if sin_angle == 0:
        return np.identity(3)

    d /= sin_angle

    eye = np.eye(3)
    ddt = np.outer(d, d)
    skew = np.array([[    0,  d[2],  -d[1]],
                  [-d[2],     0,  d[0]],
                  [d[1], -d[0],    0]], dtype=np.float64)

    M = ddt + np.sqrt(1 - sin_angle**2) * (eye - ddt) + sin_angle * skew
    return M

def pathpatch_2d_to_3d(pathpatch, z = 0, normal = 'z'):
    """
    Transforms a 2D Patch to a 3D patch using the given normal vector.

    The patch is projected into they XY plane, rotated about the origin
    and finally translated by z.
    """
    if type(normal) is str: #Translate strings to normal vectors
        index = "xyz".index(normal)
        normal = np.roll((1.0,0,0), index)

    normal = (1./np.linalg.norm(normal))*normal #Make sure the vector is normalised

    path = pathpatch.get_path() #Get the path and the associated transform
    trans = pathpatch.get_patch_transform()

    path = trans.transform_path(path) #Apply the transform

    pathpatch.__class__ = art3d.PathPatch3D #Change the class
    pathpatch._code3d = path.codes #Copy the codes
    pathpatch._facecolor3d = pathpatch.get_facecolor #Get the face color    

    verts = path.vertices #Get the vertices in 2D

    d = np.cross(normal, (0, 0, 1)) #Obtain the rotation vector   
    M = rotation_matrix(d) #Get the rotation matrix

    pathpatch._segment3d = np.array([np.dot(M, (x, y, 0)) + (0, 0, z) for x, y in verts])

def pathpatch_translate(pathpatch, delta):
    """
    Translates the 3D pathpatch by the amount delta.
    """
    pathpatch._segment3d += delta

# Plot field lines
if plot3d:
    for i in range(Nb+1): 
        from_list=solns_restr[i]
        sol=np.array(from_list)
        if (i == 0):
            # Event field line
            solCut=np.zeros((sol.shape[0],2))
            for k in range(sol.shape[0]):
                solCut[k,0]=np.dot(sol[k,:],U1)
                solCut[k,1]=np.dot(sol[k,:],U2)
            # Add field line to 3D plot
            ax1.plot3D(sol[:,0], sol[:,1], sol[:,2], 'red', lw=4)
            # Add field line to 2D plot        
            ax2.plot(solCut[:,0],solCut[:,1], 'red', lw=1)
        else:
            ax1.plot3D(sol[:,0], sol[:,1], sol[:,2], 'gray')

    ax1.view_init(elev=-120., azim=-16)
    ax1.set(xlabel="$X/R_E$ (GSM)")
    ax1.set(ylabel="$Y/R_E$ (GSM)")
    ax1.set(zlabel="$Z/R_E$ (GSM)")
    ax1.axis('square')

    ax1.set_xlim(-4, 4)
    ax1.set_ylim(-4, 4)
    ax1.set_zlim(-4, 4)

if plot3d:
    # Plot plane of field line and x-z plane asociated with y location of event.
    # define the cut plane (as span U1, U2) and the plane parallel to x-z plane
    para1, para2 = np.meshgrid(np.linspace(-1., 3., n), np.linspace(-2., 2., m))
    X_slice = U1[0]*para1+U2[0]*para2+v1[0]*np.ones((n,m))
    Y_slice = U1[1]*para1+U2[1]*para2+v1[1]*np.ones((n,m))
    Z_slice = U1[2]*para1+U2[2]*para2+v1[2]*np.ones((n,m))

    X_o = -1*para1 + 0*para2 + v1[0]*np.ones((n, m))
    Y_o =  0*para1 + 0*para2 + v1[1]*np.ones((n, m))
    Z_o =  0*para1 + 1*para2 + v1[2]*np.ones((n, m))

    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi/2., 100)
    z_day = np.outer(np.cos(u), np.sin(v))
    y_day = np.outer(np.sin(u), np.sin(v))
    x_day = np.outer(np.ones(np.size(u)), np.cos(v))
    v = np.linspace(np.pi/2., np.pi, 100)
    z_night = np.outer(np.cos(u), np.sin(v))
    y_night = np.outer(np.sin(u), np.sin(v))
    x_night = np.outer(np.ones(np.size(u)), np.cos(v))

    #plot the planes in the 3D plot
    from matplotlib import cm as cmap
    if usePatch == False:
        ax1.plot_surface(X_slice, Y_slice, Z_slice, color='g', rstride=1, cstride=1,
                                linewidth=0, antialiased=False)
        ax1.plot_surface(X_o, Y_o, Z_o, color='b', rstride=1, cstride=1,
                                linewidth=0, antialiased=False)
    ax1.plot_surface(x_day, y_day, z_day, color='w', rstride=1, cstride=1,
                           linewidth=0, antialiased=False)
    ax1.plot_surface(x_night, y_night, z_night, color='black', rstride=1, cstride=1,
                           linewidth=0, antialiased=False)
    #ax.set(xlabel='x', ylabel='y', zlabel='z')
    #p = Circle((1, 1), 1)



    if usePatch:
        p = Rectangle((-2, -4), 4, 8, color='b')
        q = Rectangle((-2, -4), 4, 8, color='g')
        ax1.add_patch(p)
        ax1.add_patch(q)
        n=np.array([0,1,0])
        m=U3
        pathpatch_2d_to_3d(p, 0., n)
        pathpatch_2d_to_3d(q, 0., m)
        pathpatch_translate(p, v1)
        pathpatch_translate(q, v1)