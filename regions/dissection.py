import numpy as np

'''
import numpy as np
import dissection as di

Q = np.array([1.5,0.5,0.])
di.dissect(Q, 2., include_flat=False)

Q = np.array([1.5,-0.5,0.])
di.dissect(Q, 2., include_flat=False)

P = np.array([-44,53, 22.])
blocks = di.dissect(P, 64.)
di.plot3dblocks(blocks)

Q = np.array([1.5,0.5,])
di.dissect(Q, 2., include_flat=True)
di.dissect(Q, 2., include_flat=False)
Q = np.array([1.5,0.,])
di.dissect(Q, 2., include_flat=True)
di.dissect(Q, 2., include_flat=False)
Q = np.array([1.5,0.5,0.])
di.dissect(Q, 2., include_flat=True)
di.dissect(Q, 2., include_flat=False)

Q = np.array([1.5,0.,0.])
di.dissect(Q, 2., include_flat=True)
di.dissect(Q, 2., include_flat=False)


Q = np.array([1.5,.5])
ret = di.dissect_2d(Q, 2.)
P = np.array([44,53, 22.])
blocks = di.dissect_3d(P, 64.)
blocks2 = di.dissect(P,64.)
ret2 = di.dissect(Q,2.)
ret2==ret
blocks2==blocks
regs = di.GetRegions(P)

di.plot3dblocks(blocks)


Q = np.array([1.5,0.5])
ret = di.dissect(Q, 2.)
ret
Q = np.array([1.5,0.])
ret = di.dissect(Q, 2.)
ret
Q = np.array([1.5,0.5,0.])
ret = di.dissect(Q, 2.)
Q = np.array([1.5,0.,0.5])
ret = di.dissect(Q, 2.)
ret
'''

def dissect_2d(p, L):
    pm = L/2.
    # want to return an array ret indexed by (a,b,c) where (ret[a,b,0], ret[a,b,1]) is the limits of the bth cartesian component for the ath block, so ret.shape==(nblocks,2,2)

    assert(p.shape==(2,))

    assert(np.all(p)>0)#temporarily

    if np.any(p>=2*pm):
        one_region = np.zeros((2,2))
        one_region[:,0] = -pm
        one_region[:,1] = +pm

        other_region = np.zeros((2,2))
        other_region[:,0] = p - pm
        other_region[:,1] = p + pm

        return np.array([one_region, other_region])

    else:
        sliver = 2*pm - p # == (0+pm) - (p-pm)
        maxdirection = np.argmax(sliver)

        left_block = np.zeros((2,2))
        right_block = np.zeros((2,2))
        center_block = np.zeros((2,2))

        left_block[:,0] = -pm
        left_block[:,1] = +pm
        right_block[:,0] = p - pm
        right_block[:,1] = p + pm

        left_block[maxdirection, 1] = left_block[maxdirection, 1] - sliver[maxdirection]
        right_block[maxdirection, 0] = right_block[maxdirection, 0] + sliver[maxdirection]

        if maxdirection==0:
            otherdirection = 1
        elif maxdirection==1:
            otherdirection = 0

        center_block[maxdirection, :] = np.array([-pm+p[maxdirection], pm])
        center_block[otherdirection, :] = np.array([-pm, pm+p[otherdirection]])

        return np.array([center_block, left_block, right_block])

def dissect_3d(p, L):
    pm = L/2.
    assert(p.shape==(3,))
    # want to return an array ret indexed by (a,b,c) where (ret[a,b,0], ret[a,b,1]) is the limits of the bth cartesian component for the ath block, so ret.shape==(nblocks,3,2)

    assert(np.all(p)>0)#temporarily


    if np.any(p>=2*pm):
        one_region = np.zeros((3,2))
        one_region[:,0] = -pm
        one_region[:,1] = +pm

        other_region = np.zeros((3,2))
        other_region[:,0] = p - pm
        other_region[:,1] = p + pm

        return np.array([one_region, other_region])

    else:
        sliver = 2*pm - p # == (0+pm) - (p-pm)
        maxdirection = np.argmax(sliver)

        left_block = np.zeros((3,2))
        right_block = np.zeros((3,2))

        left_block[:,0] = -pm
        left_block[:,1] = +pm
        right_block[:,0] = p - pm
        right_block[:,1] = p + pm

        left_block[maxdirection, 1] = left_block[maxdirection, 1] - sliver[maxdirection]
        right_block[maxdirection, 0] = right_block[maxdirection, 0] + sliver[maxdirection]

        lower_dim_equivalent_blocks = dissect_2d(np.delete(p, maxdirection), L) # will be (n_lowerblocks, 2, 2)

        n_lowerblocks = lower_dim_equivalent_blocks.shape[0]
        ret = np.zeros((n_lowerblocks+2, 3, 2))

        lower_dim_equivalent_blocks = np.insert(lower_dim_equivalent_blocks, maxdirection, 0, axis=1)
        lower_dim_equivalent_blocks[:, maxdirection, 0] = left_block[maxdirection, 1]
        lower_dim_equivalent_blocks[:, maxdirection, 1] = right_block[maxdirection, 0]

        ret[:n_lowerblocks, :, :] = lower_dim_equivalent_blocks
        ret[-2, :, :] = left_block
        ret[-1, :, :] = right_block

        return ret

def dissect(p, L, include_flat=True):
    needToFlip = p < 0
    p[needToFlip] = -p[needToFlip]

    blocks = dissect_positive(p, L, include_flat=include_flat)

    temp = blocks[:, needToFlip, :].copy() #appears setting a=b[tr] already makes copy, unlike a=b[:], but its set explicitly
    blocks[:, needToFlip, 0] = - temp[:, :, 1]
    blocks[:, needToFlip, 1] = - temp[:, :, 0]

    return blocks


def dissect_positive(p, L, include_flat=True):
    assert(np.all(p)>=0)
    pm = L/2.
    assert(len(p.shape)==1)
    n = p.size # number of dimensions of these hypercubes
    # want to return an array ret indexed by (a,b,c) where (ret[a,b,0], ret[a,b,1]) is the limits of the bth cartesian component for the ath block, so ret.shape==(nblocks,n,2)

    if np.any(p>=2*pm):
        one_region = np.zeros((n,2))
        one_region[:,0] = -pm
        one_region[:,1] = +pm

        other_region = np.zeros((n,2))
        other_region[:,0] = p - pm
        other_region[:,1] = p + pm

        return np.array([one_region, other_region])

    else:
        sliver = 2*pm - p # == (0+pm) - (p-pm)
        maxdirection = np.argmax(sliver)

        left_block = np.zeros((n,2))
        right_block = np.zeros((n,2))

        left_block[:,0] = -pm
        left_block[:,1] = +pm
        right_block[:,0] = p - pm
        right_block[:,1] = p + pm

        left_block[maxdirection, 1] = left_block[maxdirection, 1] - sliver[maxdirection]
        right_block[maxdirection, 0] = right_block[maxdirection, 0] + sliver[maxdirection]

        assert(np.all(left_block[:,0] <= left_block[:,1]))
        #print(left_block[:,0] == left_block[:,1])
        #print(right_block[:,0] == right_block[:,1])
        ignore_left = np.any(left_block[:,0] == left_block[:,1])
        ignore_right = np.any(right_block[:,0] == right_block[:,1])

        if n == 2:
            center_block = np.zeros((2,2))
            if maxdirection==0:
                otherdirection = 1
            elif maxdirection==1:
                otherdirection = 0

            center_block[maxdirection, :] = np.array([-pm+p[maxdirection], pm])
            center_block[otherdirection, :] = np.array([-pm, pm+p[otherdirection]])

            if not include_flat and ignore_left:
                left_block = np.nan*np.empty((2,2))
            if not include_flat and ignore_right:
                right_block = np.nan*np.empty((2,2))

            return np.array([center_block, left_block, right_block])

        elif n > 2:
            lower_dim_equivalent_blocks = dissect(np.delete(p, maxdirection), L, include_flat=include_flat) # will be (n_lowerblocks, n-1, 2)

            n_lowerblocks = lower_dim_equivalent_blocks.shape[0]
            ret = np.nan*np.empty((n_lowerblocks+2, n, 2))

            lower_dim_equivalent_blocks = np.insert(lower_dim_equivalent_blocks, maxdirection, 0, axis=1)
            lower_dim_equivalent_blocks[:, maxdirection, 0] = left_block[maxdirection, 1]
            lower_dim_equivalent_blocks[:, maxdirection, 1] = right_block[maxdirection, 0]

            ret[:n_lowerblocks, :, :] = lower_dim_equivalent_blocks


            if include_flat or not ignore_left:
                ret[-2, :, :] = left_block
            if include_flat or not ignore_right:
                ret[-1, :, :] = right_block

            return ret

        else:
            assert(False)


def GetRegions(point):
    globalXmin = -224.
    globalXmax = 32.
    globalYmin = -128.
    globalYmax = 128.
    globalZmin = -128.
    globalZmax = 128.

    assert(point.shape==(3,))
    pm = 31.875
    d = 0.25
    p = d*np.round(point/d)

    blocks = dissect(p, 2.*(pm + d/2.), include_flat=False)
    ret = []
    for i in range(blocks.shape[0]):
        reg = {'d' : d}
        if np.any(np.isnan(blocks[i,:,:])):
            continue
        reg['xlims'] = tuple( blocks[i,0,:] + np.array([d/2.,-d/2.]) )
        reg['ylims'] = tuple( blocks[i,1,:] + np.array([d/2.,-d/2.]) )
        reg['zlims'] = tuple( blocks[i,2,:] + np.array([d/2.,-d/2.]) )

        if reg['xlims'][1] > globalXmax:
            reg['xlims'] = ( reg['xlims'][0], d*np.floor(globalXmax/d)-d/2.)
        if reg['ylims'][1] > globalYmax:
            reg['ylims'] = ( reg['ylims'][0], d*np.floor(globalYmax/d)-d/2.)
        if reg['zlims'][1] > globalZmax:
            reg['zlims'] = ( reg['zlims'][0], d*np.floor(globalZmax/d)-d/2.)

        if reg['xlims'][0] < globalXmin:
            reg['xlims'] = ( d*np.floor(globalXmin/d)+d/2., reg['xlims'][1] )
        if reg['ylims'][0] < globalYmin:
            reg['ylims'] = ( d*np.floor(globalYmin/d)+d/2., reg['ylims'][1] )
        if reg['zlims'][0] < globalZmin:
            reg['zlims'] = ( d*np.floor(globalZmin/d)+d/2., reg['zlims'][1] )

        ret.append(reg.copy())

    return tuple(ret)



def plot3dblocks(blocks):
    #adapted from https://codereview.stackexchange.com/questions/155585/plotting-a-rectangular-prism (although didnt work)
    #https://matplotlib.org/3.1.1/gallery/mplot3d/wire3d.html
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    #ax.set_aspect("equal")

    # draw cube
    def rect_prism(x_range, y_range, z_range):
        # TODO: refactor this to use an iterator
        xx, yy = np.meshgrid(x_range, y_range)
        ax.plot_wireframe(xx, yy, z_range[0]*np.ones(xx.shape), color="r")
        ax.plot_surface(xx, yy, z_range[0]*np.ones(xx.shape), color="r", alpha=0.2)
        ax.plot_wireframe(xx, yy, z_range[1]*np.ones(xx.shape), color="r")
        ax.plot_surface(xx, yy, z_range[1]*np.ones(xx.shape), color="r", alpha=0.2)


        yy, zz = np.meshgrid(y_range, z_range)
        ax.plot_wireframe(x_range[0]*np.ones(yy.shape), yy, zz, color="r")
        ax.plot_surface(x_range[0]*np.ones(yy.shape), yy, zz, color="r", alpha=0.2)
        ax.plot_wireframe(x_range[1]*np.ones(yy.shape), yy, zz, color="r")
        ax.plot_surface(x_range[1]*np.ones(yy.shape), yy, zz, color="r", alpha=0.2)

        xx, zz = np.meshgrid(x_range, z_range)
        ax.plot_wireframe(xx, y_range[0]*np.ones(zz.shape), zz, color="r")
        ax.plot_surface(xx, y_range[0]*np.ones(zz.shape), zz, color="r", alpha=0.2)
        ax.plot_wireframe(xx, y_range[1]*np.ones(zz.shape), zz, color="r")
        ax.plot_surface(xx, y_range[1]*np.ones(zz.shape), zz, color="r", alpha=0.2)

    for i in range(blocks.shape[0]):
        rect_prism(blocks[i,0,:], blocks[i,1,:], blocks[i,2,:])
        #rect_prism(np.array([-1, 1]), np.array([-1, 1]), np.array([-0.5, 0.5]))
        #rect_prism(np.array([2, 3]), np.array([2, 3]), np.array([-2, -1]))

    plt.show()
    plt.clf()
