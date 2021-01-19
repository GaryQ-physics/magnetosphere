


def GetDel(run, time, field, points, para=False, epsilon=0.0625, debug=False):
    filename = util.time2CDFfilename(run, time)

    # field = 'b_biotsavart' ; 'b_batsrus' ; 'b1_batsrus' ; 'j_batsrus' ; 'testinterp'

    def func(i):
        point = points[i,:].copy()
        xplus = point + epsilon*np.array([1,0,0])
        xmin  = point - epsilon*np.array([1,0,0])
        yplus = point + epsilon*np.array([0,1,0])
        ymin  = point - epsilon*np.array([0,1,0])
        zplus = point + epsilon*np.array([0,0,1])
        zmin  = point - epsilon*np.array([0,0,1])
        pts = np.array([xplus, xmin, yplus, ymin, zplus, zmin])
        if debug: print(pts.shape)

        if field == 'b_biotsavart':
            #pm = 31.875
            #reg =  {'xlims': (-pm, pm),
            #        'ylims': (-pm, pm),
            #        'zlims': (-pm, pm),
            #        'd': 0.25
            #        }
            regs = di.GetRegions(point)
            Fs = bs.biot_savart_run(run, time, pts, regs, summed=True, separateRegions=False)
        elif field == 'testinterp':
            Fs = np.column_stack([Bx_testinterp(pts), By_testinterp(pts), Bz_testinterp(pts)])
        elif '_batsrus' in field:
            fi = field.split('_batsrus')[0]
            if debug: print('fi=%s'%(fi))
            Fs = probe(filename, pts, var=[fi+'x',fi+'y',fi+'z'], library='kameleon')
        else:
            raise ValueError('invalid field string')

        #print(Fs.shape)
        delF = np.nan*np.empty((3,3))
        delF[0,:] = ( Fs[0,:] - Fs[1,:] )/(2.*epsilon)
        delF[1,:] = ( Fs[2,:] - Fs[3,:] )/(2.*epsilon)
        delF[2,:] = ( Fs[4,:] - Fs[5,:] )/(2.*epsilon)
        #print(delF)

        return delF

    if para:
        assert(field != 'b_biotsavart') #would run out of memory anyway

        #paralelize

    else:
        ret = []
        for i in range(points.shape[0]):
            ret.append(func(i))

    return np.array(ret)

def GetDivergence(delF):
    divF = np.nan*np.empty(delF.shape[0])
    for i in range(delF.shape[0]):
        divF[i] = delF[i,0,0] + delF[i,1,1] + delF[i,2,2]

    return divF

def GetCurl(delF):
    curlF = np.nan*np.empty((delF.shape[0], 3))
    for i in range(delF.shape[0]):
        curlF_tens = delF[i,:,:] - delF[i,:,:].transpose()
        curlF[i,:] = np.array([ curlF_tens[1,2], curlF_tens[2,0], curlF_tens[0,1] ])

    return curlF

def GetFrobeniusNormDel(delF):
    ret = np.nan*np.empty(delF.shape[0])
    for i in range(delF.shape[0]):
        ret[i] = np.linalg.norm(delF[i,:,:], ord='fro') # frobenius norm, equivalently 2-norm of flattened vector. It's equivalent to set ord=None
    return ret

def GetOperatorNormDel(delF):
    ret = np.nan*np.empty(delF.shape[0])
    for i in range(delF.shape[0]):
        ret[i] = np.linalg.norm(delF[i,:,:], ord=2) # operator norm inherited from column vector 2-norm, equivalently largest singular value.
    return ret



