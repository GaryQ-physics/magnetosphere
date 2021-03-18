import os
import sys
import numpy as np
import pandas as pd
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf
import util
from units_and_constants import phys

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../swmf/')
import read_swmf_files as rswmf

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../../jitFORTRAN/')
import jitFORTRAN

import datetime
from native_grid2 import main

run = 'DIPTSUR2'
n_times = 2
para = False

times = list(util.get_available_slices(run)[1])
if n_times is not None:
    times = times[:n_times]


########################################################################
########################################################################
########################################################################
frobenius_script = '''!FROBENIUS_NORM
subroutine FROBENIUS_NORM(PARTIALS_T, FNORM, nBlock,nI,nJ,nK)
    implicit none

    integer, intent(in) :: nBlock,nI,nJ,nK
!f2py intent(in) :: nBlock,nI,nJ,nK
    real(4), intent(in), dimension(3,3,nK,nJ,nI,nBlock) :: PARTIALS_T
!f2py intent(in) :: PARTIALS_T
    real(4), intent(inout), dimension(nK,nJ,nI,nBlock) :: FNORM
!f2py intent(in,out) :: FNORM

    integer :: iBlock,i,j,k, a,b


    do iBlock=1,nBlock
    do i = 1,nI
    do j = 1,nJ
    do k = 1,nK
        FNORM(k,j,i,iBlock) = 0.
        do a = 1,3
        do b = 1,3
            FNORM(k,j,i,iBlock)=FNORM(k,j,i,iBlock) + PARTIALS_T(a,b, k,j,i,iBlock)**2
        end do
        end do
        FNORM(k,j,i,iBlock) = SQRT(FNORM(k,j,i,iBlock))
    end do
    end do
    end do
    end do

end subroutine FROBENIUS_NORM
'''

FROBENIUS_NORM_F = jitFORTRAN.Fortran_Subroutine(frobenius_script, 'FROBENIUS_NORM')
FROBENIUS_NORM_F.compile()
def frobenius_norm(partials):
    nBlock,nI,nJ,nK, A,B = partials.shape
    assert(A==B==3)

    fnorm = np.empty((nK,nJ,nI,nBlock), dtype=np.float32, order='F')
    fnorm[:,:,:,:] = np.nan
    FROBENIUS_NORM_F.execute(partials.transpose(), fnorm, nBlock,nI,nJ,nK)
    return fnorm.transpose()
########################################################################
########################################################################
########################################################################

if para:
    from joblib import Parallel, delayed
    import multiprocessing
    num_cores = multiprocessing.cpu_count()
    if num_cores is not None and num_cores > len(times):
        num_cores = len(times)
    print('Parallel processing {0:d} time slices using {1:d} cores'\
          .format(len(times), num_cores))
    Parallel(n_jobs=num_cores)(\
            delayed(main)(run, time, frobenius_norm=frobenius_norm) for time in times)
else:
    for time in times:
        main(run, time, frobenius_norm=frobenius_norm)


epsilons = [ 0.0625,
             0.1250,
             0.2500,
             0.5000,
             1.0000,
             2.0000,
             4.0000,
             8.0000 ]


direct = conf[run+'_derived'] + 'derivatives/native_grid/'
list_summaries = []
dtimes = []
for time in times:
    dtimes.append(datetime.datetime(time[0],time[1],time[2],time[3],time[4],time[5]))
    fname_summary = direct + '%.2d%.2d%.2dT%.2d%.2d%.2d_summary.pkl'%util.tpad(time, length=6)
    list_summaries.append(pd.read_pickle(fname_summary))

columns = list_summaries[0].columns

for epsilon in epsilons:
    SummaryTimeseries = pd.DataFrame(index=dtimes, columns=columns)

    for col, series in SummaryTimeseries.items():
        for i in range(len(times)):
            SummaryTimeseries[col][dtimes[i]] = list_summaries[i][col][epsilon]

    SummaryTimeseries.to_pickle(direct+'SummaryTimeseries_epsilon_%f.pkl'%(epsilon))


