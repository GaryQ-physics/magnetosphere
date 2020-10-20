!f2py -c deltaB_subroutine.f90 -m deltaB_subroutine

subroutine DELTAB(A, B, N, M)

    implicit none
    !will pass the transposed np array (which is in fortran order) having shape (3, N, M)

    integer, intent(in) :: N, M
    real(8), intent(in), dimension(0:2, 0:N-1, 0:M-1) :: B
!f2py intent(in) :: B
    real(8), intent(inout), dimension(0:2, 0:N-1, 0:M-1) :: A
!f2py intent(in,out) :: A

    integer :: j, k
    real(8), dimension(0:2) :: CACHE


    do k=0,M-1
        do j=0,N-1
            CACHE(2) = A(0,j,k)*B(1,j,k) - A(1,j,k)*B(0,j,k)
            CACHE(0) = A(1,j,k)*B(2,j,k) - A(2,j,k)*B(1,j,k)
            CACHE(1) = A(2,j,k)*B(0,j,k) - A(0,j,k)*B(2,j,k)
            A(:,j,k) = CACHE
        end do
    end do !del cache  ?

end subroutine DELTAB


subroutine ROUT(A, B, N, M)
    implicit none
    integer, intent(in) :: n
    integer, intent(out) :: plus1, plus2 ! Output
    plus1 = n + 4
    plus2 = n + 5
end subroutine ROUT


subroutine ARRTEST(X, Y, n, m)
    implicit none
    ! with X.shape==Y.shape==(n,m), its eq to:
    ! X *= Y    numpy arrays (order='F', otherwise like X=X*Y (copies) )
    ! note n,m are supplied automatically if left blank

    integer, intent(in) :: n, m
    real(8), intent(in), dimension(0:n-1, 0:m-1) :: Y
!f2py intent(in) :: Y
    real(8), intent(inout), dimension(0:n-1, 0:m-1) :: X
!f2py intent(in,out) :: X

    integer :: i, j


    do j=0,m-1 ! fastest with columns index as outer loop (cause fortran uses column ordering, as oposed to C (numpys defauls) row ordering)
        do i=0,n-1 ! carefull, seemd to cause strange segmentation faults at end when 0,n
            X(i,j) = X(i,j)*Y(i,j)
        end do
    end do

end subroutine ARRTEST


subroutine CROSS(A, B, N, M)
    implicit none
    !will pass the transposed np array (which is in fortran order) having shape (3, N, M)

    integer, intent(in) :: N, M
    real(8), intent(in), dimension(0:2, 0:N-1, 0:M-1) :: B
!f2py intent(in) :: B
    real(8), intent(inout), dimension(0:2, 0:N-1, 0:M-1) :: A
!f2py intent(in,out) :: A

    integer :: j, k
    real(8), dimension(0:2) :: CACHE


    do k=0,M-1
        do j=0,N-1
            CACHE(2) = A(0,j,k)*B(1,j,k) - A(1,j,k)*B(0,j,k)
            CACHE(0) = A(1,j,k)*B(2,j,k) - A(2,j,k)*B(1,j,k)
            CACHE(1) = A(2,j,k)*B(0,j,k) - A(0,j,k)*B(2,j,k)
            A(:,j,k) = CACHE
        end do
    end do

    !del cache  ?

end subroutine CROSS

