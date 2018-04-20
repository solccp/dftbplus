  !* integrate modes
  real(dp), allocatable :: eigenValues(:)
  integer  :: iCount, jCount, ii, jj, kk
  !projection variables
  real(dp), allocatable :: basis(:,:), newhess(:,:), dotpr(:,:)
!  real(dp) :: dotpr(6,6)
  real(dp), parameter :: offsetlinear=1.0e-4_dp
  real(dp) :: det, det1, det2, rnorm
  logical :: centinv, tNormalModes=.true.
  integer :: nindep, ngerade, iseed=1234235
  integer :: lwork, n3, info, licz, ll
  real(dp), dimension(:), allocatable :: work
  real(dp) :: tmp(1), diffinv, distinv, xvec(3), rr, prod, xvalplus, xvalminus
  integer, dimension(:), allocatable  :: invpartner
  real(dp), dimension(:, :), allocatable :: coordcm
  integer :: fdNormalModes
  character(*), parameter :: normalModesOut = "vibrations.molden"
  real(dp), dimension(:), allocatable :: redmas
  real(dp) :: massx
  type(ORanlux), allocatable :: randomFreq
  type(ORandomGenPool) :: randGenPool

