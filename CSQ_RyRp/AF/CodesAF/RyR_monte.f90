subroutine get_ORyR(dt, nRyR, koc, ki1i2, kic, kio, k2, k7, k6, k11, k1d, k3d, k8d, &
         Ca, Casr, q, alpha, states_RyR, ORyR)
implicit none
! Input
integer :: nRyR
double precision, intent(in) :: dt, Ca, Casr, q, alpha
double precision, intent(in) :: koc, ki1i2, kic, kio, k2, k7, k6, k11, k1d, k3d, k8d
! In - out
integer, intent(inout) :: states_RyR(nRyR)
! Output
double precision, intent(out) :: ORyR
! Internal variables
double precision, dimension(:,:), allocatable :: K
double precision, dimension(:), allocatable :: probs, cumprobs
integer :: neighbour(4,2)
! Dummy variables
integer :: ii, state, j
double precision :: prob1(1)

allocate(K(4,2), probs(2), cumprobs(2))

! Near neighbours
neighbour(1,:) = (/ 2, 3/) ! Close state connections 
neighbour(2,:) = (/ 1, 4/) ! Open state connections 
neighbour(3,:) = (/ 1, 4/) ! I2 state connections 
neighbour(4,:) = (/ 2, 3/) ! I1 state connections 


! Transition matrix with probabilities
call matrix(dt, Ca, Casr, q, alpha, &
         koc, ki1i2, kic, kio, k2, k7, k6, k11, k1d, k3d, k8d, K)

! Open states fraction
ORyR = 0d0

! Dynamics
do j = 1, nRyR
   ! Generate the transition for each RyR
   state = states_RyR(j)
   probs = K(state,:)
   
   ! Cumulative probabilities
   call cumsum(2, probs, cumprobs)
   
   ! Probabilities
   call random_number(prob1)
   ii = 1
   do while (ii .le. 2)
      if (prob1(1) .lt. cumprobs(ii)) then
         state = neighbour(state, ii)
         exit
      end if
      ii = ii + 1
   end do

!  Update state vector
   states_RyR(j) = state

!  Get ORyR
   if (state .eq. 2) then
      ORyR = ORyR + 1d0
   end if
end do

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine matrix(dt, Ca, Casr, q, alpha, &
         koc, ki1i2, kic, kio, k2, k7, k6, k11, k1d, k3d, k8d, K)
implicit none
! Inputs
double precision, intent(in) :: dt, Ca, Casr, q, alpha, koc, ki1i2, kic, kio
double precision, intent(in) :: k2, k7, k6, k11, k1d, k3d, k8d
! Output
double precision, intent(out) :: K(4,2)
! internal variables
double precision :: kCasr, kco, ki2i1, koi, kci, g, f2

g = 3d0

! Matrix with probabilities of change the state
K = 0.

! Get transition rates 
kco = alpha * 3d0 * ((k2*k1d)**g + (k11*q)**g) / (k1d**g + q**g) * Ca**3 / (1 + (Ca/20.)**3)
ki2i1 = alpha * 3d0 * ((k2*k8d)**g + (k11*q)**g) / (k8d**g + q**g) * Ca**3 / (1 + (Ca/20.)**3)
koi = 2e-4 * Ca
kci = 2e-4 * Ca

kCasr = (1 + (700d0/Casr)**30d0) / (1d0 + (Casr/1000d0)**5d0)
f2 = 1d0 + (500d0/Casr)**2d0

! Close to ...
K(1,1) = kco / kCasr
K(1,2) = kci * f2
! Open to ...
K(2,1) = koc
K(2,2) = koi * f2
! I2 to ...
K(3,1) = kic
K(3,2) = ki2i1 / kCasr
! I1 to ...
K(4,1) = kio
K(4,2) = ki1i2

K = K * dt

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine cumsum(N, x, y)
implicit none
! Inputs
integer, intent(in) :: N
double precision, dimension(N), intent(in) :: x
! Output
double precision, dimension(N), intent(out) :: y
! Internal variables
integer i

y = 0.

do i = 1, N
   y(i) = sum(x(1:i))
end do

end subroutine
