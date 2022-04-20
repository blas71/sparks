subroutine get_OLCC(dt, nLCC, V, KLCC, Ca, states_LCC, OLCC)
implicit none
! Input 1
integer, intent(in) :: nLCC
double precision, intent(in) :: dt, V, KLCC, Ca
double precision, intent(out) :: OLCC
! In - out
integer, intent(inout) :: states_LCC(nLCC)

! INTERNAL VARIABLES
integer :: state
double precision :: A(5,5)
double precision, dimension(:), allocatable :: probs, cumprobs
! Dummy variables
integer :: i, ii
double precision :: prob1(1)

allocate(probs(5), cumprobs(5))

! Get transition rates
call transitions_rates_LCC(V, KLCC, Ca, A)
A = A * dt

! Open states fraction
OLCC = 0d0

! Dynamics
do i = 1, nLCC

!  Condition for LCC(i)
   state = states_LCC(i)

!  PROBABILITIES
   probs = A(state,:)

!  Cumulative probabilities
   call cumsum(5, probs, cumprobs)

!  Probabilities
   !call ranmar(1, prob1)
   call random_number(prob1)
   ii = 1
   do while (ii .le. 5)
      if (prob1(1) .lt. cumprobs(ii)) then
         state = ii
         exit
      end if
      ii = ii + 1
   end do

!  Update state vector
   states_LCC(i) = state

!  Get OLCC
   if (state .eq. 3) then
      OLCC = OLCC + 1d0 / nLCC
   end if
end do

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine transitions_rates_LCC(V, KLCC, cd, A)
implicit none
! Input 1
double precision, intent(in) :: V, KLCC, cd
! Ouput
double precision, intent(out) :: A(5,5)
! Internal variables
double precision :: pinf, p0, rho0, p0f, tau0, pif, fca

pinf = 1d0 / (1d0 + exp(- (V - 15d0) / 8d0))
p0 = 1d0 / (1d0 + exp(- (V + 40d0) / 10d0))
rho0 = 10d0 + 495d0 * exp(V / 15.6d0)
p0f = 1d0 - 1d0 / (1d0 + exp(- (V + 40d0) / 4d0))
tau0 = (rho0 - 450d0) * p0f + 450d0
pif = 1d0 / (1d0 + exp(- (V + 40d0) / 3d0))
fca = 1d0 / (1d0 + (KLCC / cd)**3d0)

A = 0d0

A(2,3) = 0.5d0
A(3,2) = 2d0
A(4,2) = 2.24e-3

A(1,2) = pinf
A(2,1) = 1d0 - pinf
A(1,5) = p0 / tau0
A(5,1) = (1d0 - p0) / tau0
A(4,5) = (1d0 - pif) / 3d0
A(2,4) = 0.00413d0 + 0.024d0 * fca
A(3,4) = 0.00195d0 + 0.01826d0 * fca
A(2,3) = 0.5d0
A(3,2) = 2d0
A(4,2) = 2.24e-3
A(4,3) = A(3,4) * A(2,3) * A(4,2) / (A(3,2) * A(2,4))
A(5,4) = A(4,5) * A(5,1) * A(2,4) * A(1,2) / (A(1,5) * A(4,2) * A(2,1))

end subroutine
