!****************************************************************
!*               subroutine CONGRUENTIAL                        *
!****************************************************************

subroutine congruential(a, b, N, seed, z)
implicit none 

! Function arguments
integer N
integer*8 a, b, seed
double precision z(N) 

! Local variables
integer*8 zaux(N)
integer i
      
zaux(1) = mod(a * seed, b)

do i = 1, N-1
   zaux(i+1) = mod(a * zaux(i), b)
end do

z = zaux / (dfloat(b) - 1.d0)

return
end 

!****************************************************************
!*               subroutine INIT_RANMAR                         *
!****************************************************************

subroutine init_ranmar()
implicit none

! Common block containing state of the prng
integer jmin, jmax            ! j is pointing to last position in x
double precision x(97), y
common /ranmar_state/ jmin, jmax, x, y

! Variables to call congruential prng
integer*8 a, b
parameter(a = 16871, b = 2147483647)
integer*8 seed
parameter(seed = 2324456)
integer N
parameter (N = 2000)
double precision x_i(N)

! Other local variables
integer i

call congruential(a, b, N, seed, x_i)
      
do i = 1, 97
   x(i) = x_i(i + N - 97)
end do

y = 6701772.d0 / (2.d0**31.d0 - 1.d0)

jmin = 33
jmax = 97
      
return 
end

!****************************************************************
!*                 subroutine RANMAR                            *
!****************************************************************

subroutine ranmar(N, z)
implicit none

! Subroutine arguments
integer N
double precision z(N)

! Common block containing state of the prng
integer jmin, jmax            ! j is pointing to last position in x
double precision x(97), y
common /ranmar_state/ jmin, jmax, x, y

! Parameters of RANMAR
double precision c, d
parameter (c = 7654321.d0 / 16777216.d0)
parameter (d = 16777213.d0 / 16777216.d0)

! Local variables
integer i

do i = 1, N

!  First random number
   y = y - c
   if (y .lt. 0.d0) then          
      y = y + d
   end if 

!  Second random number
   x(jmax) = x(jmax) - x(jmin)
   if (x(jmax) .lt. 0.d0) then
      x(jmax) = x(jmax) + 1.d0
   end if

!  Desired random number
   z(i) = x(jmax) - y
   if (z(i) .lt. 0.d0) then 
      z(i) = z(i) + 1.d0
   end if

!  Increase j
   jmin = jmin + 1
   if (jmin .gt. 97) then
      jmin = 1
   end if

   jmax = jmax + 1
   if (jmax .gt. 97) then
      jmax = 1
   end if

end do
      
return 
end

!****************************************************
!            subroutine NORMAL_DIST                 *
!****************************************************

subroutine normal_dist(N, sigma, z)
implicit none

! Subroutine arguments
integer N
double precision sigma, z(N)
      
!     Local variables
double precision  pi, u1(N/2 + 1), u2(N/2 + 1)
integer i, imax
parameter (pi = 3.14159265358979323846d0)
      
call ranmar(N/2 + 1, u1)
call ranmar(N/2 + 1, u2)

imax = mod(N, 2) + N/2

do i = 1, N/2
   z(2*i) = sigma * sqrt(-2 * log(u1(i))) * sin(2 * pi * u2(i))
end do

do i = 1, imax
   z(2*i-1) = sigma * sqrt(-2 * log(u1(i))) * cos(2 * pi * u2(i))
end do

return
end
