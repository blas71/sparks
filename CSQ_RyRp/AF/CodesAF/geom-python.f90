program cell2D
implicit none
! GEOMETRY
integer :: Mx, My
integer :: nRyR_int, nRyR_sarco, nRyRtot, nLCC_sarco, nNCX, n_SR, n_TnC, nTT, nAT, nLCCtot, n_Back
integer, dimension(:,:), allocatable :: ind_RR_int, ind_RR_sarco, ind_LCC, ind_NCX, ind_TT, ind_AT, ind_Back
integer, dimension(:), allocatable :: nRyRPixel
integer, dimension(:,:), allocatable :: ind_RR, ind_SR, ind_TnC, ind_LCC_sarco
double precision :: Lx, Ly, dx
! CALCIUM
double precision, dimension(:,:), allocatable :: Ca, CaSr, CaSrTot, casave, cabsave
double precision, dimension(:,:), allocatable :: CabCM, CaTnC, CabSR, CabRhod
double precision, dimension(:,:), allocatable :: lapCa, lapCaSr
! currents
double precision, dimension(:,:), allocatable :: J_currents_Ca, J_currents_CaSr
double precision, dimension(:,:), allocatable :: J_currents_CabCM, J_currents_CabRhod
double precision, dimension(:,:), allocatable :: J_currents_CabTnC, J_currents_CabSR
double precision :: JNaCa, Jup, Jrel, JCaL
double precision :: JNaCaT, JCaLT
double precision :: JupT, JrelT
double precision :: JbuffCT, JbuffC
double precision :: Temp, F, R 
double precision :: gNaCa, Ca0, Na0, Nai, Kda, ksat, nu, KmNa0, KmNai, KmCai, KmCa0
double precision :: gup, Ki, Ksr
double precision :: grel, koc, ki1i2, kic, kio
double precision :: gCaL, KLCC
double precision :: k7, k6, k11, k1d, k3d, k8d, k2
! Buffers
double precision :: konCM, koffCM, BCAM, BSR
double precision :: kon, koff, konSR, koffSR, BT_val
double precision :: Ksrtot, Bsrtot0, Bsrtot1
double precision :: konRhod, koffRhod, BTRhod
double precision, dimension(:,:), allocatable :: bcoeff, ccoeff, BsrTot, BT
! RyR and LCC
double precision :: ORyR, OLCC, P0, P1
integer, dimension(:,:), allocatable :: states_RyR, states_LCC
double precision, dimension(:,:), allocatable :: alpha
integer, dimension(:,:,:), allocatable :: tmatrix
integer :: nRyR, nRyR_membrane, nLCC, openRyRsi(4)
! Diffusion coefficient
double precision, dimension(:,:), allocatable :: Di, DSR, Di_full, Dsr_full
double precision :: Dcoef_i, Dcoef_sr
double precision :: vfrac_SR, vfrac_ci
double precision, dimension(:,:), allocatable :: Vi, Vsr, ViVsr
! Potential parameters
double precision :: V, Vmax, Vrest
! Time variables
integer :: Nt, Mt
double precision :: dt, t, tt, Ts, ADP
! save variables
integer :: nplot, nsave, nfilm, nopen, nexp, nplot2, isave
! Dummy variables
integer :: i, j, k, l
!Recording files for VTK
character*40 Calc, Cab_name
integer :: umiler, nplot_cent, cent, nplot_des, des, nplot_uni

! read parameters
open(unit=1, file='../data/parameters.data', status='old')
read(1,*) Lx
read(1,*) Ly
read(1,*) dx
read(1,*) nRyR_int
read(1,*) nRyR_sarco
read(1,*) nLCC_sarco
read(1,*) nTT
read(1,*) nAT
read(1,*) nNCX
read(1,*) n_SR
read(1,*) n_TnC
read(1,*) n_Back
close(1)

nRyRtot = nRyR_int + nRyR_sarco
nLCCtot = nLCC_sarco + nTT + nAT

include 'params.f90'

! spatial dimenions
Mx = int(Lx / dx)
My = int(Ly / dx)

! set dimensions
allocate(ind_RR_int(nRyR_int,2), ind_RR_sarco(nRyR_sarco,2), ind_LCC(nLCCtot,2))
allocate(nRyRPixel(nRyRtot))
allocate(ind_RR(nRyRtot,2), ind_NCX(nNCX,2), ind_LCC_sarco(nLCC_sarco,2), ind_TT(nTT,2), ind_AT(nAT,2))
allocate(ind_SR(n_SR,2), ind_TnC(n_TnC,2), ind_Back(n_Back,2))
allocate(Ca(Mx,My), CaSr(Mx,My), lapCa(Mx,My), lapCaSr(Mx,My))
allocate(CabCM(Mx,My), CaTnC(Mx,My), CabSR(Mx,My), CaSrTot(Mx,My))
allocate(CabRhod(Mx,My), casave(Mx,My), cabsave(Mx,My), alpha(Mx,My))
allocate(bcoeff(Mx,My), ccoeff(Mx,My), BsrTot(Mx,My), BT(Mx,My))
allocate(J_currents_Ca(Mx,My), J_currents_CaSr(Mx,My), J_currents_CabCM(Mx,My))
allocate(J_currents_CabTnC(Mx,My), J_currents_CabSR(Mx,My))
allocate(J_currents_CabRhod(Mx,My))
allocate(states_RyR(nRyRtot,nRyR), states_LCC(nLCCtot,nLCC))
allocate(tmatrix(nRyRtot,nRyR,2))
allocate(Vi(Mx,My), Vsr(Mx,My), ViVsr(Mx,My))
allocate(Di(Mx,My), DSR(Mx,My))
allocate(Di_full(0:Mx+1,0:My+1), Dsr_full(0:Mx+1,0:My+1))

! read indices
open(unit=1, file='../data/posRyR.data', status='old')
do i = 1, nRyR_int
   read(1,*) ind_RR_int(i,:)
end do
close(1)
do i = 1, nRyR_int
   ind_RR_int(i,:) = ind_RR_int(i,:) + 1
   ind_RR(i,:) = ind_RR_int(i,:)
end do

open(unit=1, file='../data/posRyR-sarco.data', status='old')
do i = 1, nRyR_sarco
   read(1,*) ind_RR_sarco(i,:)
end do
close(1)
do i = 1, nRyR_sarco
   ind_RR_sarco(i,:) = ind_RR_sarco(i,:) + 1
   ind_RR(i+nRyR_int,:) = ind_RR_sarco(i,:)
end do

open(unit=1, file='../data/nRyRpixel.data', status='old')
do i = 1, nRyR_int
   read(1,*) nRyRPixel(i)
end do
close(1)
do i = 1, nRyR_sarco
   nRyRPixel(i+nRyR_int) = nRyR_membrane
end do

open(unit=1, file='../data/posLCC.data', status='old')
do i = 1, nLCC_sarco
   read(1,*) ind_LCC_sarco(i,:) 
end do
close(1)
do i = 1, nLCC_sarco
   ind_LCC_sarco(i,:) = ind_LCC_sarco(i,:) + 1
   ind_LCC(i,:) = ind_LCC_sarco(i,:)
end do

open(unit=1, file='../data/posTT-LCC.data', status='old')
do i = 1, nTT
   read(1,*) ind_TT(i,:) 
end do
close(1)
do i = 1, nTT
   ind_TT(i,:) = ind_TT(i,:) + 1
   ind_LCC(i+nLCC_sarco,:) = ind_TT(i,:)
end do

open(unit=1, file='../data/posAT-LCC.data', status='old')
do i = 1, nAT
   read(1,*) ind_AT(i,:) 
end do
close(1)
do i = 1, nAT
   ind_AT(i,:) = ind_AT(i,:) + 1
   ind_LCC(i+nLCC_sarco+nTT,:) = ind_AT(i,:)
end do

open(unit=1, file='../data/posNCX.data', status='old')
do i = 1, nNCX
   read(1,*) ind_NCX(i,:) 
end do
close(1)
do i = 1, nNCX
   ind_NCX(i,:) = ind_NCX(i,:) + 1
end do

BT = BT_val

open(unit=1, file='../data/posSR.data', status='old')
do i = 1, n_SR
   read(1,*) ind_SR(i,:) 
end do
close(1)
do i = 1, n_SR
   ind_SR(i,:) = ind_SR(i,:) + 1
   BT(ind_SR(i,1), ind_SR(i,2)) = BT_val / 5d0
end do

open(unit=1, file='../data/posTnC.data', status='old')
do i = 1, n_TnC
   read(1,*) ind_TnC(i,:) 
end do
close(1)
do i = 1, n_TnC
   ind_TnC(i,:) = ind_TnC(i,:) + 1
end do

open(unit=1, file='../data/posBackground.data', status='old')
do i = 1, n_Back
   read(1,*) ind_Back(i,:) 
end do
close(1)
do i = 1, n_Back
   ind_Back(i,:) = ind_Back(i,:) + 1
end do

! --- VOLUME FACTORS and DIFFUSION ---
open(unit=70, file='../data/diff-coef.data', status='unknown')
! cytosol zone: it's the same as TnC zone
call diffusion_cytosol(vfrac_ci, Dcoef_i)
call diffusion_SR(vfrac_ci, Dcoef_sr)
do i = 1, n_TnC
   Vi(ind_TnC(i,1), ind_TnC(i,2)) = 1d0 - vfrac_ci
   Vsr(ind_TnC(i,1), ind_TnC(i,2)) = vfrac_ci
   Di(ind_TnC(i,1), ind_TnC(i,2)) = 0.2d0 * Vi(ind_TnC(i,1), ind_TnC(i,2))
   Dsr(ind_TnC(i,1), ind_TnC(i,2)) = 0.09d0 * Vsr(ind_TnC(i,1), ind_TnC(i,2))
end do
write(70,*) 'cytosol zone: Di = ', Dcoef_i
write(70,*) 'cytosol zone: Dsr = ', Dcoef_sr

! SR zone
call diffusion_cytosol(vfrac_SR, Dcoef_i)
call diffusion_SR(vfrac_SR, Dcoef_sr)
do i = 1, n_SR
   Vi(ind_SR(i,1), ind_SR(i,2)) = 1d0 - vfrac_SR
   Vsr(ind_SR(i,1), ind_SR(i,2)) = vfrac_SR
   Di(ind_SR(i,1), ind_SR(i,2)) = 0.18d0 * Vi(ind_SR(i,1), ind_SR(i,2))
   Dsr(ind_SR(i,1), ind_SR(i,2)) = 0.1d0 * Vsr(ind_SR(i,1), ind_SR(i,2))
end do
write(70,*) 'SR zone: Di = ', Dcoef_i
write(70,*) 'SR zone: Dsr = ', Dcoef_sr
close(70)

! Volume ratio
ViVsr = Vi / Vsr

! extend matrix with the diffusion BC 
call diff_bc(Mx, My, Di, Di_full)
call diff_bc(Mx, My, Dsr, Dsr_full)

deallocate(Di, Dsr)

call CSQ_fun(Mx, My, BsrTot0, BsrTot1, BsrTot)
call RyR_p_fun(Mx, My, P0, P1, alpha)

! initial conditions
Ca = 0.08d0
CaSr = 900d0
CabCM = Ca * BCAM / (koffCM/konCM + Ca)
CabSR = Ca * BSR / (koffSR/konSR + Ca)
CaTnC = Ca * BT / (koff/kon + Ca)
CabRhod = Ca * BTRhod / (koffRhod/konRhod + Ca)

states_RyR = 1
states_LCC = 1

! Initial condition for Ca
! open(unit=50, file='../data/final-state-Ca.data', status='old')
! read(50,*) Ca
! read(50,*) CaSr
! read(50,*) CaTnC
! read(50,*) CabCM
! read(50,*) CabSR
! read(50,*) CabRhod
! close(50)
! 
! open(unit=60, file='../data/final-state-RR-LCC.data', status='old')
! read(60,*) states_RyR
! read(60,*) states_LCC
! close(60)

CaSrTot = CaSr + BsrTot*CaSr / (Ksrtot + CaSr)

! Random number generator
call init_ranmar()

! files
open(unit=2, file='../data/currents-cy.data', status='unknown')
open(unit=10, file='../data/calcium.data', status='unknown')
open(unit=28, file='../data/openRyR-tot.data', status='unknown')
open(unit=30, file='../data/timeRyR-open.data', status='unknown')

tmatrix = 0

!Name of the VTK file, goes from 0000-1999
Calc='Calcxxxx.data'
Cab_name='Cabcxxxx.data'
! check there are enough files
if (int(Nt/nexp) .gt. 9999) then
   write(*,*) 'WARNING: not enough Calcxxxx.data files!!'
end if

!-- DYNAMIC ---
t = 0d0
V = Vrest

casave = 0d0
cabsave = 0d0
isave = 0

do i = 1, Nt
!  POTENTIAL
   Mt = int(t / Ts)
   tt = t - Mt * Ts
   !if (tt .lt. ADP .and. t .lt. 8 * Ts) then
   !if (tt .lt. ADP .and. t .gt. Ts) then
   !if (tt .lt. ADP) then
   !   V = (Vmax - Vrest) * dsqrt(1 - (tt / ADP)**0.5) + Vrest
   !   !V = Vmax ! Square potential
   !else
   !   V = Vrest
   !end if

   if (i .gt. 416666) then
        casave = casave + Ca
        cabsave = cabsave + CabRhod
        isave = isave + 1
        if(mod(i,nexp).eq.1) then
           umiler = int(nplot2/1000)
           nplot_cent = mod(nplot2,1000)
           cent = int(nplot_cent/100)
           nplot_des = mod(nplot_cent,100)
           des = int(nplot_des/10)
           nplot_uni = mod(nplot_des,10)

           Calc(5:5)=char(ichar('0')+umiler)
           Calc(6:6)=char(ichar('0')+cent)
           Calc(7:7)=char(ichar('0')+des)
           Calc(8:8)=char(ichar('0')+nplot_uni)

           Cab_name(5:5)=char(ichar('0')+umiler)
           Cab_name(6:6)=char(ichar('0')+cent)
           Cab_name(7:7)=char(ichar('0')+des)
           Cab_name(8:8)=char(ichar('0')+nplot_uni)

           open(70,file=Calc, status='unknown')
           write(70,*) casave / isave
           close(70)

           open(70,file=Cab_name, status='unknown')
           write(70,*) cabsave / isave
           close(70)

           casave = 0d0
           cabsave = 0d0
           isave = 0
           nplot2=nplot2+1
        end if
   end if

   ! save snapshoots for the films
   !if(mod(i,nfilm).eq.1) then
   !   nplot=nplot+1
   !   write(30000,*) nplot, t
   !   write(1000+nplot,*) Ca
   !   !write(10000+nplot,*) CaSr
   !   close(1000+nplot)
   !   !close(10000+nplot)
   !end if

   if(mod(i,nsave).eq.1) then
      write(10,*) t, sum(Vi*Ca)/sum(Vi), sum(Vsr*Casr)/sum(Vsr), &
                 sum(CaTnC)/n_TnC, sum(Vi*CabCM)/sum(Vi), &
                 sum(Vi*CabSR)/sum(Vi), sum(Vsr*CaSrTot)/sum(Vsr), &
                 sum(Vi*CabRhod)/sum(Vi)
   end if

   bcoeff = Ksrtot+Bsrtot - CaSrTot
   ccoeff = 4d0 * CaSrTot*Ksrtot
   CaSr = 0.5d0 * (-bcoeff + dsqrt(bcoeff*bcoeff + ccoeff))

   !if (mod(i,500) .eq. 0) write(*,*) float(i) / Nt * 100

   ! Laplacian matrix
   call laplacian(dx, dx, Mx, My, Ca, Di_full, lapCa)
   call laplacian(dx, dx, Mx, My, CaSr, Dsr_full, lapCaSr)

!  INIT MATRIX FIELDS
   J_currents_Ca = 0
   J_currents_CaSr = 0
   J_currents_CabTnC = 0
   J_currents_CabCM = 0
   J_currents_CabRhod = 0
   J_currents_CabSR = 0

   JNaCaT = 0
   JupT = 0
   JrelT = 0
   JCaLT = 0
   JbuffCT = 0

!  -----    Jup CURRENT    ------
!     concentrations: Ca and CaSr
!     Fluxes Ca: -Jup
!     Fluxes Sr: Jup

!$omp parallel do shared(Ca, CasR, J_currents_Ca, J_currents_CaSr, &
!$omp ViVsr, gup, Ki, Ksr) private(j, k, Jup)
   do j = 1, Mx
      do k = 1, My
         call current_Jup(Ca(j,k), CaSr(j,k), gup, Ki, Ksr, Jup)
         ! CURRENTS
         J_currents_Ca(j,k) = J_currents_Ca(j,k) - Jup 
         J_currents_CaSr(j,k) = J_currents_CaSr(j,k) + ViVsr(j,k) * Jup

         JupT = JupT - Vi(j,k)*Jup
      end do
   end do
!$omp end parallel do

! -----   DYNAMIC ON THE SURFACE    ------
!     concentrations: Ca
!     Fluxes Ca: JNaCa

!$omp parallel do shared(Ca, J_currents_Ca,  &
!$omp V, F, R, Temp, Na0, Nai, Ca0, KmNa0, KmCai, &
!$omp KmCa0, KmNai, gNaCa, Kda, nu, ksat) private(j, JNaCa)
   do j = 1, nNCX
      call current_JNaCa(V, F, R, Temp, Na0, Nai, Ca0, KmNa0, KmCai, KmCa0, &
      KmNai, gNaCa, Kda, nu, ksat, Ca(ind_NCX(j,1), ind_NCX(j,2)), JNaCa)
      ! CURRENTS
      J_currents_Ca(ind_NCX(j,1), ind_NCX(j,2)) = J_currents_Ca(ind_NCX(j,1), ind_NCX(j,2)) + JNaCa 
      JNaCaT = JNaCaT + Vi(ind_NCX(j,1), ind_NCX(j,2)) * JNaCa
   end do
!$omp end parallel do

   do j = 1, n_Back
      J_currents_Ca(ind_Back(j,1), ind_Back(j,2)) = J_currents_Ca(ind_Back(j,1), ind_Back(j,2)) + 0.0975d0 
      JNaCaT = JNaCaT + Vi(ind_Back(j,1), ind_Back(j,2)) * 0.0975d0
   end do

!  -----    DYNAMIC FOR RyR    -----
!     concentrations: Ca and CaSr
!     Fluxes Ca: Jrel
!     Fluxes Sr: -Jrel
   
!$omp parallel do shared(Ca, CaSr, J_currents_Ca, J_currents_CaSr, &
!$omp ViVsr, states_RyR, ind_RR, koc, ki1i2, kic, kio, ka, kb, &
!$omp kco, ki2i1, koi, kci, grel) private(j, Jrel)
   do j = 1, nRyRtot
!     Update states of RyR
      call get_ORyR(dt, nRyRPixel(j), koc, ki1i2, kic, kio, k2, k7, k6, k11, k1d, k3d, k8d, &
         Ca(ind_RR(j,1), ind_RR(j,2)),  CaSr(ind_RR(j,1), ind_RR(j,2)), &
         BsrTot(ind_RR(j,1), ind_RR(j,2))*Ksrtot / (Ksrtot + CaSr(ind_RR(j,1), ind_RR(j,2))), &
         alpha(ind_RR(j,1), ind_RR(j,2)), states_RyR(j,:), ORyR)
!     Currents
      call current_Jrel(grel, ORyR/nRyRPixel(j), CaSr(ind_RR(j,1), ind_RR(j,2)), & 
         Ca(ind_RR(j,1), ind_RR(j,2)), Jrel)
      ! update CURRENTS
      J_currents_Ca(ind_RR(j,1), ind_RR(j,2)) = J_currents_Ca(ind_RR(j,1), ind_RR(j,2)) &
         + Jrel 
      J_currents_CaSr(ind_RR(j,1), ind_RR(j,2)) = J_currents_CaSr(ind_RR(j,1), ind_RR(j,2)) &
         - ViVsr(ind_RR(j,1), ind_RR(j,2)) * Jrel 


      JrelT = JrelT + Vi(ind_RR(j,1), ind_RR(j,2)) * Jrel
   end do
!$omp end parallel do

!  -----    DYNAMIC LCC -----
!     concentrations: Ca and CaSr
!     Fluxes Ca: -JCal

!$omp parallel do shared(Ca, J_currents_Ca, &
!$omp states_LCC, ind_LCC, V, KLCC, F, R, Temp, gCaL) &
!$omp private(j, JCaL, OLCC)
   do j = 1, nLCCtot
!     Update states of LCC
      call get_OLCC(dt, nLCC, V, KLCC, Ca(ind_LCC(j,1), ind_LCC(j,2)), states_LCC(j,:), OLCC)
      call current_JCaL(V, F, R, Temp, gCaL, Ca0, Ca(ind_LCC(j,1), ind_LCC(j,2)), OLCC, JCaL)
      ! update CURRENTS
      J_currents_Ca(ind_LCC(j,1), ind_LCC(j,2)) = J_currents_Ca(ind_LCC(j,1), ind_LCC(j,2)) &
         - JCaL 
      JCaLT = JCaLT - Vi(ind_LCC(j,1), ind_LCC(j,2)) * JCaL
   end do
!$omp end parallel do

! ----------   BUFFERS  -------------

!  DYNAMIC OF TnC BUFFER

!$omp parallel do shared(Ca, CaTnC, J_currents_Ca, J_currents_CabbTnC, ind_TnC, &
!$omp kon, koff, BT) &
!$omp private(j, JbuffC)
   do j = 1, n_TnC
      call current_TnC(Ca(ind_TnC(j,1), ind_TnC(j,2)), CaTnC(ind_TnC(j,1), ind_TnC(j,2)), &
         kon, koff, BT, JbuffC)
      J_currents_CabTnC(ind_TnC(j,1), ind_TnC(j,2)) = J_currents_CabTnC(ind_TnC(j,1), ind_TnC(j,2)) &
         + JbuffC
      JbuffCT = JbuffCT - JbuffC
   end do
!$omp end parallel do

!  buffers Ca: Rhod2

!$omp parallel do shared(Ca, CabSR, J_currents_CabRhod, BTRhod, konRhod, koffRhod) &
!$omp private(j,k)
   do j = 1, Mx
      do k = 1, My
         call current_TnC(Ca(j,k), CabRhod(j,k), konRhod, koffRhod, BTRhod, J_currents_CabRhod(j,k))
      end do
   end do
!$omp end parallel do

!  buffers Ca: SR & calmodulina

!$omp parallel do shared(Ca, CabCM, J_currents_CabCM, BCAM, konCM, koffCM) &
!$omp private(j,k)
   do j = 1, Mx
      do k = 1, My
         call current_TnC(Ca(j,k), CabCM(j,k), konCM, koffCM, BCAM, J_currents_CabCM(j,k))
      end do
   end do
!$omp end parallel do

!$omp parallel do shared(Ca, CabSR, J_currents_CabSR, BSR, konSR, koffSR) &
!$omp private(j,k)
   do j = 1, Mx
      do k = 1, My
         call current_TnC(Ca(j,k), CabSR(j,k), konSR, koffSR, BSR, J_currents_CabSR(j,k))
      end do
   end do
!$omp end parallel do

   J_currents_Ca = J_currents_Ca - J_currents_CabCM - J_currents_CabSR - &
                   J_currents_CabTnC - J_currents_CabRhod

   if (mod(i,nsave) .eq. 0) then
   ! ---   save currents ---
      write(2,*) t, JNaCaT, JrelT, JupT, JCaLT, JbuffCT

      do j = 1, nRyRtot
         do k = 1, nRyRPixel(j)
            if (states_RyR(j,k) .ne. 2) then
               if (tmatrix(j,k,2)-tmatrix(j,k,1) .ne. 0) then
                  write(30,*) j, k, tmatrix(j,k,2), tmatrix(j,k,1)
               end if
               tmatrix(j,k,:) = i
            else
               tmatrix(j,k,2) = i
            end if
         end do
      end do

      ! open RyRs
      openRyRsi = 0
      do j = 1, nRyRtot
         do l = 1, nRyRPixel(j)
            openRyRsi(states_RyR(j,l)) = openRyRsi(states_RyR(j,l)) + 1
         end do
      end do
      write(28,*) t, openRyRsi
   end if


! ----------   INTEGRATION  -------------

!$omp parallel do shared(Ca, CaSr, CabCM, J_currents_Ca, J_currents_CaSr, &
!$omp J_currents_CabCM) private(k, j)
   do j = 1, Mx
      do k = 1, My
         Ca(j,k) = Ca(j,k) + dt * (J_currents_Ca(j,k) + lapCa(j,k)/Vi(j,k))
         CaSrTot(j,k) = CaSrTot(j,k) + dt * (J_currents_CaSr(j,k) + lapCaSr(j,k)/ Vsr(j,k))
         CaTnC(j,k) = CaTnC(j,k) + dt * J_currents_CabTnC(j,k)
         CabSR(j,k) = CabSR(j,k) + dt * J_currents_CabSR(j,k)
         CabCM(j,k) = CabCM(j,k) + dt * J_currents_CabCM(j,k)
         CabRhod(j,k) = CabRhod(j,k) + dt * J_currents_CabRhod(j,k)
      end do
   end do
!$omp end parallel do

   !if (mod(i,Nt/20) .eq. 0 .or. i .eq. Nt) then
   if (mod(i-1,int(Ts/dt)) .eq. 0 .or. i .eq. Nt) then
      open(unit=50, file='../data/final-state-Ca.data', status='unknown')
      open(unit=60, file='../data/final-state-RR-LCC.data', status='unknown')
      write(50,*) Ca
      write(50,*) CaSr
      write(50,*) CaTnC
      write(50,*) CabCM
      write(50,*) CabSR
      write(50,*) CabRhod
      write(60,*) states_RyR
      write(60,*) states_LCC
      close(50)
      close(60)
   end if

   t = t + dt

end do

close(2)
close(3)
close(10)
close(28)
close(30)

end program

! -------------

subroutine laplacian(dx, dy, Nx, Ny, aux, D, lap)
implicit none
integer, intent(in) :: Nx, Ny
double precision, intent(in) :: dx, dy, aux(Nx,Ny), D(0:Nx+1,0:Ny+1)
double precision, intent(out) :: lap(Nx,Ny)
integer :: i, j
double precision :: gradf(2,2), Dm(2,2), lapx, lapy
double precision :: fun(0:Nx+1,0:Ny+1)

fun(1:Nx, 1:Ny) = aux

! Non-flux conditions
fun(0,:) = fun(1,:)
fun(Nx+1,:) = fun(Nx,:)
fun(:,0) = fun(:,1)
fun(:,Ny+1) = fun(:,Ny)

fun(0,0) = fun(1,1)
fun(Nx+1,0) = fun(Nx,1)
fun(0,Ny+1) = fun(1,Ny)
fun(Nx+1,Ny+1) = fun(Nx,Ny)

!lap = 0

!$omp parallel do shared(fun, D, dx, dy, lap) private(i, j, gradf, Dm, lapx, lapy)
do i = 1, Nx
   do j = 1, Ny
      call grad1D(fun(i,j), fun(i+1,j), gradf(1,1)) 
      call grad1D(fun(i-1,j), fun(i,j), gradf(1,2)) 
      call grad1D(fun(i,j), fun(i,j+1), gradf(2,1)) 
      call grad1D(fun(i,j-1), fun(i,j), gradf(2,2)) 

      call middle_point(D(i,j), D(i+1,j), Dm(1,1)) 
      call middle_point(D(i-1,j), D(i,j), Dm(1,2)) 
      call middle_point(D(i,j), D(i,j+1), Dm(2,1)) 
      call middle_point(D(i,j-1), D(i,j), Dm(2,2)) 
     
      call grad1D(gradf(1,2) * Dm(1,2), gradf(1,1) * Dm(1,1), lapx) 
      call grad1D(gradf(2,2) * Dm(2,2), gradf(2,1) * Dm(2,1), lapy) 

      lap(i,j) = lapx / (dx*dx) + lapy / (dy*dy)
   end do
end do
!$omp end parallel do

end subroutine

! -------------

subroutine grad1D(f1, f2, gradf)
implicit none
double precision, intent(in) :: f1, f2
double precision, intent(out) :: gradf

gradf = f2 - f1

end subroutine

! -------------

subroutine middle_point(f1, f2, fm)
implicit none
double precision, intent(in) :: f1, f2
double precision, intent(out) :: fm

fm = 0.5d0 * (f2 + f1)

end subroutine

! -------------

subroutine diffusion_cytosol(vfrac, D)
implicit none
double precision, intent(in) :: vfrac
double precision, intent(out) :: D
! Internal variable
double precision :: D0, D1, m

D0 = 0.25d0 ! micro m**2 / ms - diffusion when the volume fraction of SR is 0
D1 = 0.09d0 ! micro m**2 / ms - diffusion when the volume fraction of SR is 1

m = D1 - D0 ! slope

D = m * vfrac + d0 ! diffusion interpolation

end subroutine

! -------------

subroutine diffusion_SR(vfrac, D)
implicit none
double precision, intent(in) :: vfrac
double precision, intent(out) :: D
! Internal variable
double precision :: D0, D1, m

D0 = 0.09d0 ! micro m**2 / ms - diffusion when the volume fraction of SR is 0
D1 = 0.25d0 ! micro m**2 / ms - diffusion when the volume fraction of SR is 1

m = D1 - D0 ! slope

D = m * vfrac + d0 ! diffusion interpolation

end subroutine

! ----------------

subroutine diff_bc(Nx, Ny, D_aux, D)
implicit none
integer, intent(in) :: Nx, Ny
double precision, intent(in) :: D_aux(Nx,Ny)
double precision, intent(out) :: D(0:Nx+1,0:Ny+1) 

D(1:Nx, 1:Ny) = D_aux

! Non-flux conditions. EDGES
D(0,:) = D(1,:)
D(Nx+1,:) = D(Nx,:)
D(:,0) = D(:,1)
D(:,Ny+1) = D(:,Ny)

! Non-flux conditions. VERTEXS
D(0,0) = D(1,1)
D(Nx+1,0) = D(Nx,1)
D(0,Ny+1) = D(1,Ny)
D(Nx+1,Ny+1) = D(Nx,Ny)

end subroutine

! ----------------

subroutine CSQ_fun(Nx, Ny, BsrTot0, BsrTot1, BsrTot)
implicit none
integer, intent(in) :: Nx, Ny
double precision, intent(in) :: BsrTot0, BsrTot1
double precision, intent(out) :: BsrTot(Nx,Ny) 
integer :: i, j, d1, d2, d3, d4, d

do i = 1, Nx
   do j = 1, Ny
      d1 = i
      d2 = Nx-i
      d3 = j
      d4 = Ny-j
      d = min(d1, d2, d3, d4)
      BsrTot(i,j) = BsrTot0 + (BsrTot1-BsrTot0)/(Nx/2) * d
   end do
end do

end subroutine

! ----------------

subroutine RyR_p_fun(Nx, Ny, P0, P1, alpha)
implicit none
integer, intent(in) :: Nx, Ny
double precision, intent(in) :: P0, P1
double precision, intent(out) :: alpha(Nx,Ny) 
integer :: i, j, d1, d2, d3, d4, d
double precision :: Pij

do i = 1, Nx
   do j = 1, Ny
      d1 = i
      d2 = Nx-i
      d3 = j
      d4 = Ny-j
      d = min(d1, d2, d3, d4)
      Pij = P0 + (P1-P0)/(Nx/2) * d
      alpha(i,j) = 1 + 15 * Pij**12 / (Pij**12 + 1.15**12)
      alpha(i,j) = alpha(i,j)/8d0
   end do
end do

end subroutine
