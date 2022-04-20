subroutine current_Jup(Ca, CaSr, gup, Ki, Ksr, Jup)
implicit none
double precision, intent(in) :: Ca, CaSr, gup, Ki, Ksr
double precision, intent(out) :: Jup

Jup = gup * ((Ca / Ki)**2 - (CaSr / Ksr)**2) / (1d0 + (Ca / Ki)**2 + (CaSr / Ksr)**2)

!Jup=0.

end subroutine

! -------

subroutine current_JNaCa(V, F, R, Temp, Na0, Nai, Ca0, KmNa0, KmCai, KmCa0, KmNai, gNaCa, Kda, nu, &
      ksat, Ca, JNaCa)
implicit none
double precision, intent(in) :: V, F, R, Temp, Na0, Nai, Ca0, KmNa0, KmCai, KmCa0, KmNai
double precision, intent(in) :: gNaCa, Kda, nu, ksat, Ca
double precision, intent(out) :: JNaCa
double precision :: z, Scs

z = V * F / (R * Temp)

Scs = Na0**3 * Ca + KmNa0**3 * Ca * (1 + Ca / KmCai) + KmCa0 * Nai**3 + &
   Nai**3 * Ca0 + KmCai * Na0**3 * (1 + (Nai / KmNai)**3)

JNaCa = gNaCa / (1 + (Kda / Ca)**3) * (exp(nu * z) * Nai**3 * Ca0 - exp((nu - 1 ) &
   * z) * Na0**3 * Ca) / (Scs * (1 + ksat * exp((nu - 1) * z)))

end subroutine

! -------

subroutine current_Jrel(grel, ORyR, CaSr, Ca, Jrel)
implicit none
double precision, intent(in) :: grel, ORyR, CaSr, Ca
double precision, intent(out) :: Jrel

Jrel = grel * ORyR * (CaSr - Ca)

!Jrel = 0

end subroutine

! -------

subroutine current_JCaL(V, F, R, Temp, gCaL, Ca0, Ca, OLCC, JCaL)
implicit none
double precision, intent(in) :: V, F, R, Temp, gCaL, Ca0, Ca, OLCC
double precision, intent(out) :: JCaL
double precision :: z, zm

z = V * F / (R * Temp)
zm = 0.341d0 * z * F

JCaL = gCaL * 4 * OLCC * zm * (exp(2*z)*Ca - Ca0) / (exp(2*z) - 1)

!JCaL=0.

end subroutine

! -----

subroutine current_TnC(Ca, CaTnC, kon, koff, BT, JbuffC)
implicit none
double precision, intent(in) :: Ca, CaTnC, kon, koff, BT
double precision, intent(out) :: JbuffC

JbuffC = kon * Ca * (BT - CaTnc) - koff * CaTnC

!JbuffC=0.

end subroutine
   
