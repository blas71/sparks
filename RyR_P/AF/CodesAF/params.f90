Nt = 3000000
dt = 0.012d0 ! ms
Ts = 1000d0 ! ms

! Set number of RyR and LCC
nRyR = 9
nRyR_membrane = 7
nLCC = 5

!Di = 0.08d0 ! microns**2 / ms
!DSR = 0.08d0 ! microns**2 / ms

Vrest = -85d0 ! mV
Vmax = 10d0 ! mV
ADP = 350d0 ! ms

Temp = 308d0 ! K
F = 96485e-21 ! C / mini mol !!! 1 mini mol = 10**21 * 1microM * 1micro m**3
R = 8.31e-18 ! mV * C / (K mini mol)

Ksrtot = 2000d0 ! micro M
Bsrtot0 = 6000d0 ! micro M
Bsrtot1 = 6000d0 ! micro M

P0 = 1.2d0
P1 = 0.9d0

kon = 0.0327d0 / 10d0 ! 1 / (micro M * ms)
koff = 0.0196d0 / 10d0 ! 1 /  ms
BT_val = 70d0 ! micro M
BCAM = 24d0 ! micro M
BSR = 0.5d0 * 47d0 ! micro M
BTRhod = 5d0

konCM = 0.03d0 / 10d0 ! 1 / (micro M * ms)
koffCM = 0.2d0 / 10d0 ! 1 /  ms
konSR = 0.1d0 / 10d0 ! 1 / (micro M * ms)
koffSR = 0.06d0 / 10d0 ! 1 /  ms
konRhod = 0.069d0
koffRhod = 0.13d0

gNaCa = 9d0 ! micro M / ms
Ca0 = 1800d0 ! micro M
Na0 = 136000d0 ! micro M
Nai = 10000d0 ! micro M
Kda = 0.1d0 ! micro M
ksat = 0.27d0

nu = 0.35d0
KmNa0 = 87500d0 ! micro M
KmCa0 = 1300d0 ! micro M
KmNai = 12300d0 ! micro M
KmCai = 3.6d0 ! micro M

gup = 0.1d0 ! micro M / ms
Ki = 0.123d0 ! micro M
Ksr = 1300d0 ! micro M

grel = 80d0 / 9d0 ! 1 / ms
koc = 0.027d0 ! 1 / ms
ki1i2 = koc ! 1 / ms
kic = 1d0/200d0 ! 1 / ms
kio = kic ! 1 / ms

k2 = 0.1 ! 1 / (micro M ^2 * ms)
k6 = 0.033 ! 1 / (micro M * ms)
k7 = 0.001666 ! 1 / (micro M * ms)
k11 = 0.00667 ! 1 / (micro M ^2 * ms)
k1d = 2000d0 ! micro M
k3d = k1d * k11 / k2
k8d = k1d * k6 / k7

gCaL = 1.5d0 * 5.24d15 / 4d0 ! mini mol / (C * ms)
KLCC = 5d0 ! micro M

! volume factor SR
vfrac_ci = 0.005d0
vfrac_SR = 0.05d0 

! save parameters
nsave = 28
nplot = 1
nplot2 = 1
nfilm = 142
nopen = 0
nexp = 833
