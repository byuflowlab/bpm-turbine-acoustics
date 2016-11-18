! Subroutine to use the BPM equations for turbine acoustics

! Directivity function for high-frequency noise
! (not high-angle separation or turbulent inflow noise);
! becomes inaccurate for theta_e approaching 180 deg
subroutine Dhfunc(theta_e,phi_e,M,Mc,Dh)
  implicit none
  integer, parameter :: dp = kind(0.d0)
  ! in
  real(dp), intent(in) :: theta_e,phi_e,M,Mc
  ! out
  real(dp), intent(out) :: Dh
  ! local
  real(dp) :: theta,phi,pi
  intrinsic sin
  intrinsic cos
  ! constants
  pi = 3.1415926535897932_dp

  theta = theta_e*pi/180.0_dp
  phi = phi_e*pi/180.0_dp

  Dh = (2.0_dp*(sin(theta/2.0_dp))**2*(sin(phi))**2)/((1.0_dp&
  +M*cos(theta))*(1.0_dp+(M-Mc)*cos(theta))**2)

end subroutine Dhfunc

! Directivity function for low-frequency noise
subroutine Dlfunc(theta_e,phi_e,M,Dl)
  implicit none
  integer, parameter :: dp = kind(0.d0)
  ! in
  real(dp), intent(in) :: theta_e,phi_e,M
  ! out
  real(dp), intent(out) :: Dl
  ! local
  real(dp) :: theta,phi,pi
  intrinsic sin
  intrinsic cos
  ! constants
  pi = 3.1415926535897932_dp

  theta = theta_e*pi/180.0_dp
  phi = phi_e*pi/180.0_dp

  Dl = ((sin(theta))**2*(sin(phi))**2)/(1.0_dp+M*cos(theta))**4

end subroutine Dlfunc

! Spectral Function A
subroutine Afunc(ain,Re,Aspec)
  implicit none
  integer, parameter :: dp = kind(0.d0)
  ! in
  real(dp), intent(in) :: ain,Re
  ! out
  real(dp), intent(out) :: Aspec
  ! local
  real(dp) :: a,Amin,Amax,a0,Amin0,Amax0,AR
  intrinsic log10
  intrinsic abs
  intrinsic sqrt

  a = abs(log10(ain))

  ! Calculating Amin
  if(a < 0.204_dp) then
    Amin = sqrt(67.552_dp-886.788_dp*a**2)-8.219_dp
  else if (a >= 0.204_dp .and. a <= 0.244_dp) then
    Amin = -32.665_dp*a+3.981_dp
  else
    Amin = -142.795_dp*a**3+103.656_dp*a**2-57.757_dp*a+6.006_dp
  end if

  ! Calculating Amax
  if (a < 0.13_dp) then
    Amax = sqrt(67.552_dp-886.788_dp*a**2)-8.219_dp
  else if(a >= 0.13_dp .and. a <= 0.321_dp) then
    Amax = -15.901_dp*a+1.098_dp
  else
    Amax = -4.669_dp*a**3+3.491_dp*a**2-16.699_dp*a+1.149_dp
  end if

  ! Calculating a0
  if (Re < 9.52e4_dp) then
    a0 = 0.57_dp
  else if (Re >= 9.52e4_dp .and. Re <= 8.57e5_dp) then
    a0 = -9.57e-13_dp*(Re-8.57e5_dp)**2+1.13_dp
  else
    a0 = 1.13_dp
  end if

  ! Calculating Amin(a0)
  if (a0 < 0.204_dp) then
    Amin0 = sqrt(67.552_dp-886.788_dp*a0**2)-8.219_dp
  else if (a0 >= 0.204 .and. a0 <= 0.244) then
    Amin0 = -32.665_dp*a0+3.981_dp
  else
    Amin0 = -142.795_dp*a0**3+103.656_dp*a0**2-57.757_dp*a0+6.006_dp
  end if

  ! Calculating Amax(a0)
  if (a0 < 0.13_dp) then
    Amax0 = sqrt(67.552_dp-886.788_dp*a0**2)-8.219_dp
  else if (a0 >= 0.13_dp .and. a0 <= 0.321_dp) then
    Amax0 = -15.901_dp*a0+1.098_dp
  else
    Amax0 = -4.669_dp*a0**3+3.491_dp*a0**2-16.699_dp*a0+1.149_dp
  end if

  AR = (-20.0_dp-Amin0)/(Amax0-Amin0)

  Aspec = Amin+AR*(Amax-Amin)

end subroutine Afunc

! Spectral Function B
subroutine Bfunc(bin,Re,Bspec)
  implicit none
  integer, parameter :: dp = kind(0.d0)
  ! in
  real(dp), intent(in) :: bin,Re
  ! out
  real(dp), intent(out) :: Bspec
  ! local
  real(dp) :: b,Bmin,Bmax,b0,Bmin0,Bmax0,BR
  intrinsic log10
  intrinsic abs
  intrinsic sqrt

  b = abs(log10(bin))

  ! Calculating Bmin
  if (b < 0.13_dp) then
    Bmin = sqrt(16.888_dp-886.788_dp*b**2)-4.109_dp
  else if (b >= 0.13_dp .and. b <= 0.145_dp) then
    Bmin = -83.607_dp*b+8.138_dp
  else
    Bmin = -817.810_dp*b**3+355.210_dp*b**2-135.024_dp*b+10.619_dp
  end if

  ! Calculating Bmax
  if (b < 0.10_dp) then
    Bmax = sqrt(16.888_dp-886.788_dp*b**2)-4.109_dp
  else if (b >= 0.10_dp .and. b <= 0.187_dp) then
    Bmax = -31.330_dp*b+1.854_dp
  else
    Bmax = -80.541_dp*b**3+44.174_dp*b**2-39.381_dp*b+2.344_dp
  end if

  ! Calculating b0
  if (Re < 9.52e4_dp) then
    b0 = 0.30_dp
  else if (Re >= 9.52e4_dp .and. Re <= 8.57e5_dp) then
    b0 = -4.48e-13_dp*(Re-8.57e5_dp)**2+0.56_dp
  else
    b0 = 0.56_dp
  end if

  ! Calculating Bmin(b0)
  if (b0 < 0.13_dp) then
    Bmin0 = sqrt(16.888_dp-886.788_dp*b0**2)-4.109_dp
  else if (b0 >= 0.13_dp .and. b0 <= 0.145_dp) then
    Bmin0 = -83.607_dp*b0+8.138_dp
  else
    Bmin0 = -817.810_dp*b0**3+355.210_dp*b0**2-135.024_dp*b0+10.619_dp
  end if

  ! Calculating Bmax(b0)
  if (b0 < 0.10_dp) then
    Bmax0 = sqrt(16.888_dp-886.788_dp*b0**2)-4.109_dp
  else if (b0 >= 0.10_dp .and. b0 <= 0.187_dp) then
    Bmax0 = -31.330_dp*b+1.854_dp
  else
    Bmax0 = -80.541_dp*b0**3+44.174_dp*b0**2-39.381_dp*b0+2.344_dp
  end if

  BR = (-20_dp-Bmin0)/(Bmax0-Bmin0)

  Bspec =  Bmin+BR*(Bmax-Bmin)

end subroutine Bfunc

subroutine G1func(e,G1)
  implicit none
  integer, parameter :: dp = kind(0.d0)
  ! in
  real(dp), intent(in) :: e
  ! out
  real(dp), intent(out) :: G1
  ! local
  intrinsic log10
  intrinsic sqrt

  if (e <= 0.5974_dp) then
    G1 = 39.8_dp*log10(e)-11.12_dp
  else if (e <= 0.8545_dp .and. e > 0.5974_dp) then
    G1 = 98.409_dp*log10(e)+2.0_dp
  else if (e <= 1.17_dp .and. e > 0.8545_dp) then
    G1 = sqrt(2.484_dp-506.25_dp*(log10(e))**2)-5.076_dp
  else if (e <= 1.674_dp .and. e > 1.17_dp) then
    G1 = -98.409_dp*log10(e)+2.0_dp
  else
    G1 = -39.8_dp*log10(e)-11.12_dp
  end if

end subroutine G1func

subroutine G2func(d,G2)
  implicit none
  integer, parameter :: dp = kind(0.d0)
  ! in
  real(dp), intent(in) :: d
  ! out
  real(dp), intent(out) :: G2
  ! local
  intrinsic log10

  if (d <= 0.3237_dp) then
    G2 = 77.852_dp*log10(d)+15.328_dp
  else if (d <= 0.5689_dp .and. d > 0.3237_dp) then
    G2 = 65.188_dp*log10(d)+9.125_dp
  else if (d <= 1.7579_dp .and. d > 0.5689_dp) then
    G2 = -114.052_dp*(log10(d))**2
  else if (d <= 3.0889_dp .and. d > 1.7579_dp) then
    G2 = -65.188_dp*log10(d)+9.125_dp
  else
    G2 = -77.852_dp*log10(d)+15.328_dp
  end if

end subroutine G2func

subroutine G3func(alpha,G3)
  implicit none
  integer, parameter :: dp = kind(0.d0)
  ! in
  real(dp), intent(in) :: alpha
  ! out
  real(dp), intent(out) :: G3

  G3 = 171.04_dp-3.03_dp*alpha

end subroutine G3func

subroutine G4func(hdav,psi,G4)
  implicit none
  integer, parameter :: dp = kind(0.d0)
  ! in
  real(dp), intent(in) :: hdav,psi
  ! out
  real(dp), intent(out) :: G4
  ! local
  intrinsic log10

  if (hdav <= 5.0_dp) then
    G4 = 17.5_dp*log10(hdav)+157.5_dp-1.114_dp*psi
  else
    G4 = 169.7_dp-1.114_dp*psi
  end if

end subroutine G4func

subroutine G5func(hdav,psi,StSt_peak,G5)
  implicit none
  integer, parameter :: dp = kind(0.d0)
  ! in
  real(dp), intent(in) :: hdav,psi,StSt_peak
  ! out
  real(dp), intent(out) :: G5
  ! local
  real(dp) :: eta,mu,m,eta_0,k,G14,hdav_prime,eta0,mu0,m0,eta_00,k0,G0
  intrinsic log10
  intrinsic sqrt

  ! finding G5 at phi = 14 deg
  eta = log10(StSt_peak)

  if (hdav < 0.25_dp) then
    mu = 0.1221_dp
  else if (hdav < 0.62_dp .and. hdav >= 0.25_dp) then
    mu = -0.2175_dp*hdav+0.1755_dp
  else if (hdav < 1.15_dp .and. hdav >= 0.62_dp) then
    mu = -0.308_dp*hdav+0.0596_dp
  else
    mu = 0.0242_dp
  end if

  if(hdav <= 0.02_dp) then
    m = 0.0_dp
  else if (hdav <= 0.5_dp .and. hdav > 0.02_dp) then
    m = 68.724_dp*(hdav)-1.35_dp
  else if (hdav <= 0.62_dp .and. hdav > 0.5_dp) then
    m = 308.475_dp*hdav-121.23_dp
  else if (hdav <= 1.15_dp .and. hdav > 0.62_dp) then
    m = 224.811_dp*hdav-69.35_dp
  else if (hdav <= 1.2 .and. hdav > 1.15) then
    m = 1583.28_dp*hdav-1631.59_dp
  else
    m = 268.344_dp
  end if

  eta_0 = -sqrt((m**2*mu**4)/(6.25_dp+m**2*mu**2))
  k = 2.5_dp*sqrt(1.0_dp-(eta_0/mu)**2)-2.5_dp-m*eta_0

  if (eta < eta_0) then
    G14 = m*eta+k
  else if (eta < 0.0_dp .and. eta >= eta_0) then
    G14 = 2.5_dp*sqrt(1.0_dp-(eta/mu)**2)-2.5_dp
  else if (eta < 0.03616_dp .and. eta >= 0.0_dp) then
    G14 = sqrt(1.5625_dp-1194.99_dp*eta**2)-1.25_dp
  else
    G14 = -155.543_dp*eta+4.375_dp
  end if

  ! finding G5 at psi = 0 deg
  hdav_prime = 6.724_dp*hdav**2-4.019_dp*hdav+1.107_dp

  eta0 = log10(StSt_peak)

  if (hdav_prime < 0.25_dp) then
    mu0 = 0.1221_dp
  else if (hdav_prime < 0.62_dp .and. hdav_prime >= 0.25_dp) then
    mu0 = -0.2175_dp*hdav_prime+0.1755_dp
  else if (hdav_prime < 1.15_dp .and. hdav_prime >= 0.62_dp) then
    mu0 = -0.308_dp*hdav_prime+0.0596_dp
  else
    mu0 = 0.0242_dp
  end if

  if (hdav_prime <= 0.02_dp) then
    m0 = 0.0_dp
  else if (hdav_prime <= 0.5_dp .and. hdav_prime > 0.02_dp) then
    m0 = 68.724_dp*hdav_prime-1.35_dp
  else if (hdav_prime <= 0.62_dp .and. hdav_prime > 0.5_dp) then
    m0 = 308.475_dp*hdav_prime-121.23_dp
  else if (hdav_prime <= 1.15_dp .and. hdav_prime > 0.62_dp) then
    m0 = 224.811_dp*hdav_prime-69.35_dp
  else if (hdav_prime <= 1.2_dp .and. hdav_prime > 1.15_dp) then
    m0 = 1583.28_dp*hdav_prime-1631.59_dp
  else
    m0 = 268.344_dp
  end if

  eta_00 = -sqrt((m0**2*mu0**4)/(6.25_dp+m0**2*mu0**2))
  k0 = 2.5_dp*sqrt(1.0_dp-(eta_00/mu0)**2)-2.5_dp-m0*eta_00

  if (eta0 < eta_00) then
    G0 = m0*eta0+k
  else if (eta0 < 0.0_dp .and. eta0 >= eta_00) then
    G0 = 2.5_dp*sqrt(1.0_dp-(eta0/mu0)**2)-2.5_dp
  else if (eta0 < 0.03616_dp .and. eta0 >= 0.0_dp) then
    G0 = sqrt(1.5625_dp-1194.99_dp*eta0**2)-1.25_dp
  else
    G0 = -155.543_dp*eta0+4.375_dp
  end if

  G5 = G0+0.0714_dp*psi*(G14-G0)

end subroutine G5func

! Turbulent Boundary Layer Trailing Edge Noise
subroutine TBLTEfunc(f,V,L,c,r,theta_e,phi_e,alpha,nu,conv,trip,TBLTE)
  implicit none
  integer, parameter :: dp = kind(0.d0)
  ! in
  real(dp), intent(in) :: f,V,L,c,r,theta_e,phi_e,alpha,nu,conv
  logical, intent(in) :: trip
  ! out
  real(dp), intent(out) :: TBLTE
  ! local
  real(dp) :: c0,M,Mc,Re,d0,d0_d,ds_d,dpr,dp_d,Dh,Dl,Stp,Sts,St1,St2,St_bar
  real(dp) :: apre,asuc,bang,gamma,gamma0,beta,beta0,K1,K2,Re_dp,DeltaK1
  real(dp) :: Ap,As,B,SPLp,SPLs,SPLa,A,rc
  intrinsic log10
  intrinsic sqrt
  intrinsic abs
  ! constants
  c0 = 343.2_dp ! speed of sound (m/s)
  M = V/c0
  Mc = conv*M
  Re = (V*c)/nu

  if (trip .eqv. .false.) then
    ! UNTRIPPED boundary layer at 0 deg- thickness, displacement thickness
    d0 = c*(10.0_dp**(1.6569_dp-0.9045_dp*log10(Re)+0.0596_dp*(log10(Re))**2))
    d0_d = c*(10.0_dp**(3.0187_dp-1.5397_dp*log10(Re)+0.1059_dp*(log10(Re))**2))
  else
    ! TRIPPED boundary layer at 0 deg- displacement thickness
    d0 = c*(10.0_dp**(1.892_dp-0.9045_dp*log10(Re)+0.0596_dp*(log10(Re))**2))
    if (Re <= 0.3e6_dp) then
      d0_d = c*0.0601_dp*Re**(-0.114_dp)
    else
      d0_d = c*(10.0_dp**(3.411_dp-1.5397_dp*log10(Re)+0.1059_dp*(log10(Re))**2))
    end if
  end if

  ! boundary layer on pressure side- thickness, displacement thickness
  dpr = d0*(10.0_dp**(-0.04175_dp*alpha+0.00106_dp*alpha**2))
  dp_d = d0_d*(10.0_dp**(-0.0432_dp*alpha+0.00113_dp*alpha**2))

  if (trip .eqv. .false.) then
    ! UNTRIPPED boundary layer on suction side- displacement thickness
    if (alpha <= 7.5_dp .and. alpha >= 0.0_dp) then
      ds_d = d0_d*10.0_dp**(0.0679_dp*alpha)
    else if (alpha <= 12.5_dp .and. alpha > 7.5_dp) then
      ds_d = d0_d*0.0162_dp*10.0_dp**(0.3066_dp*alpha)
    !else if (alpha <= 25.0_dp .and. alpha > 12.5_dp) then
    else
      ds_d = d0_d*52.42_dp*10.0_dp**(0.0258_dp*alpha)
    end if
  else
    ! TRIPPED boundary layer on suction side- displacement thickness
    if (alpha <= 5.0_dp .and. alpha >= 0.0_dp) then
      ds_d = d0_d*10.0_dp**(0.0679_dp*alpha)
    else if (alpha <= 12.5_dp .and. alpha > 5.0_dp) then
      ds_d = d0_d*0.381_dp*10.0_dp**(0.1516_dp*alpha)
    !else if (alpha <= 25.0_dp .and. alpha > 12.5_dp) then
    else
      ds_d = d0_d*14.296_dp*10.0_dp**(0.0258_dp*alpha)
    end if
  end if

  call Dhfunc(theta_e,phi_e,M,Mc,Dh)
  call Dlfunc(theta_e,phi_e,M,Dl)

  Stp = (f*dp_d)/V
  Sts = (f*ds_d)/V

  St1 = 0.02_dp*M**(-0.6_dp)

  if (alpha < 1.33_dp) then
    St2 = 1.0_dp
  else if (alpha < 12.5_dp .and. alpha >= 1.33_dp) then
    St2 = 10.0_dp**(0.0054_dp*(alpha-1.33_dp)**2)
  else
    St2 = 4.72_dp
  end if

  St_bar = (St1+St2)/2.0_dp

  apre = abs(Stp/St1)
  asuc = abs(Sts/St_bar)
  bang = abs(Sts/St2)

  gamma = 27.094_dp*M+3.31_dp
  gamma0 = 23.43_dp*M+4.651_dp
  beta = 72.65_dp*M+10.74_dp
  beta0 = -34.19_dp*M-13.82_dp

  if (Re < 2.47e5_dp) then
    K1 = -4.31_dp*log10(Re)+156.3_dp
  else if (Re >= 2.47e5_dp .and. Re <= 8.0e5_dp) then
    K1 = -9.0_dp*log10(Re)+181.6_dp
  else
    K1 = 128.5_dp
  end if

  if (alpha < (gamma0-gamma)) then
    K2 = K1-1000.0_dp
  else if (alpha >= (gamma0-gamma) .and. alpha <= (gamma0+gamma)) then
    K2 = K1+sqrt(beta**2-(beta/gamma)**2*(alpha-gamma0)**2)+beta0
  else
    K2 = K1-12.0_dp
  end if

  Re_dp = (V*dp_d)/nu

  if (Re_dp <= 5000.0_dp) then
    DeltaK1 = alpha*(1.43_dp*log10(Re_dp))-5.29_dp
  else
    DeltaK1 = 0.0_dp
  end if

  if (r < 1e-8_dp) then
    rc = 1e-8_dp
  else
    rc = r
  end if

  if(alpha <= 12.5_dp) then
    call Afunc(apre,Re,Ap)
    call Afunc(asuc,Re,As)
    call Bfunc(bang,Re,B)

    SPLp = 10.0_dp*log10((dp_d*M**5*L*Dh)/rc**2)+Ap+(K1-3.0_dp)+DeltaK1
    SPLs = 10.0_dp*log10((ds_d*M**5*L*Dh)/rc**2)+As+(K1-3.0_dp)
    SPLa = 10.0_dp*log10((ds_d*M**5*L*Dh)/rc**2)+B+K2

    TBLTE =  10.0_dp*log10(10.0_dp**(SPLp/10.0_dp)+10.0_dp**(SPLs/10.0_dp)&
    +10.0_dp**(SPLa/10.0_dp))

  else
    ! Turbulent Boundary Layer Separation Stall Noise (TBLSS); this is were the airfoil is stalling and stall noise dominates
    ! SPLp = -infinity; 10**(SPLp/10) = 0
    ! SPLs = -infinity; 10**(SPLs/10) = 0

    call Afunc(bang,3.0_dp*Re,A)

    SPLa = 10.0_dp*log10((ds_d*M**5*L*Dl)/rc**2)+A+K2

    TBLTE = 10.0_dp*log10(10.0_dp**(SPLa/10.0_dp))

  end if

end subroutine TBLTEfunc

! Turbulent Boundary Layer Tip Vortex Noise
subroutine TBLTVfunc(f,V,c,r,theta_e,phi_e,atip,conv,tipflat,TBLTV)
  implicit none
  integer, parameter :: dp = kind(0.d0)
  ! in
  real(dp), intent(in) :: f,V,c,r,theta_e,phi_e,atip,conv
  logical, intent(in) :: tipflat
  ! out
  real(dp), intent(out) :: TBLTV
  ! local
  real(dp) :: c0,M,Mc,Mmax,Dh,l,St,rc
  intrinsic log10
  ! constants
  c0 = 343.2_dp ! speed of sound (m/s)
  M = V/c0
  Mc = conv*M
  Mmax = M*(1.0_dp+0.036*atip)

  call Dhfunc(theta_e,phi_e,M,Mc,Dh)

  if (tipflat .eqv. .false.) then
    ! rounded tip
    l = 0.008_dp*c*atip
  else
    ! flat tip
    if (atip <= 2.0_dp .and. atip >= 0.0_dp) then
      l = c*(0.0230_dp+0.0169_dp*atip)
    else
      l = c*(0.0378_dp+0.0095_dp*atip)
    end if
  end if

  St = (f*l)/(V*(1.0_dp+0.036_dp*atip))

  if (r < 1e-8_dp) then
    rc = 1e-8_dp
  else
    rc = r
  end if

  TBLTV =  10.0_dp*log10((M**2*Mmax**3*l**2*Dh)/rc**2)&
  -30.5_dp*(log10(St)+0.3_dp)**2+126.0_dp

end subroutine TBLTVfunc

! Laminar Boundary Layer Vortex Shedding
subroutine LBLVSfunc(f,V,L,c,r,theta_e,phi_e,alpha,nu,conv,trip,LBLVS)
  implicit none
  integer, parameter :: dp = kind(0.d0)
  ! in
  real(dp), intent(in) :: f,V,L,c,r,theta_e,phi_e,alpha,nu,conv
  logical, intent(in) :: trip
  ! out
  real(dp), intent(out) :: LBLVS
  ! local
  real(dp) :: c0,M,Mc,Re,d0,dpr,St,Dh,St1,St_peak,e,G1,Re0,d,G2,G3,rc
  intrinsic log10
  ! constants
  c0 = 343.2_dp ! speed of sound (m/s)
  M = V/c0
  Mc = conv*M
  Re = (V*c)/nu

  if (trip .eqv. .false.) then
    ! UNTRIPPED boundary layer at 0 deg- thickness
    d0 = c*(10.0_dp**(1.6569_dp-0.9045_dp*log10(Re)+0.0596_dp*(log10(Re))**2))
  else
    ! TRIPPED boundary layer at 0 deg- thickness
    d0 = c*(10.0_dp**(1.892_dp-0.9045_dp*log10(Re)+0.0596_dp*(log10(Re))**2))
  end if
  ! boundary layer on pressure side- thickness
  dpr = d0*(10.0_dp**(-0.04175_dp*alpha+0.00106_dp*alpha**2))

  St = (f*dpr)/V

  call Dhfunc(theta_e,phi_e,M,Mc,Dh)

  if (Re <= 1.3e5_dp) then
    St1 = 0.18_dp
  else if (Re <= 4.0e5_dp .and. Re > 1.3e5_dp) then
    St1 = 0.001756_dp*Re**0.3931_dp
  else
    St1 = 0.28_dp
  end if

  St_peak = St1*10.0_dp**(-0.04_dp*alpha)

  e = St/St_peak

  call G1func(e,G1)

  if (alpha <= 3.0_dp) then
    Re0 = 10.0_dp**(0.215_dp*alpha+4.978_dp)
  else
    Re0 = 10.0_dp**(0.12_dp*alpha+5.263_dp)
  end if

  d = Re/Re0

  call G2func(d,G2)
  call G3func(alpha,G3)

  if (r < 1e-8_dp) then
    rc = 1e-8_dp
  else
    rc = r
  end if

  LBLVS = 10.0_dp*log10((dpr*M**5*L*Dh)/rc**2)+G1+G2+G3

end subroutine LBLVSfunc

! Trailing Edge Bluntness Vortex Shedding Noise
subroutine TEBVSfunc(f,V,L,c,h,r,psi,theta_e,phi_e,alpha,nu,conv,trip,TEBVS)
  implicit none
  integer, parameter :: dp = kind(0.d0)
  ! in
  real(dp), intent(in) :: f,V,L,c,h,r,psi,theta_e,phi_e,alpha,nu,conv
  logical, intent(in) :: trip
  ! out
  real(dp), intent(out) :: TEBVS
  ! local
  real(dp) :: c0,M,Mc,Re,d0,d0_d,dpr,dp_d,ds_d,Dh,St,dav,hdav
  real(dp) :: St_peak,StSt_peak,G4,G5,rc
  intrinsic log10
  ! constants
  c0 = 343.2_dp ! speed of sound (m/s)
  M = V/c0
  Mc = conv*M
  Re = (V*c)/nu

  if (trip .eqv. .false.) then
    ! UNTRIPPED boundary layer at 0 deg- thickness, displacement thickness
    d0 = c*(10.0_dp**(1.6569_dp-0.9045_dp*log10(Re)+0.0596_dp*(log10(Re))**2))
    d0_d = c*(10.0_dp**(3.0187_dp-1.5397_dp*log10(Re)+0.1059_dp*(log10(Re))**2))
  else
    ! TRIPPED boundary layer at 0 deg- thickness, displacement thickness
    d0 = c*(10.0_dp**(1.892_dp-0.9045_dp*log10(Re)+0.0596_dp*(log10(Re))**2))
    if (Re <= 0.3e6_dp) then
      d0_d = c*0.0601_dp*Re**(-0.114_dp)
    else
      d0_d = c*(10.0_dp**(3.411_dp-1.5397_dp*log10(Re)+0.1059_dp*(log10(Re))**2))
    end if
  end if

  ! boundary layer on pressure side- thickness, displacement thickness
  dpr = d0*(10.0_dp**(-0.04175_dp*alpha + 0.00106_dp*alpha**2))
  dp_d = d0_d*(10.0_dp**(-0.0432_dp*alpha + 0.00113_dp*alpha**2))

  if (trip .eqv. .false.) then
    ! UNTRIPPED boundary layer on suction side- displacement thickness
    if (alpha <= 7.5_dp .and. alpha >= 0.0_dp) then
      ds_d = d0_d*10.0_dp**(0.0679_dp*alpha)
    else if (alpha <= 12.5_dp .and. alpha > 7.5_dp) then
      ds_d = d0_d*0.0162_dp*10.0_dp**(0.3066_dp*alpha)
    !else if (alpha <= 25_dp .and. alpha > 12.5_dp) then
    else
      ds_d = d0_d* 52.42_dp* 10.0_dp**(0.0258_dp*alpha)
    end if
  else
    ! TRIPPED boundary layer on suction side- displacement thickness
    if (alpha <= 5.0_dp .and. alpha >= 0.0_dp) then
      ds_d = d0_d*10.0_dp**(0.0679_dp*alpha)
    else if (alpha <= 12.5_dp .and. alpha > 5.0_dp) then
      ds_d = d0_d*0.381_dp*10.0_dp**(0.1516_dp*alpha)
    !else if (alpha <= 25.0_dp .and. alpha > 12.5_dp) then
    else
      ds_d = d0_d*14.296_dp*10.0_dp**(0.0258_dp*alpha)
    end if
  end if

  call Dhfunc(theta_e,phi_e,M,Mc,Dh)
  St = (f*h)/V
  dav = (dp_d+ds_d)/2.0_dp

  hdav = h/dav

  if (hdav >= 0.2_dp) then
    St_peak = (0.212_dp-0.0045_dp*psi)/(1.0_dp+0.235_dp*(hdav)**(-1)-0.00132_dp*(hdav)**(-2))
  else
    St_peak = 0.1_dp*(hdav)+0.095_dp-0.00243_dp*psi
  end if

  StSt_peak = St/St_peak

  call G4func(hdav,psi,G4)
  call G5func(hdav,psi,StSt_peak,G5)

  if (r < 1e-8_dp) then
    rc = 1e-8_dp
  else
    rc = r
  end if

  TEBVS = 10.0_dp* log10((h*M**(5.5_dp)*L*Dh)/rc**2)+G4+G5

end subroutine TEBVSfunc

! Computing the overall sound pressure level (OASPL) of a turbine defined below (in dB)
subroutine OASPL(n,r,theta_e,phi_e,rpm,wind,B,rad,c,alpha,SPLoa)
  implicit none
  integer, parameter :: dp = kind(0.d0)
  ! in
  integer, intent(in) :: n
  real(dp), intent(in) :: r,theta_e,phi_e,rpm,wind,B
  real(dp), dimension(n), intent(in) :: rad
  real(dp), dimension(n-1), intent(in) :: c,alpha
  ! out
  real(dp), intent(out) :: SPLoa
  ! local
  integer :: nf,i,j,k
  real(dp) :: pi,nu,conv,omega,atip,psi,refPres,TBLTE,TBLTV,LBLVS,TEBVS
  real(dp) :: TE_t,TV_t,BLVS_t,BVS_t
  real(dp), dimension(n-1) :: wide,rad_mid,vel,L,V,h
  real(dp), dimension(n-1) :: TE_temp,TV_temp,BLVS_temp,BVS_temp
  real(dp), dimension(34) :: f,TE,TV,BLVS,BVS,SPLg,AdB,SPLp
  logical :: trip,tipflat
  intrinsic sqrt
  intrinsic sum
  intrinsic log10
  ! constants
  pi = 3.1415926535897932_dp
  nf = 34

  ! Using untripped or tripped boundary layer specficiation
  trip = .false. ! untripped
  ! trip = .true. ! tripped

  ! Tip specfication
  tipflat = .false. ! round
  ! tipflat = .true. ! flat

  nu = 1.78e-5_dp  ! kinematic viscosity of air (m^2/s)
  conv = 0.8_dp  ! convection factor for speed

  ! Parameters of the wind turbine (f,V,L,c,h,alpha,atip)
  omega = (rpm*2.0_dp*pi)/60.0_dp  ! angular velocity (rad/sec)

  do i = 1,n-1
    wide(i) = rad(i+1)-rad(i) ! length of each radial section (m)
    rad_mid(i) = (wide(i)/2.0_dp)+rad(i) ! midpoint of each radial section (m)
    vel(i) = sqrt((omega*rad_mid(i))**2+wind**2) ! wind speed over the blade (m/s)
    L(i) = wide(i)
    V(i) = vel(i)
  end do

  h(1:n-1) = 0.01_dp*c(1:n-1)  ! trailing edge thickness; 1% of chord length (m)
  atip = alpha(n-1)  ! angle of attack of the tip region (deg)
  psi = 14.0_dp ! solid angle between both airfoil surfaces just upstream of the trailing edge (NACA 0012) (deg)

  ! One-third octave band frequencies (Hz)
  f(1) = 10.0_dp
  f(2) = 12.5_dp
  f(3) = 16.0_dp
  f(4) = 20.0_dp
  f(5) = 25.0_dp
  f(6) = 31.5_dp
  f(7) = 40.0_dp
  f(8) = 50.0_dp
  f(9) = 63.0_dp
  f(10) = 80.0_dp
  f(11) = 100.0_dp
  f(12) = 125.0_dp
  f(13) = 160.0_dp
  f(14) = 200.0_dp
  f(15) = 250.0_dp
  f(16) = 315.0_dp
  f(17) = 400.0_dp
  f(18) = 500.0_dp
  f(19) = 630.0_dp
  f(20) = 800.0_dp
  f(21) = 1000.0_dp
  f(22) = 1250.0_dp
  f(23) = 1600.0_dp
  f(24) = 2000.0_dp
  f(25) = 2500.0_dp
  f(26) = 3150.0_dp
  f(27) = 4000.0_dp
  f(28) = 5000.0_dp
  f(29) = 6300.0_dp
  f(30) = 8000.0_dp
  f(31) = 10000.0_dp
  f(32) = 12500.0_dp
  f(33) = 16000.0_dp
  f(34) = 20000.0_dp

  refPres = 2e-5_dp

  ! Calculating sound pressure (Pa) for each noise source at each frequency and radial position
  do j=1,nf
    do k=1,n-1
      call TBLTEfunc(f(j),V(k),L(k),c(k),r,theta_e,phi_e,alpha(k),nu,conv,&
      trip,TBLTE)
      call TBLTVfunc(f(j),V(k),c(k),r,theta_e,phi_e,atip,conv,tipflat,TBLTV)
      call LBLVSfunc(f(j),V(k),L(k),c(k),r,theta_e,phi_e,alpha(k),nu,conv,&
      trip,LBLVS)
      call TEBVSfunc(f(j),V(k),L(k),c(k),h(k),r,psi,theta_e,phi_e,alpha(k),&
      nu,conv,trip,TEBVS)

      TE_temp(k) = refPres*10.0_dp**(TBLTE/20.0_dp)
      TV_temp(k) = refPres*10.0_dp**(TBLTV/20.0_dp)
      BLVS_temp(k) = refPres*10.0_dp**(LBLVS/20.0_dp)
      BVS_temp(k) = refPres*10.0_dp**(TEBVS/20.0_dp)
    end do

    ! Adding incoherent noises
    TE_t = sum(TE_temp)
    TV_t = sum(TV_temp)
    BLVS_t = sum(BLVS_temp)
    BVS_t = sum(BVS_temp)

    ! Converting to sound pressure levels (SPL) (dB)
    TE(j) = 10.0_dp*log10((TE_t/refPres)**2)
    TV(j) = 10.0_dp*log10((TV_t/refPres)**2)
    BLVS(j) = 10.0_dp*log10((BLVS_t/refPres)**2)
    BVS(j) = 10.0_dp*log10((BVS_t/refPres)**2)

    ! multiplied by number of blades (B)
    SPLg(j) = 10.0_dp*log10(B*(10.0_dp**(TE(j)/10.0_dp)+10.0_dp**(TV(j)/&
    10.0_dp)+10.0_dp**(BLVS(j)/10.0_dp)+10.0_dp**(BVS(j)/10.0_dp)))
  end do

  ! A-weighting curve (dBA) for sound perception correction
  AdB(1) = -70.4_dp
  AdB(2) = -63.4_dp
  AdB(3) = -56.7_dp
  AdB(4) = -50.5_dp
  AdB(5) = -44.7_dp
  AdB(6) = -39.4_dp
  AdB(7) = -34.6_dp
  AdB(8) = -30.2_dp
  AdB(9) = -26.2_dp
  AdB(10) = -22.5_dp
  AdB(11) = -19.1_dp
  AdB(12) = -16.1_dp
  AdB(13) = -13.4_dp
  AdB(14) = -10.9_dp
  AdB(15) = -8.6_dp
  AdB(16) = -6.6_dp
  AdB(17) = -4.8_dp
  AdB(18) = -3.2_dp
  AdB(19) = -1.9_dp
  AdB(20) = -0.8_dp
  AdB(21) = 0.0_dp
  AdB(22) = 0.6_dp
  AdB(23) = 1.0_dp
  AdB(24) = 1.2_dp
  AdB(25) = 1.3_dp
  AdB(26) = 1.2_dp
  AdB(27) = 1.0_dp
  AdB(28) = 0.5_dp
  AdB(29) = -0.1_dp
  AdB(30) = -1.1_dp
  AdB(31) = -2.5_dp
  AdB(32) = -4.3_dp
  AdB(33) = -6.6_dp
  AdB(34) = -9.3_dp

  SPLg(1:nf) = SPLg(1:nf) + AdB(1:nf)

  ! Converting to sound pressure (Pa) for incoherent noise addition
  SPLp(1:nf) = refPres*10.0_dp**(SPLg(1:nf)/20.0_dp)

  ! Converting back to SPL (dB)
  SPLoa = 10.0_dp*log10((sum(SPLp)/refPres)**2)

end subroutine OASPL

! Placing a turbine in a specified location and finding the OASPL of the turbine with reference to an observer
subroutine turbinepos(nturb,nnrel,nobs,x,y,obs,wind,rpm,windvel,B,h,rad,c,alpha,corr,turbinenoise)
  implicit none
  integer, parameter :: dp = kind(0.d0)
  ! in
  integer, intent(in) :: nturb,nnrel,nobs
  real(dp), intent(in) :: wind,B,h,corr
  real(dp), dimension(nturb), intent(in) :: x,y,rpm,windvel
  real(dp), dimension(nobs), intent(in) :: obs
  real(dp), dimension(nnrel), intent(in) :: rad
  real(dp), dimension(nnrel-1), intent(in) :: c,alpha
  ! out
  real(dp), intent(out) :: turbinenoise
  ! local
  integer :: i,j,k
  real(dp) :: pi,windrad,xi,yi,rxy,ang,turbine_noise,SPLoa,turbnoise
  real(dp), dimension(nturb) :: r,theta_e,phi_e,t,theta_er,phi_er
  real(dp), dimension(nobs) :: obsi
  intrinsic sqrt
  intrinsic sin
  intrinsic cos
  intrinsic atan2
  intrinsic abs
  ! constants
  pi = 3.1415926535897932_dp

  windrad = wind*pi/180.0_dp

  do i = 1,nturb
    r(i) = sqrt((obs(1)-x(i))**2+(obs(2)-y(i))**2+(obs(3)-h)**2)

    obsi(1) = obs(1)-x(i)
    obsi(2) = obs(2)-y(i)
    obsi(3) = obs(3)
    xi = 0.0_dp
    yi = 0.0_dp

    ! Adjusting the coordinates for the wind direction
    rxy = sqrt((obsi(1))**2+(obsi(2))**2)
    ang = atan2(obsi(2),obsi(1))+windrad

    obsi(1) = rxy*cos(ang)
    obsi(2) = rxy*sin(ang)

    phi_e(i) = abs((atan2((obsi(2)-yi),(obsi(1)-xi)))*180.0_dp/pi)
    theta_e(i) = (atan2(abs(h-obsi(3)),abs(obsi(2)-yi)))*180.0_dp/pi
    if (phi_e(i) < 5.0_dp) then
      phi_e(i) = (1.0_dp/10.0_dp)*phi_e(i)**2 + 2.5_dp
    end if
    if (phi_e(i) > 175.0_dp) then
      phi_er(i) = 180.0_dp - phi_e(i)
      phi_er(i) = (1.0_dp/10.0_dp)*phi_er(i)**2 + 2.5_dp
      phi_e(i) = 180.0_dp - phi_er(i)
    end if

    if (theta_e(i) > 85.0_dp) then
      theta_er(i) = 90.0_dp - theta_e(i)
      theta_er(i) = (1.0_dp/10.0_dp)*theta_er(i)**2 + 2.5_dp
      theta_e(i) = 90.0_dp - theta_er(i)
    end if

    ! if (phi_e(i) == 0.0_dp .or. phi_e(i) == 180.0_dp) then
    !   r(i) = r(i) - 0.2 ! directivity adjustment based on work by Luis Vargas (Wind Turbine Noise Prediction)
    ! end if
  end do

  ! Calculating the OASPL of each of the turbines in a for loop of "n" turbines (SPL function)
  do j = 1,nturb
    call OASPL(nnrel,r(j),theta_e(j),phi_e(j),rpm(j),windvel(j),B,rad,c,alpha,SPLoa)
    t(j) = SPLoa
  end do

  ! Calculating the preliminary turbine_noise (10**(t1/10) + 10**(t2/10) + ...)
  turbine_noise = 0.0_dp
  do k = 1,nturb
    turbine_noise = turbine_noise+10.0_dp**(t(k)/10.0_dp)
  end do

  turbnoise = 10.0_dp*log10(turbine_noise)

  turbinenoise = turbnoise+corr ! correction based on Rosiere valdiation (243.84 m, 800 ft should be 47 dB)

end subroutine turbinepos


! To build for Python interface:
! f2py -c  --opt=-O2 -m _bpmacoustic BPM_Acoustic_Model.f90
! python C:\Python27\Scripts\f2py.py -c --opt=-O2 --compiler=mingw32 --fcompiler=gfortran -m _bpmacoustic BPM_Acoustic_Model.f90
