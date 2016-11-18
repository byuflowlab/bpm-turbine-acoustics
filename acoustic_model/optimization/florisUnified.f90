! Flow field calculations have been intentionally left out to save development time.
! The flow field can be calculated using the pure python version of floris 

! This implementation is fully smooth and differentiable with the exception of a 
! discontinuity at the hub of each turbine. The discontinuity only presents issues if
! turbines are place within 1E-15 * rotor diameter of one another, which is extremely 
! unlikely during optimization if the user does not explicitly place them there.

    
subroutine Hermite_Spline(x, x0, x1, y0, dy0, y1, dy1, y)
    !    This function produces the y and dy values for a hermite cubic spline
    !    interpolating between two end points with known slopes
    !
    !    :param x: x position of output y
    !    :param x0: x position of upwind endpoint of spline
    !    :param x1: x position of downwind endpoint of spline
    !    :param y0: y position of upwind endpoint of spline
    !    :param dy0: slope at upwind endpoint of spline
    !    :param y1: y position of downwind endpoint of spline
    !    :param dy1: slope at downwind endpoint of spline
    !
    !    :return: y: y value of spline at location x
    
    implicit none
        
    ! define precision to be the standard for a double precision ! on local system
    integer, parameter :: dp = kind(0.d0)
    
    ! in
    real(dp), intent(in) :: x, x0, x1, y0, dy0, y1, dy1
    
    ! out
    real(dp), intent(out) :: y !, dy_dx
    
    ! local
    real(dp) :: c3, c2, c1, c0

    ! initialize coefficients for parametric cubic spline
    c3 = (2.0_dp*(y1))/(x0**3 - 3.0_dp*x0**2*x1 + 3.0_dp*x0*x1**2 - x1**3) - &
         (2.0_dp*(y0))/(x0**3 - 3.0_dp*x0**2*x1 + 3.0_dp*x0*x1**2 - x1**3) + &
         (dy0)/(x0**2 - 2.0_dp*x0*x1 + x1**2) + &
         (dy1)/(x0**2 - 2.0_dp*x0*x1 + x1**2)
         
    c2 = (3.0_dp*(y0)*(x0 + x1))/(x0**3 - 3.0_dp*x0**2*x1 + 3.0_dp*x0*x1**2 - x1**3) - &
         ((dy1)*(2.0_dp*x0 + x1))/(x0**2 - 2.0_dp*x0*x1 + x1**2) - ((dy0)*(x0 + &
         2.0_dp*x1))/(x0**2 - 2.0_dp*x0*x1 + x1**2) - (3.0_dp*(y1)*(x0 + x1))/(x0**3 - &
         3.0_dp*x0**2*x1 + 3.0_dp*x0*x1**2 - x1**3)
         
    c1 = ((dy0)*(x1**2 + 2.0_dp*x0*x1))/(x0**2 - 2.0_dp*x0*x1 + x1**2) + ((dy1)*(x0**2 + &
         2.0_dp*x1*x0))/(x0**2 - 2.0_dp*x0*x1 + x1**2) - (6.0_dp*x0*x1*(y0))/(x0**3 - &
         3.0_dp*x0**2*x1 + 3.0_dp*x0*x1**2 - x1**3) + (6.0_dp*x0*x1*(y1))/(x0**3 - &
         3.0_dp*x0**2*x1 + 3.0_dp*x0*x1**2 - x1**3)
         
    c0 = ((y0)*(- x1**3 + 3.0_dp*x0*x1**2))/(x0**3 - 3.0_dp*x0**2*x1 + 3.0_dp*x0*x1**2 - &
         x1**3) - ((y1)*(- x0**3 + 3.0_dp*x1*x0**2))/(x0**3 - 3.0_dp*x0**2*x1 + &
         3.0_dp*x0*x1**2 - x1**3) - (x0*x1**2*(dy0))/(x0**2 - 2.0_dp*x0*x1 + x1**2) - &
         (x0**2*x1*(dy1))/(x0**2 - 2.0_dp*x0*x1 + x1**2)
!    print *, 'c3 = ', c3
!    print *, 'c2 = ', c2
!    print *, 'c1 = ', c1
!    print *, 'c0 = ', c0
    ! Solve for y and dy values at the given point
    y = c3*x**3 + c2*x**2 + c1*x + c0
    !dy_dx = c3*3*x**2 + c2*2*x + c1

end subroutine Hermite_Spline


subroutine calcOverlapAreas(nTurbines, turbineX, turbineY, rotorDiameter, wakeDiameters, &
                            wakeCenters, wakeOverlapTRel_mat)
!    calculate overlap of rotors and wake zones (wake zone location defined by wake 
!    center and wake diameter)
!   turbineX,turbineY is x,y-location of center of rotor
!
!    wakeOverlap(TURBI,TURB,ZONEI) = overlap area of zone ZONEI of wake of turbine 
!     TURB with rotor of downstream turbine
!    TURBI

    implicit none
        
    ! define precision to be the standard for a double precision ! on local system
    integer, parameter :: dp = kind(0.d0)
    
    ! in
    integer, intent(in) :: nTurbines
    real(dp), dimension(nTurbines), intent(in) :: turbineX, turbineY, rotorDiameter
    real(dp), dimension(nTurbines, nTurbines, 3), intent(in) :: wakeDiameters
    real(dp), dimension(nTurbines, nTurbines), intent(in) :: wakeCenters
    
    ! out    
    real(dp), dimension(nTurbines, nTurbines, 3), intent(out) :: wakeOverlapTRel_mat
    
    ! local
    integer :: turb, turbI, zone
    real(dp), parameter :: pi = 3.141592653589793_dp, tol = 0.000001_dp
    real(dp) :: OVdYd, OVr, OVRR, OVL, OVz
    real(dp), dimension(nTurbines, nTurbines, 3) :: wakeOverlap
        
    wakeOverlapTRel_mat = 0.0_dp
    wakeOverlap = 0.0_dp
    
    do turb = 1, nTurbines
        do turbI = 1, nTurbines
            if (turbineX(turbI) > turbineX(turb)) then
                OVdYd = wakeCenters(turbI, turb)-turbineY(turbI)    ! distance between wake center and rotor center
                OVr = rotorDiameter(turbI)/2                        ! rotor diameter
                do zone = 1, 3
                    OVRR = wakeDiameters(turbI, turb, zone)/2.0_dp        ! wake diameter
                    OVdYd = abs(OVdYd)
                    if (OVdYd >= 0.0_dp + tol) then
                        ! calculate the distance from the wake center to the vertical line between
                        ! the two circle intersection points
                        OVL = (-OVr*OVr+OVRR*OVRR+OVdYd*OVdYd)/(2.0_dp*OVdYd)
                    else
                        OVL = 0.0_dp
                    end if

                    OVz = OVRR*OVRR-OVL*OVL

                    ! Finish calculating the distance from the intersection line to the outer edge of the wake zone
                    if (OVz > 0.0_dp + tol) then
                        OVz = sqrt(OVz)
                    else
                        OVz = 0.0_dp
                    end if

                    if (OVdYd < (OVr+OVRR)) then ! if the rotor overlaps the wake zone

                        if (OVL < OVRR .and. (OVdYd-OVL) < OVr) then
                            wakeOverlap(turbI, turb, zone) = OVRR*OVRR*dacos(OVL/OVRR) + OVr*OVr*dacos((OVdYd-OVL)/OVr) - OVdYd*OVz
                        else if (OVRR > OVr) then
                            wakeOverlap(turbI, turb, zone) = pi*OVr*OVr
                        else
                            wakeOverlap(turbI, turb, zone) = pi*OVRR*OVRR
                        end if
                    else
                        wakeOverlap(turbI, turb, zone) = 0.0_dp
                    end if
                    
                end do
                
            end if
            
        end do
        
    end do


    do turb = 1, nTurbines
    
        do turbI = 1, nTurbines
    
            wakeOverlap(turbI, turb, 3) = wakeOverlap(turbI, turb, 3)-wakeOverlap(turbI, turb, 2)
            wakeOverlap(turbI, turb, 2) = wakeOverlap(turbI, turb, 2)-wakeOverlap(turbI, turb, 1)
    
        end do
    
    end do
    
    wakeOverlapTRel_mat = wakeOverlap

    do turbI = 1, nTurbines
            wakeOverlapTRel_mat(turbI, :, :) = wakeOverlapTRel_mat(turbI, :, &
                                                         :)/((pi*rotorDiameter(turbI) &
                                                       *rotorDiameter(turbI))/4.0_dp)
    end do
    
    ! do turbI = 1, nTurbines
!         do turb = 1, nTurbines
!             do zone = 1, 3
!                 print *, "wakeOverlapTRel_mat[", turbI, ", ", turb, ", ", zone, "] = ", wakeOverlapTRel_mat(turbI, turb, zone)
!             end do
!         end do
!     end do
        
   
                                    
end subroutine calcOverlapAreas


subroutine CTtoAxialInd(CT, nTurbines, axial_induction)
    
    implicit none
    
    ! define precision to be the standard for a double precision ! on local system
    integer, parameter :: dp = kind(0.d0)

    ! in
    integer, intent(in) :: nTurbines
    real(dp), dimension(nTurbines), intent(in) :: CT

    ! local
    integer :: i

    ! out
    real(dp), dimension(nTurbines), intent(out) :: axial_induction

    axial_induction = 0.0_dp

    ! execute
    do i = 1, nTurbines
        if (CT(i) > 0.96) then  ! Glauert condition
            axial_induction(i) = 0.143_dp + sqrt(0.0203_dp-0.6427_dp*(0.889_dp - CT(i)))
        else
            axial_induction(i) = 0.5_dp*(1.0_dp-sqrt(1.0_dp-CT(i)))
        end if
    end do
    
end subroutine CTtoAxialInd
    

subroutine floris_unified(nTurbines, turbineXw, turbineYw, yaw_deg, rotorDiameter, Vinf, &
                          & Ct, a_in, ke_in, kd, me, initialWakeDisplacement, bd, MU, &
                          & aU, bU, initialWakeAngle, cos_spread, keCorrCT, Region2CT, &
                          & keCorrArray, useWakeAngle, adjustInitialWakeDiamToYaw, & 
                          & axialIndProvided, useaUbU, velocitiesTurbines, &
                          & wakeCentersYT_vec, wakeDiametersT_vec, wakeOverlapTRel_vec)
    
    ! independent variables: yaw_deg Ct turbineXw turbineYw rotorDiameter a_in    
    ! dependent variables: velocitiesTurbines
    
    implicit none
    
    ! define precision to be the standard for a double precision ! on local system
    integer, parameter :: dp = kind(0.d0)
    
    ! in
    integer, intent(in) :: nTurbines
    real(dp), intent(in) :: kd, initialWakeDisplacement, initialWakeAngle, ke_in
    real(dp), intent(in) :: keCorrCT, Region2CT, bd, cos_spread, Vinf, keCorrArray
    real(dp), dimension(nTurbines), intent(in) :: yaw_deg, Ct, a_in, turbineXw, turbineYw
    real(dp), dimension(nTurbines), intent(in) :: rotorDiameter
    real(dp), dimension(3), intent(in) :: me, MU
    real(dp), intent(in) :: aU, bU
    logical, intent(in) :: useWakeAngle, adjustInitialWakeDiamToYaw, axialIndProvided, &
                           & useaUbU
                           
    ! local (General)
    real(dp), dimension(nTurbines) :: ke, yaw
    real(dp) :: deltax
    Integer :: turb, turbI, zone   
    real(dp), parameter :: pi = 3.141592653589793_dp
    
    ! local (Wake centers and diameters)
    real(dp) :: spline_bound ! in rotor diameters    
    real(dp) :: wakeAngleInit, zeroloc
    real(dp) :: factor, displacement, x, x1, x2, y1, y2, dy1, dy2
    real(dp) :: wakeDiameter0
    real(dp), dimension(nTurbines, nTurbines, 3) :: wakeDiametersT_mat
    real(dp), dimension(nTurbines, nTurbines) :: wakeCentersYT_mat
    
    ! local (Wake overlap)
    real(dp) :: rmax
    real(dp), dimension(nTurbines, nTurbines, 3) :: wakeOverlapTRel_mat

    ! local (Velocity)
    real(dp), dimension(nTurbines) :: a, keArray
    real(dp), dimension(3) :: mmU
    real(dp) :: s, cosFac, wakeEffCoeff, wakeEffCoeffPerZone
    
    ! out
    real(dp), dimension(nTurbines), intent(out) :: velocitiesTurbines
    real(dp), dimension(nTurbines*nTurbines), intent(out) :: wakeCentersYT_vec
    real(dp), dimension(3*nTurbines*nTurbines), intent(out) :: wakeDiametersT_vec    
    real(dp), dimension(3*nTurbines*nTurbines), intent(out) :: wakeOverlapTRel_vec
    
    intrinsic cos
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!! Wake Centers and Diameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    spline_bound = 1.0_dp
       
    yaw = yaw_deg*pi/180.0_dp
        
    ! calculate y-locations of wake centers in wind ref. frame
    wakeCentersYT_mat = 0.0_dp
    
    do turb = 1, nTurbines
        wakeAngleInit = 0.5_dp*sin(yaw(turb))*Ct(turb)
        
        if (useWakeAngle) then
            wakeAngleInit = wakeAngleInit + initialWakeAngle*pi/180.0_dp
        end if
        
        do turbI = 1, nTurbines            
            
            if (turbineXw(turb) < turbineXw(turbI)) then
                deltax = turbineXw(turbI) - turbineXw(turb)
                factor = (2.0_dp*kd*deltax/rotorDiameter(turb)) + 1.0_dp
                wakeCentersYT_mat(turbI, turb) = turbineYw(turb)
                
                displacement = wakeAngleInit*(wakeAngleInit* &
                                                 & wakeAngleInit + 15.0_dp*factor*factor* &
                                                 factor*factor)/((30.0_dp*kd/ & 
                                                 rotorDiameter(turb))*(factor*factor* &
                                                 & factor*factor*factor))
                                                 
                displacement = displacement - &
                                                 & wakeAngleInit*(wakeAngleInit* &
                                                 & wakeAngleInit + 15.0_dp)/(30.0_dp*kd/ &
                                                 rotorDiameter(turb))

                wakeCentersYT_mat(turbI, turb) = wakeCentersYT_mat(turbI, turb)+ &
                                                 & initialWakeDisplacement + displacement
  
                if (useWakeAngle .eqv. .false.) then
                    wakeCentersYT_mat(turbI, turb) = wakeCentersYT_mat(turbI, turb) + bd*(deltax)
                end if
                
            end if

        end do
        
    end do
    
    !adjust k_e to C_T, adjusted to yaw
    ke = ke_in + keCorrCT*(Ct-Region2CT)
    
    wakeDiametersT_mat = 0.0_dp

    do turb = 1, nTurbines
        
        if (adjustInitialWakeDiamToYaw) then
            wakeDiameter0 = rotorDiameter(turb)*cos(yaw(turb))
        else
            wakeDiameter0 = rotorDiameter(turb)        
        end if
        
        do turbI = 1, nTurbines
        
            ! turbine separation
            deltax = turbineXw(turbI) - turbineXw(turb)            
            
            ! x position of interest
            x = turbineXw(turbI)                         
              
            zone = 1
            ! define centerpoint of spline
            zeroloc = turbineXw(turb) - wakeDiameter0/(2.0_dp*ke(turb)*me(zone))
            
            if (zeroloc + spline_bound*rotorDiameter(turb) < turbineXw(turbI)) then ! check this
                wakeDiametersT_mat(turbI, turb, zone) = 0.0_dp
            
            else if (zeroloc - spline_bound*rotorDiameter(turb) < turbineXw(turbI)) then !check this
                               
                !!!!!!!!!!!!!!!!!!!!!! calculate spline values !!!!!!!!!!!!!!!!!!!!!!!!!!
                
                ! position of upwind point
                x1 = zeroloc - spline_bound*rotorDiameter(turb)
                
                ! diameter of upwind point
                y1 = wakeDiameter0+2.0_dp*ke(turb)*me(zone)*(x1 - turbineXw(turb))
                                
                ! slope at upwind point
                dy1 = 2.0_dp*ke(turb)*me(zone)
                
                ! position of downwind point
                x2 = zeroloc+spline_bound*rotorDiameter(turb)             

                ! diameter at downwind point
                y2 = 0.0_dp
                
                ! slope at downwind point
                dy2 = 0.0_dp
                
                ! solve for the wake zone diameter and its derivative w.r.t. the downwind
                ! location at the point of interest
                call Hermite_Spline(x, x1, x2, y1, dy1, y2, dy2, wakeDiametersT_mat(turbI, turb, zone))
            
            else if (turbineXw(turb) < turbineXw(turbI)) then
                wakeDiametersT_mat(turbI, turb, zone) = wakeDiameter0+2.0_dp*ke(turb)*me(zone)*deltax            
            end if
            
                        
            if (turbineXw(turb) < turbineXw(turbI)) then
                zone = 2
                wakeDiametersT_mat(turbI, turb, zone) = wakeDiameter0 + 2.0_dp*ke(turb)*me(zone)*deltax                   
                zone = 3
                wakeDiametersT_mat(turbI, turb, zone) = wakeDiameter0+2.0_dp*ke(turb)*me(zone)*deltax
            end if      
            
        end do
        
    end do
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Wake Overlap !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! calculate relative overlap
    call calcOverlapAreas(nTurbines, turbineXw, turbineYw, rotorDiameter, &
                          & wakeDiametersT_mat, wakeCentersYT_mat, wakeOverlapTRel_mat)
!     
!     ! calculate cosine factor TODO: put this into the same loop as velocity
!     do turbI = 1, nTurbines
!         do turb = 1, nTurbines
!             if (turbineXw(turb) < turbineXw(turbI)) then
!                 do zone = 1, 3
!                     rmax = cos_spread*0.5_dp*(wakeDiametersT_mat(turbI, turb, 3) + rotorDiameter(turbI))
!                     cosFac_mat(turbI, turb, zone) = 0.5_dp*(1.0_dp + cos(pi*dabs( &
!                         wakeCentersYT_mat(turbI, turb)-turbineYw(turbI))/rmax))
!                     !cosFac_mat(turbI, turb, zone) = 1.0_dp
!             
!                 end do
!             else
!                 cosFac_mat(turbI, turb, :) = 1.0_dp               
!             end if            
!         end do
!     end do
        
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Velocity !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
    ! convert yaw from degrees to radians
    yaw = yaw_deg*pi/180.0_dp

    if (axialIndProvided) then
        a = a_in
    else
        call CTtoAxialInd(Ct, nTurbines, a)
    end if
    
    ! adjust ke to Ct as adjusted to yaw
    ke = ke_in + keCorrCT*(Ct-Region2CT)

    do turb = 1, nTurbines
        s = sum(wakeOverlapTRel_mat(turb, :, 1) + wakeOverlapTRel_mat(turb, :, 2))
        keArray(turb) = ke(turb)*(1+s*keCorrArray) 
    end do
    
    ! find effective wind speeds at downstream turbines, then predict the power of 
    ! downstream turbine    
    velocitiesTurbines = Vinf
    do turbI = 1, nTurbines
        wakeEffCoeff = 0.0_dp
        
        ! find overlap-area weighted effect of each wake zone
        do turb = 1, nTurbines
            wakeEffCoeffPerZone = 0.0_dp
            deltax = turbineXw(turbI) - turbineXw(turb)
            
            if (useaUbU) then
                mmU = MU/cos(aU*pi/180.0_dp + bU*yaw(turb))
            end if
            
            if (deltax > 0 .and. turbI /= turb) then
                do zone = 1, 3
                
                    rmax = cos_spread*0.5_dp*(wakeDiametersT_mat(turbI, turb, 3) + rotorDiameter(turbI))
                    cosFac = 0.5_dp*(1.0_dp + cos(pi*dabs(wakeCentersYT_mat(turbI, turb) &
                                     & - turbineYw(turbI))/rmax))
                                
                    if (useaUbU) then
                        wakeEffCoeffPerZone = wakeEffCoeffPerZone + &
                        (((cosFac*rotorDiameter(turb))/(rotorDiameter(turb)+2.0_dp*keArray(turb) &
                        *mmU(zone)*deltax))**2)*wakeOverlapTRel_mat(turbI, turb, zone)   
                    else
                        wakeEffCoeffPerZone = wakeEffCoeffPerZone + &
                        (((cosFac*rotorDiameter(turb))/(rotorDiameter(turb)+2.0_dp*keArray(turb) &
                        *MU(zone)*deltax))**2)*wakeOverlapTRel_mat(turbI, turb, zone)   
                    end if                     
                            
                end do
                wakeEffCoeff = wakeEffCoeff + (a(turb)*wakeEffCoeffPerZone)**2
            end if
        end do
        wakeEffCoeff = 1.0_dp - 2.0_dp*sqrt(wakeEffCoeff)
        
        ! multiply the inflow speed with the wake coefficients to find effective wind 
        ! speed at turbine
        velocitiesTurbines(turbI) = velocitiesTurbines(turbI)*wakeEffCoeff
    end do
    
    ! pack desired matrices into vectors for output
    do turbI = 1, nTurbines
        ! wake centers
        wakeCentersYT_vec(nTurbines*(turbI-1)+1:nTurbines*(turbI-1)+nTurbines) &
                                     = wakeCentersYT_mat(turbI, :)
                                     
        ! wake diameters
        wakeDiametersT_vec(3*nTurbines*(turbI-1)+1:3*nTurbines*(turbI-1)+nTurbines) &
                                 = wakeDiametersT_mat(turbI, :, 1)
        wakeDiametersT_vec(3*nTurbines*(turbI-1)+nTurbines+1:3*nTurbines*(turbI-1) &
                                   +2*nTurbines) = wakeDiametersT_mat(turbI, :, 2)
        wakeDiametersT_vec(3*nTurbines*(turbI-1)+2*nTurbines+1:nTurbines*(turbI-1) &
                                   +3*nTurbines) = wakeDiametersT_mat(turbI, :, 3) 
        
        ! relative wake overlap
        wakeOverlapTRel_vec(3*nTurbines*(turbI-1)+1:3*nTurbines*(turbI-1)+nTurbines) &
                             = wakeOverlapTRel_mat(turbI, :, 1)
        wakeOverlapTRel_vec(3*nTurbines*(turbI-1)+nTurbines+1:3*nTurbines*(turbI-1) &
                               +2*nTurbines) = wakeOverlapTRel_mat(turbI, :, 2)
        wakeOverlapTRel_vec(3*nTurbines*(turbI-1)+2*nTurbines+1:3*nTurbines*(turbI-1) &
                               +3*nTurbines) = wakeOverlapTRel_mat(turbI, :, 3)
        
       
    end do
    
end subroutine floris_unified




!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.11 (r5902M) - 15 Dec 2015 09:00
!
!  Differentiation of floris_unified in reverse (adjoint) mode:
!   gradient     of useful results: velocitiesturbines
!   with respect to varying inputs: rotordiameter turbinexw yaw_deg
!                velocitiesturbines turbineyw ct a_in
!   RW status of diff variables: rotordiameter:out turbinexw:out
!                yaw_deg:out velocitiesturbines:in-zero turbineyw:out
!                ct:out a_in:out
SUBROUTINE FLORIS_UNIFIED_BV(nturbines, turbinexw, turbinexwb, turbineyw&
& , turbineywb, yaw_deg, yaw_degb, rotordiameter, rotordiameterb, vinf, &
& ct, ctb, a_in, a_inb, ke_in, kd, me, initialwakedisplacement, bd, mu, &
& au, bu, initialwakeangle, cos_spread, kecorrct, region2ct, kecorrarray&
& , usewakeangle, adjustinitialwakediamtoyaw, axialindprovided, useaubu&
& , velocitiesturbinesb, nbdirs)
  
!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
! in
  INTEGER, INTENT(IN) :: nturbines
  REAL(dp), INTENT(IN) :: kd, initialwakedisplacement, initialwakeangle&
& , ke_in
  REAL(dp), INTENT(IN) :: kecorrct, region2ct, bd, cos_spread, vinf, &
& kecorrarray
  REAL(dp), DIMENSION(nturbines), INTENT(IN) :: yaw_deg, ct, a_in, &
& turbinexw, turbineyw
  REAL(dp), DIMENSION(nbdirs, nturbines), intent(out) :: yaw_degb, ctb, a_inb, &
& turbinexwb, turbineywb
  REAL(dp), DIMENSION(nturbines), INTENT(IN) :: rotordiameter
  REAL(dp), DIMENSION(nbdirs, nturbines), intent(out) :: rotordiameterb
  REAL(dp), DIMENSION(3), INTENT(IN) :: me, mu
  REAL(dp), INTENT(IN) :: au, bu
  LOGICAL, INTENT(IN) :: usewakeangle, adjustinitialwakediamtoyaw, &
& axialindprovided, useaubu
! local (General)
  REAL(dp), DIMENSION(nturbines) :: ke, yaw
  REAL(dp), DIMENSION(nbdirs, nturbines) :: keb, yawb
  REAL(dp) :: deltax
  REAL(dp), DIMENSION(nbdirs) :: deltaxb
  INTEGER :: turb, turbi, zone
  REAL(dp), PARAMETER :: pi=3.141592653589793_dp
! local (Wake centers and diameters)
! in rotor diameters    
  REAL(dp) :: spline_bound
  REAL(dp) :: wakeangleinit, zeroloc
  REAL(dp), DIMENSION(nbdirs) :: wakeangleinitb, zerolocb
  REAL(dp) :: factor, displacement, x, x1, x2, y1, y2, dy1, dy2
  REAL(dp), DIMENSION(nbdirs) :: factorb, displacementb, xb, x1b, x2b&
& , y1b, dy1b
  REAL(dp) :: wakediameter0
  REAL(dp), DIMENSION(nbdirs) :: wakediameter0b
  REAL(dp), DIMENSION(nturbines, nturbines, 3) :: wakediameterst_mat
  REAL(dp), DIMENSION(nbdirs, nturbines, nturbines, 3) :: &
& wakediameterst_matb
  REAL(dp), DIMENSION(nturbines, nturbines) :: wakecentersyt_mat
  REAL(dp), DIMENSION(nbdirs, nturbines, nturbines) :: &
& wakecentersyt_matb
! local (Wake overlap)
  REAL(dp) :: rmax
  REAL(dp), DIMENSION(nbdirs) :: rmaxb
  REAL(dp), DIMENSION(nturbines, nturbines, 3) :: wakeoverlaptrel_mat
  REAL(dp), DIMENSION(nbdirs, nturbines, nturbines, 3) :: &
& wakeoverlaptrel_matb
! local (Velocity)
  REAL(dp), DIMENSION(nturbines) :: a, kearray
  REAL(dp), DIMENSION(nbdirs, nturbines) :: ab, kearrayb
  REAL(dp), DIMENSION(3) :: mmu
  REAL(dp), DIMENSION(nbdirs, 3) :: mmub
  REAL(dp) :: s, cosfac, wakeeffcoeff, wakeeffcoeffperzone
  REAL(dp), DIMENSION(nbdirs) :: sb, cosfacb, wakeeffcoeffb, &
& wakeeffcoeffperzoneb
! out
  REAL(dp), DIMENSION(nturbines) :: velocitiesturbines
  REAL(dp), DIMENSION(nbdirs, nturbines) :: velocitiesturbinesb
  REAL(dp), DIMENSION(nturbines*nturbines) :: &
& wakecentersyt_vec
  REAL(dp), DIMENSION(3*nturbines*nturbines) :: &
& wakediameterst_vec
  REAL(dp), DIMENSION(3*nturbines*nturbines) :: &
& wakeoverlaptrel_vec
  INTRINSIC COS
  INTRINSIC KIND
  INTRINSIC SIN
  INTRINSIC SUM
  INTRINSIC DABS
  INTRINSIC SQRT
  INTEGER :: nd
  INTEGER :: branch
  INTEGER :: nbdirs
  REAL(dp) :: temp3
  REAL(dp) :: temp2
  REAL(dp) :: temp1
  REAL(dp) :: temp0
  REAL(dp) :: tempb9(nbdirs)
  REAL(dp) :: tempb8(nbdirs)
  REAL :: tempb7(nbdirs)
  REAL(dp) :: tempb6(nbdirs)
  REAL(dp) :: tempb5(nbdirs)
  REAL(dp) :: tempb4(nbdirs)
  REAL(dp) :: tempb3(nbdirs)
  REAL(dp) :: tempb2(nbdirs)
  REAL(dp) :: tempb1(nbdirs)
  REAL(dp) :: tempb0(nbdirs)
  REAL(dp) :: tempb13(nbdirs)
  REAL(dp) :: tempb12(nbdirs)
  REAL(dp) :: tempb11(nbdirs)
  REAL(dp) :: tempb10(nbdirs)
  DOUBLE PRECISION :: dabs0b(nbdirs)
  REAL(dp) :: temp11
  REAL(dp) :: temp10
  REAL(dp) :: tempb(nbdirs)
  DOUBLE PRECISION :: dabs0
  REAL(dp) :: temp
  REAL(dp) :: temp9
  REAL(dp) :: temp8
  REAL(dp) :: temp7
  REAL(dp) :: temp6
  REAL(dp) :: temp5
  REAL(dp) :: temp4
!!!!!!!!!!!!!!!!!!!!!!!!!!!! Wake Centers and Diameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
  spline_bound = 1.0_dp
  yaw = yaw_deg*pi/180.0_dp
! calculate y-locations of wake centers in wind ref. frame
  wakecentersyt_mat = 0.0_dp
  DO turb=1,nturbines
    CALL PUSHREAL4ARRAY(wakeangleinit, dp/4)
    wakeangleinit = 0.5_dp*SIN(yaw(turb))*ct(turb)
    IF (usewakeangle) wakeangleinit = wakeangleinit + initialwakeangle*&
&       pi/180.0_dp
    DO turbi=1,nturbines
      IF (turbinexw(turb) .LT. turbinexw(turbi)) THEN
        CALL PUSHREAL4ARRAY(deltax, dp/4)
        deltax = turbinexw(turbi) - turbinexw(turb)
        factor = 2.0_dp*kd*deltax/rotordiameter(turb) + 1.0_dp
        wakecentersyt_mat(turbi, turb) = turbineyw(turb)
        displacement = wakeangleinit*(wakeangleinit*wakeangleinit+&
&         15.0_dp*factor*factor*factor*factor)/(30.0_dp*kd/rotordiameter&
&         (turb)*(factor*factor*factor*factor*factor))
        displacement = displacement - wakeangleinit*(wakeangleinit*&
&         wakeangleinit+15.0_dp)/(30.0_dp*kd/rotordiameter(turb))
        wakecentersyt_mat(turbi, turb) = wakecentersyt_mat(turbi, turb) &
&         + initialwakedisplacement + displacement
        IF (usewakeangle .EQV. .false.) THEN
          wakecentersyt_mat(turbi, turb) = wakecentersyt_mat(turbi, turb&
&           ) + bd*deltax
          CALL PUSHCONTROL2B(2)
        ELSE
          CALL PUSHCONTROL2B(1)
        END IF
      ELSE
        CALL PUSHCONTROL2B(0)
      END IF
    END DO
  END DO
!adjust k_e to C_T, adjusted to yaw
  ke = ke_in + kecorrct*(ct-region2ct)
  wakediameterst_mat = 0.0_dp
  DO turb=1,nturbines
    IF (adjustinitialwakediamtoyaw) THEN
      CALL PUSHREAL4ARRAY(wakediameter0, dp/4)
      wakediameter0 = rotordiameter(turb)*COS(yaw(turb))
      CALL PUSHCONTROL1B(1)
    ELSE
      CALL PUSHREAL4ARRAY(wakediameter0, dp/4)
      wakediameter0 = rotordiameter(turb)
      CALL PUSHCONTROL1B(0)
    END IF
    DO turbi=1,nturbines
! turbine separation
      CALL PUSHREAL4ARRAY(deltax, dp/4)
      deltax = turbinexw(turbi) - turbinexw(turb)
! x position of interest
      x = turbinexw(turbi)
      CALL PUSHINTEGER4(zone)
      zone = 1
! define centerpoint of spline
      zeroloc = turbinexw(turb) - wakediameter0/(2.0_dp*ke(turb)*me(zone&
&       ))
      IF (zeroloc + spline_bound*rotordiameter(turb) .LT. turbinexw(&
&         turbi)) THEN
! check this
        wakediameterst_mat(turbi, turb, zone) = 0.0_dp
        CALL PUSHCONTROL2B(0)
      ELSE IF (zeroloc - spline_bound*rotordiameter(turb) .LT. turbinexw&
&         (turbi)) THEN
!check this
!!!!!!!!!!!!!!!!!!!!!! calculate spline values !!!!!!!!!!!!!!!!!!!!!!!!!!
! position of upwind point
        x1 = zeroloc - spline_bound*rotordiameter(turb)
! diameter of upwind point
        y1 = wakediameter0 + 2.0_dp*ke(turb)*me(zone)*(x1-turbinexw(turb&
&         ))
! slope at upwind point
        dy1 = 2.0_dp*ke(turb)*me(zone)
! position of downwind point
        x2 = zeroloc + spline_bound*rotordiameter(turb)
! diameter at downwind point
        y2 = 0.0_dp
! slope at downwind point
        dy2 = 0.0_dp
! solve for the wake zone diameter and its derivative w.r.t. the downwind
! location at the point of interest
        CALL HERMITE_SPLINE(x, x1, x2, y1, dy1, y2, dy2, &
&                     wakediameterst_mat(turbi, turb, zone))
        CALL PUSHCONTROL2B(1)
      ELSE IF (turbinexw(turb) .LT. turbinexw(turbi)) THEN
        wakediameterst_mat(turbi, turb, zone) = wakediameter0 + 2.0_dp*&
&         ke(turb)*me(zone)*deltax
        CALL PUSHCONTROL2B(2)
      ELSE
        CALL PUSHCONTROL2B(3)
      END IF
      IF (turbinexw(turb) .LT. turbinexw(turbi)) THEN
        CALL PUSHINTEGER4(zone)
        zone = 2
        wakediameterst_mat(turbi, turb, zone) = wakediameter0 + 2.0_dp*&
&         ke(turb)*me(zone)*deltax
        zone = 3
        wakediameterst_mat(turbi, turb, zone) = wakediameter0 + 2.0_dp*&
&         ke(turb)*me(zone)*deltax
        CALL PUSHCONTROL1B(1)
      ELSE
        CALL PUSHCONTROL1B(0)
      END IF
    END DO
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Wake Overlap !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate relative overlap
  CALL CALCOVERLAPAREAS(nturbines, turbinexw, turbineyw, rotordiameter, &
&                 wakediameterst_mat, wakecentersyt_mat, &
&                 wakeoverlaptrel_mat)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Velocity !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! convert yaw from degrees to radians
  CALL PUSHREAL4ARRAY(yaw, dp*nturbines/4)
  yaw = yaw_deg*pi/180.0_dp
  IF (axialindprovided) THEN
    a = a_in
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL CTTOAXIALIND(ct, nturbines, a)
    CALL PUSHCONTROL1B(1)
  END IF
! adjust ke to Ct as adjusted to yaw
  ke = ke_in + kecorrct*(ct-region2ct)
  DO turb=1,nturbines
    CALL PUSHREAL4ARRAY(s, dp/4)
    s = SUM(wakeoverlaptrel_mat(turb, :, 1) + wakeoverlaptrel_mat(turb, &
&     :, 2))
    kearray(turb) = ke(turb)*(1+s*kecorrarray)
  END DO
! find effective wind speeds at downstream turbines, then predict the power of 
! downstream turbine    
  velocitiesturbines = vinf
  DO turbi=1,nturbines
    CALL PUSHREAL4ARRAY(wakeeffcoeff, dp/4)
    wakeeffcoeff = 0.0_dp
! find overlap-area weighted effect of each wake zone
    DO turb=1,nturbines
      CALL PUSHREAL4ARRAY(wakeeffcoeffperzone, dp/4)
      wakeeffcoeffperzone = 0.0_dp
      CALL PUSHREAL4ARRAY(deltax, dp/4)
      deltax = turbinexw(turbi) - turbinexw(turb)
      IF (useaubu) THEN
        CALL PUSHREAL4ARRAY(mmu, dp*3/4)
        mmu = mu/COS(au*pi/180.0_dp+bu*yaw(turb))
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
      IF (deltax .GT. 0 .AND. turbi .NE. turb) THEN
        CALL PUSHINTEGER4(zone)
        DO zone=1,3
          rmax = cos_spread*0.5_dp*(wakediameterst_mat(turbi, turb, 3)+&
&           rotordiameter(turbi))
          IF (wakecentersyt_mat(turbi, turb) - turbineyw(turbi) .GE. 0.&
&         ) THEN
            CALL PUSHREAL8(dabs0)
            dabs0 = wakecentersyt_mat(turbi, turb) - turbineyw(turbi)
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHREAL8(dabs0)
            dabs0 = -(wakecentersyt_mat(turbi, turb)-turbineyw(turbi))
            CALL PUSHCONTROL1B(1)
          END IF
          CALL PUSHREAL4ARRAY(cosfac, dp/4)
          cosfac = 0.5_dp*(1.0_dp+COS(pi*dabs0/rmax))
          IF (useaubu) THEN
            wakeeffcoeffperzone = wakeeffcoeffperzone + (cosfac*&
&             rotordiameter(turb)/(rotordiameter(turb)+2.0_dp*kearray(&
&             turb)*mmu(zone)*deltax))**2*wakeoverlaptrel_mat(turbi, &
&             turb, zone)
            CALL PUSHCONTROL1B(1)
          ELSE
            wakeeffcoeffperzone = wakeeffcoeffperzone + (cosfac*&
&             rotordiameter(turb)/(rotordiameter(turb)+2.0_dp*kearray(&
&             turb)*mu(zone)*deltax))**2*wakeoverlaptrel_mat(turbi, turb&
&             , zone)
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
        wakeeffcoeff = wakeeffcoeff + (a(turb)*wakeeffcoeffperzone)**2
        CALL PUSHCONTROL1B(1)
      ELSE
        CALL PUSHCONTROL1B(0)
      END IF
    END DO
    CALL PUSHREAL4ARRAY(wakeeffcoeff, dp/4)
    wakeeffcoeff = 1.0_dp - 2.0_dp*SQRT(wakeeffcoeff)
! multiply the inflow speed with the wake coefficients to find effective wind 
! speed at turbine
  END DO
  DO nd=1,nbdirs
    rotordiameterb(nd, :) = 0.0
    turbinexwb(nd, :) = 0.0
    turbineywb(nd, :) = 0.0
    kearrayb(nd, :) = 0.0
    yawb(nd, :) = 0.0
    wakeoverlaptrel_matb(nd, :, :, :) = 0.0
    wakediameterst_matb(nd, :, :, :) = 0.0
    wakecentersyt_matb(nd, :, :) = 0.0
    mmub(nd, :) = 0.0
    ab(nd, :) = 0.0
  END DO
  DO turbi=nturbines,1,-1
    DO nd=1,nbdirs
      wakeeffcoeffb(nd) = velocitiesturbines(turbi)*velocitiesturbinesb(&
&       nd, turbi)
      velocitiesturbinesb(nd, turbi) = wakeeffcoeff*velocitiesturbinesb(&
&       nd, turbi)
    END DO
    CALL POPREAL4ARRAY(wakeeffcoeff, dp/4)
    DO nd=1,nbdirs
      IF (wakeeffcoeff .EQ. 0.0) THEN
        wakeeffcoeffb(nd) = 0.0
      ELSE
        wakeeffcoeffb(nd) = -(2.0_dp*wakeeffcoeffb(nd)/(2.0*SQRT(&
&         wakeeffcoeff)))
      END IF
    END DO
    DO turb=nturbines,1,-1
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        DO nd=1,nbdirs
          deltaxb(nd) = 0.0
        END DO
      ELSE
        DO nd=1,nbdirs
          tempb13(nd) = 2*a(turb)*wakeeffcoeffperzone*wakeeffcoeffb(nd)
          ab(nd, turb) = ab(nd, turb) + wakeeffcoeffperzone*tempb13(nd)
          wakeeffcoeffperzoneb(nd) = a(turb)*tempb13(nd)
        END DO
        deltax = turbinexw(turbi) - turbinexw(turb)
        DO nd=1,nbdirs
          deltaxb(nd) = 0.0
        END DO
        DO zone=3,1,-1
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            temp11 = 2.0_dp*mu(zone)
            temp8 = rotordiameter(turb) + temp11*kearray(turb)*deltax
            temp10 = wakeoverlaptrel_mat(turbi, turb, zone)
            temp9 = cosfac**2*rotordiameter(turb)**2
            DO nd=1,nbdirs
              tempb11(nd) = wakeeffcoeffperzoneb(nd)/temp8**2
              tempb12(nd) = -(temp9*temp10*2*tempb11(nd)/temp8)
              cosfacb(nd) = rotordiameter(turb)**2*temp10*2*cosfac*&
&               tempb11(nd)
              rotordiameterb(nd, turb) = rotordiameterb(nd, turb) + &
&               tempb12(nd) + cosfac**2*temp10*2*rotordiameter(turb)*&
&               tempb11(nd)
              wakeoverlaptrel_matb(nd, turbi, turb, zone) = &
&               wakeoverlaptrel_matb(nd, turbi, turb, zone) + temp9*&
&               tempb11(nd)
              kearrayb(nd, turb) = kearrayb(nd, turb) + temp11*deltax*&
&               tempb12(nd)
              deltaxb(nd) = deltaxb(nd) + temp11*kearray(turb)*tempb12(&
&               nd)
            END DO
          ELSE
            temp5 = rotordiameter(turb) + 2.0_dp*kearray(turb)*deltax*&
&             mmu(zone)
            temp7 = wakeoverlaptrel_mat(turbi, turb, zone)
            temp6 = cosfac**2*rotordiameter(turb)**2
            DO nd=1,nbdirs
              tempb8(nd) = wakeeffcoeffperzoneb(nd)/temp5**2
              tempb9(nd) = -(temp6*temp7*2*tempb8(nd)/temp5)
              tempb10(nd) = 2.0_dp*mmu(zone)*tempb9(nd)
              cosfacb(nd) = rotordiameter(turb)**2*temp7*2*cosfac*tempb8&
&               (nd)
              rotordiameterb(nd, turb) = rotordiameterb(nd, turb) + &
&               tempb9(nd) + cosfac**2*temp7*2*rotordiameter(turb)*&
&               tempb8(nd)
              wakeoverlaptrel_matb(nd, turbi, turb, zone) = &
&               wakeoverlaptrel_matb(nd, turbi, turb, zone) + temp6*&
&               tempb8(nd)
              kearrayb(nd, turb) = kearrayb(nd, turb) + deltax*tempb10(&
&               nd)
              deltaxb(nd) = deltaxb(nd) + kearray(turb)*tempb10(nd)
              mmub(nd, zone) = mmub(nd, zone) + 2.0_dp*kearray(turb)*&
&               deltax*tempb9(nd)
            END DO
          END IF
          rmax = cos_spread*0.5_dp*(wakediameterst_mat(turbi, turb, 3)+&
&           rotordiameter(turbi))
          CALL POPREAL4ARRAY(cosfac, dp/4)
          DO nd=1,nbdirs
            tempb7(nd) = -(pi*SIN(pi*(dabs0/rmax))*0.5_dp*cosfacb(nd)/&
&             rmax)
            dabs0b(nd) = tempb7(nd)
            rmaxb(nd) = -(dabs0*tempb7(nd)/rmax)
          END DO
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            CALL POPREAL8(dabs0)
            DO nd=1,nbdirs
              wakecentersyt_matb(nd, turbi, turb) = wakecentersyt_matb(&
&               nd, turbi, turb) + dabs0b(nd)
              turbineywb(nd, turbi) = turbineywb(nd, turbi) - dabs0b(nd)
            END DO
          ELSE
            CALL POPREAL8(dabs0)
            DO nd=1,nbdirs
              turbineywb(nd, turbi) = turbineywb(nd, turbi) + dabs0b(nd)
              wakecentersyt_matb(nd, turbi, turb) = wakecentersyt_matb(&
&               nd, turbi, turb) - dabs0b(nd)
            END DO
          END IF
          DO nd=1,nbdirs
            tempb6(nd) = cos_spread*0.5_dp*rmaxb(nd)
            wakediameterst_matb(nd, turbi, turb, 3) = &
&             wakediameterst_matb(nd, turbi, turb, 3) + tempb6(nd)
            rotordiameterb(nd, turbi) = rotordiameterb(nd, turbi) + &
&             tempb6(nd)
          END DO
        END DO
        CALL POPINTEGER4(zone)
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        CALL POPREAL4ARRAY(mmu, dp*3/4)
        temp4 = au*pi/180.0_dp + bu*yaw(turb)
        temp3 = COS(temp4)
        DO nd=1,nbdirs
          yawb(nd, turb) = yawb(nd, turb) - SIN(temp4)*bu*SUM(-(mu*mmub(&
&           nd, :)/temp3))/temp3
        END DO
        DO nd=1,nbdirs
          mmub(nd, :) = 0.0
        END DO
      END IF
      CALL POPREAL4ARRAY(deltax, dp/4)
      CALL POPREAL4ARRAY(wakeeffcoeffperzone, dp/4)
      DO nd=1,nbdirs
        turbinexwb(nd, turbi) = turbinexwb(nd, turbi) + deltaxb(nd)
        turbinexwb(nd, turb) = turbinexwb(nd, turb) - deltaxb(nd)
      END DO
    END DO
    CALL POPREAL4ARRAY(wakeeffcoeff, dp/4)
  END DO
  DO nd=1,nbdirs
    keb(nd, :) = 0.0
  END DO
  DO turb=nturbines,1,-1
    DO nd=1,nbdirs
      keb(nd, turb) = keb(nd, turb) + (kecorrarray*s+1)*kearrayb(nd, &
&       turb)
      sb(nd) = ke(turb)*kecorrarray*kearrayb(nd, turb)
      kearrayb(nd, turb) = 0.0
      wakeoverlaptrel_matb(nd, turb, :, 1) = wakeoverlaptrel_matb(nd, &
&       turb, :, 1) + sb(nd)
      wakeoverlaptrel_matb(nd, turb, :, 2) = wakeoverlaptrel_matb(nd, &
&       turb, :, 2) + sb(nd)
    END DO
    CALL POPREAL4ARRAY(s, dp/4)
  END DO
  DO nd=1,nbdirs
    ctb(nd, :) = 0.0
    ctb(nd, :) = kecorrct*keb(nd, :)
  END DO
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) THEN
    DO nd=1,nbdirs
      a_inb(nd, :) = 0.0
      a_inb(nd, :) = ab(nd, :)
    END DO
  ELSE
    CALL CTTOAXIALIND_BV(ct, ctb, nturbines, a, ab, nbdirs)
    DO nd=1,nbdirs
      a_inb(nd, :) = 0.0
    END DO
  END IF
  DO nd=1,nbdirs
    yaw_degb(nd, :) = 0.0
    yaw_degb(nd, :) = pi*yawb(nd, :)/180.0_dp
  END DO
  CALL POPREAL4ARRAY(yaw, dp*nturbines/4)
  CALL CALCOVERLAPAREAS_BV(nturbines, turbinexw, turbineyw, turbineywb, &
&                    rotordiameter, rotordiameterb, wakediameterst_mat, &
&                    wakediameterst_matb, wakecentersyt_mat, &
&                    wakecentersyt_matb, wakeoverlaptrel_mat, &
&                    wakeoverlaptrel_matb, nbdirs)
  ke = ke_in + kecorrct*(ct-region2ct)
  DO nd=1,nbdirs
    yawb(nd, :) = 0.0
    keb(nd, :) = 0.0
  END DO
  DO turb=nturbines,1,-1
    DO nd=1,nbdirs
      wakediameter0b(nd) = 0.0
    END DO
    DO turbi=nturbines,1,-1
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        DO nd=1,nbdirs
          deltaxb(nd) = 0.0
        END DO
      ELSE
        deltax = turbinexw(turbi) - turbinexw(turb)
        zone = 3
        DO nd=1,nbdirs
          tempb4(nd) = me(zone)*2.0_dp*wakediameterst_matb(nd, turbi, &
&           turb, zone)
          wakediameter0b(nd) = wakediameter0b(nd) + wakediameterst_matb(&
&           nd, turbi, turb, zone)
          wakediameterst_matb(nd, turbi, turb, zone) = 0.0
        END DO
        zone = 2
        DO nd=1,nbdirs
          tempb5(nd) = me(zone)*2.0_dp*wakediameterst_matb(nd, turbi, &
&           turb, zone)
          keb(nd, turb) = keb(nd, turb) + deltax*tempb5(nd) + deltax*&
&           tempb4(nd)
          deltaxb(nd) = ke(turb)*tempb5(nd) + ke(turb)*tempb4(nd)
          wakediameter0b(nd) = wakediameter0b(nd) + wakediameterst_matb(&
&           nd, turbi, turb, zone)
          wakediameterst_matb(nd, turbi, turb, zone) = 0.0
        END DO
        CALL POPINTEGER4(zone)
      END IF
      CALL POPCONTROL2B(branch)
      IF (branch .LT. 2) THEN
        IF (branch .EQ. 0) THEN
          DO nd=1,nbdirs
            wakediameterst_matb(nd, turbi, turb, zone) = 0.0
          END DO
          DO nd=1,nbdirs
            xb(nd) = 0.0
            zerolocb(nd) = 0.0
          END DO
        ELSE
          zone = 1
          dy1 = 2.0_dp*ke(turb)*me(zone)
          zeroloc = turbinexw(turb) - wakediameter0/(2.0_dp*ke(turb)*me(&
&           zone))
          x1 = zeroloc - spline_bound*rotordiameter(turb)
          y1 = wakediameter0 + 2.0_dp*ke(turb)*me(zone)*(x1-turbinexw(&
&           turb))
          dy2 = 0.0_dp
          y2 = 0.0_dp
          x = turbinexw(turbi)
          x2 = zeroloc + spline_bound*rotordiameter(turb)
          CALL HERMITE_SPLINE_BV(x, xb, x1, x1b, x2, x2b, y1, y1b, dy1, &
&                          dy1b, y2, dy2, wakediameterst_mat(turbi, turb&
&                          , zone), wakediameterst_matb(1, turbi, turb, &
&                          zone), nbdirs)
          DO nd=1,nbdirs
            tempb2(nd) = me(zone)*2.0_dp*y1b(nd)
            x1b(nd) = x1b(nd) + ke(turb)*tempb2(nd)
            wakediameterst_matb(nd, turbi, turb, zone) = 0.0
            zerolocb(nd) = x1b(nd) + x2b(nd)
            rotordiameterb(nd, turb) = rotordiameterb(nd, turb) + &
&             spline_bound*x2b(nd) - spline_bound*x1b(nd)
            keb(nd, turb) = keb(nd, turb) + (x1-turbinexw(turb))*tempb2(&
&             nd) + me(zone)*2.0_dp*dy1b(nd)
            wakediameter0b(nd) = wakediameter0b(nd) + y1b(nd)
            turbinexwb(nd, turb) = turbinexwb(nd, turb) - ke(turb)*&
&             tempb2(nd)
          END DO
        END IF
      ELSE
        IF (branch .EQ. 2) THEN
          DO nd=1,nbdirs
            tempb3(nd) = me(zone)*2.0_dp*wakediameterst_matb(nd, turbi, &
&             turb, zone)
            wakediameter0b(nd) = wakediameter0b(nd) + &
&             wakediameterst_matb(nd, turbi, turb, zone)
            keb(nd, turb) = keb(nd, turb) + deltax*tempb3(nd)
            deltaxb(nd) = deltaxb(nd) + ke(turb)*tempb3(nd)
            wakediameterst_matb(nd, turbi, turb, zone) = 0.0
          END DO
        END IF
        DO nd=1,nbdirs
          xb(nd) = 0.0
          zerolocb(nd) = 0.0
        END DO
      END IF
      temp2 = 2.0_dp*me(zone)*ke(turb)
      DO nd=1,nbdirs
        turbinexwb(nd, turb) = turbinexwb(nd, turb) + zerolocb(nd)
        wakediameter0b(nd) = wakediameter0b(nd) - zerolocb(nd)/temp2
        keb(nd, turb) = keb(nd, turb) + wakediameter0*2.0_dp*me(zone)*&
&         zerolocb(nd)/temp2**2
        turbinexwb(nd, turbi) = turbinexwb(nd, turbi) + deltaxb(nd) + xb&
&         (nd)
        turbinexwb(nd, turb) = turbinexwb(nd, turb) - deltaxb(nd)
      END DO
      CALL POPINTEGER4(zone)
      CALL POPREAL4ARRAY(deltax, dp/4)
    END DO
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      CALL POPREAL4ARRAY(wakediameter0, dp/4)
      yaw = yaw_deg*pi/180.0_dp
      DO nd=1,nbdirs
        rotordiameterb(nd, turb) = rotordiameterb(nd, turb) + &
&         wakediameter0b(nd)
      END DO
    ELSE
      yaw = yaw_deg*pi/180.0_dp
      CALL POPREAL4ARRAY(wakediameter0, dp/4)
      DO nd=1,nbdirs
        rotordiameterb(nd, turb) = rotordiameterb(nd, turb) + COS(yaw(&
&         turb))*wakediameter0b(nd)
        yawb(nd, turb) = yawb(nd, turb) - rotordiameter(turb)*SIN(yaw(&
&         turb))*wakediameter0b(nd)
      END DO
    END IF
  END DO
  DO nd=1,nbdirs
    ctb(nd, :) = ctb(nd, :) + kecorrct*keb(nd, :)
  END DO
  DO turb=nturbines,1,-1
    DO nd=1,nbdirs
      wakeangleinitb(nd) = 0.0
    END DO
    DO turbi=nturbines,1,-1
      CALL POPCONTROL2B(branch)
      IF (branch .NE. 0) THEN
        IF (branch .EQ. 1) THEN
          DO nd=1,nbdirs
            deltaxb(nd) = 0.0
          END DO
        ELSE
          deltax = turbinexw(turbi) - turbinexw(turb)
          DO nd=1,nbdirs
            deltaxb(nd) = bd*wakecentersyt_matb(nd, turbi, turb)
          END DO
        END IF
        factor = 2.0_dp*kd*deltax/rotordiameter(turb) + 1.0_dp
        temp1 = 30.0_dp*kd*factor**5
        temp0 = wakeangleinit*rotordiameter(turb)
        temp = wakeangleinit**2 + 15.0_dp*factor**4
        DO nd=1,nbdirs
          displacementb(nd) = wakecentersyt_matb(nd, turbi, turb)
          tempb(nd) = -((wakeangleinit**2+15.0_dp)*displacementb(nd)/(&
&           30.0_dp*kd))
          tempb0(nd) = displacementb(nd)/temp1
          wakeangleinitb(nd) = wakeangleinitb(nd) + (temp*rotordiameter(&
&           turb)+temp0*2*wakeangleinit)*tempb0(nd) + rotordiameter(turb&
&           )*tempb(nd) - wakeangleinit**2*rotordiameter(turb)*2*&
&           displacementb(nd)/(30.0_dp*kd)
          factorb(nd) = (15.0_dp*temp0*4*factor**3-30.0_dp*kd*temp*temp0&
&           *5*factor**4/temp1)*tempb0(nd)
          turbineywb(nd, turb) = turbineywb(nd, turb) + &
&           wakecentersyt_matb(nd, turbi, turb)
          wakecentersyt_matb(nd, turbi, turb) = 0.0
          tempb1(nd) = kd*2.0_dp*factorb(nd)/rotordiameter(turb)
          rotordiameterb(nd, turb) = rotordiameterb(nd, turb) + temp*&
&           wakeangleinit*tempb0(nd) - deltax*tempb1(nd)/rotordiameter(&
&           turb) + wakeangleinit*tempb(nd)
          deltaxb(nd) = deltaxb(nd) + tempb1(nd)
          turbinexwb(nd, turbi) = turbinexwb(nd, turbi) + deltaxb(nd)
          turbinexwb(nd, turb) = turbinexwb(nd, turb) - deltaxb(nd)
        END DO
        CALL POPREAL4ARRAY(deltax, dp/4)
      END IF
    END DO
    CALL POPREAL4ARRAY(wakeangleinit, dp/4)
    DO nd=1,nbdirs
      yawb(nd, turb) = yawb(nd, turb) + ct(turb)*0.5_dp*COS(yaw(turb))*&
&       wakeangleinitb(nd)
      ctb(nd, turb) = ctb(nd, turb) + 0.5_dp*SIN(yaw(turb))*&
&       wakeangleinitb(nd)
    END DO
  END DO
  DO nd=1,nbdirs
    yaw_degb(nd, :) = yaw_degb(nd, :) + pi*yawb(nd, :)/180.0_dp
  END DO
  DO nd=1,nbdirs
    velocitiesturbinesb(nd, :) = 0.0
  END DO
END SUBROUTINE FLORIS_UNIFIED_BV

!  Differentiation of hermite_spline in reverse (adjoint) mode:
!   gradient     of useful results: y
!   with respect to varying inputs: x x0 x1 dy0 y0
! Flow field calculations have been intentionally left out to save development time.
! The flow field can be calculated using the pure python version of floris 
! This implementation is fully smooth and differentiable with the exception of a 
! discontinuity at the hub of each turbine. The discontinuity only presents issues if
! turbines are place within 1E-15 * rotor diameter of one another, which is extremely 
! unlikely during optimization if the user does not explicitly place them there.
SUBROUTINE HERMITE_SPLINE_BV(x, xb, x0, x0b, x1, x1b, y0, y0b, dy0, dy0b&
& , y1, dy1, y, yb, nbdirs)
  
!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
!dy_dx = c3*3*x**2 + c2*2*x + c1
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
! in
  REAL(dp), INTENT(IN) :: x, x0, x1, y0, dy0, y1, dy1
  REAL(dp), DIMENSION(nbdirs) :: xb, x0b, x1b, y0b, dy0b
! out
!, dy_dx
  REAL(dp) :: y
  REAL(dp), DIMENSION(nbdirs) :: yb
! local
  REAL(dp) :: c3, c2, c1, c0
  REAL(dp), DIMENSION(nbdirs) :: c3b, c2b, c1b, c0b
  INTRINSIC KIND
  INTEGER :: nd
  INTEGER :: nbdirs
  REAL(dp) :: temp3
  REAL(dp) :: temp2
  REAL(dp) :: temp1
  REAL(dp) :: temp0
  REAL(dp) :: tempb9(nbdirs)
  REAL(dp) :: tempb8(nbdirs)
  REAL(dp) :: tempb7(nbdirs)
  REAL(dp) :: tempb6(nbdirs)
  REAL(dp) :: tempb5(nbdirs)
  REAL(dp) :: tempb4(nbdirs)
  REAL(dp) :: tempb19(nbdirs)
  REAL(dp) :: tempb3(nbdirs)
  REAL(dp) :: tempb18(nbdirs)
  REAL(dp) :: tempb2(nbdirs)
  REAL(dp) :: tempb17(nbdirs)
  REAL(dp) :: tempb1(nbdirs)
  REAL(dp) :: tempb16(nbdirs)
  REAL(dp) :: tempb0(nbdirs)
  REAL(dp) :: tempb15(nbdirs)
  REAL(dp) :: tempb14(nbdirs)
  REAL(dp) :: tempb13(nbdirs)
  REAL(dp) :: tempb12(nbdirs)
  REAL(dp) :: tempb11(nbdirs)
  REAL(dp) :: tempb10(nbdirs)
  REAL(dp) :: temp16
  REAL(dp) :: temp15
  REAL(dp) :: temp14
  REAL(dp) :: temp13
  REAL(dp) :: temp12
  REAL(dp) :: temp11
  REAL(dp) :: temp10
  REAL(dp) :: tempb(nbdirs)
  REAL(dp) :: tempb30(nbdirs)
  REAL(dp) :: tempb29(nbdirs)
  REAL(dp) :: tempb28(nbdirs)
  REAL(dp) :: tempb27(nbdirs)
  REAL(dp) :: tempb26(nbdirs)
  REAL(dp) :: tempb25(nbdirs)
  REAL(dp) :: temp
  REAL(dp) :: tempb24(nbdirs)
  REAL(dp) :: tempb23(nbdirs)
  REAL(dp) :: tempb22(nbdirs)
  REAL(dp) :: temp9
  REAL(dp) :: tempb21(nbdirs)
  REAL(dp) :: temp8
  REAL(dp) :: tempb20(nbdirs)
  REAL(dp) :: temp7
  REAL(dp) :: temp6
  REAL(dp) :: temp5
  REAL(dp) :: temp4
! initialize coefficients for parametric cubic spline
  c3 = 2.0_dp*y1/(x0**3-3.0_dp*x0**2*x1+3.0_dp*x0*x1**2-x1**3) - 2.0_dp*&
&   y0/(x0**3-3.0_dp*x0**2*x1+3.0_dp*x0*x1**2-x1**3) + dy0/(x0**2-2.0_dp&
&   *x0*x1+x1**2) + dy1/(x0**2-2.0_dp*x0*x1+x1**2)
  c2 = 3.0_dp*y0*(x0+x1)/(x0**3-3.0_dp*x0**2*x1+3.0_dp*x0*x1**2-x1**3) -&
&   dy1*(2.0_dp*x0+x1)/(x0**2-2.0_dp*x0*x1+x1**2) - dy0*(x0+2.0_dp*x1)/(&
&   x0**2-2.0_dp*x0*x1+x1**2) - 3.0_dp*y1*(x0+x1)/(x0**3-3.0_dp*x0**2*x1&
&   +3.0_dp*x0*x1**2-x1**3)
  c1 = dy0*(x1**2+2.0_dp*x0*x1)/(x0**2-2.0_dp*x0*x1+x1**2) + dy1*(x0**2+&
&   2.0_dp*x1*x0)/(x0**2-2.0_dp*x0*x1+x1**2) - 6.0_dp*x0*x1*y0/(x0**3-&
&   3.0_dp*x0**2*x1+3.0_dp*x0*x1**2-x1**3) + 6.0_dp*x0*x1*y1/(x0**3-&
&   3.0_dp*x0**2*x1+3.0_dp*x0*x1**2-x1**3)
!    print *, 'c3 = ', c3
!    print *, 'c2 = ', c2
!    print *, 'c1 = ', c1
!    print *, 'c0 = ', c0
! Solve for y and dy values at the given point
  temp13 = x0**3 - 3.0_dp*x0**2*x1 + 3.0_dp*x0*x1**2 - x1**3
  temp12 = 3.0_dp*x0*x1**2 - x1**3
  temp14 = x0**3 - 3.0_dp*x0**2*x1 + 3.0_dp*x0*x1**2 - x1**3
  temp15 = x0**2 - 2.0_dp*x0*x1 + x1**2
  temp16 = x0**2 - 2.0_dp*x0*x1 + x1**2
  temp8 = x0**2 - 2.0_dp*x0*x1 + x1**2
  temp7 = x1**2 + 2.0_dp*x0*x1
  temp9 = x0**2 - 2.0_dp*x0*x1 + x1**2
  temp10 = x0**3 - 3.0_dp*x0**2*x1 + 3.0_dp*x0*x1**2 - x1**3
  temp11 = x0**3 - 3.0_dp*x0**2*x1 + 3.0_dp*x0*x1**2 - x1**3
  temp3 = x0**3 - 3.0_dp*x0**2*x1 + 3.0_dp*x0*x1**2 - x1**3
  temp4 = x0**2 - 2.0_dp*x0*x1 + x1**2
  temp5 = x0**2 - 2.0_dp*x0*x1 + x1**2
  temp6 = x0**3 - 3.0_dp*x0**2*x1 + 3.0_dp*x0*x1**2 - x1**3
  temp = x0**3 - 3.0_dp*x0**2*x1 + 3.0_dp*x0*x1**2 - x1**3
  temp0 = x0**3 - 3.0_dp*x0**2*x1 + 3.0_dp*x0*x1**2 - x1**3
  temp1 = x0**2 - 2.0_dp*x0*x1 + x1**2
  temp2 = x0**2 - 2.0_dp*x0*x1 + x1**2
  DO nd=1,nbdirs
    c3b(nd) = x**3*yb(nd)
    xb(nd) = (c1+c2*2*x+c3*3*x**2)*yb(nd)
    c2b(nd) = x**2*yb(nd)
    c1b(nd) = x*yb(nd)
    c0b(nd) = yb(nd)
    tempb(nd) = c0b(nd)/temp13
    tempb0(nd) = y0*tempb(nd)
    tempb1(nd) = -(y0*temp12*tempb(nd)/temp13)
    tempb2(nd) = -(y1*c0b(nd)/temp14)
    tempb3(nd) = -((3.0_dp*(x1*x0**2)-x0**3)*tempb2(nd)/temp14)
    tempb4(nd) = -(c0b(nd)/temp15)
    tempb5(nd) = x1**2*tempb4(nd)
    tempb6(nd) = -(x1**2*x0*dy0*tempb4(nd)/temp15)
    tempb7(nd) = -(dy1*c0b(nd)/temp16)
    tempb8(nd) = -(x0**2*x1*tempb7(nd)/temp16)
    tempb30(nd) = c1b(nd)/temp8
    tempb12(nd) = dy0*tempb30(nd)
    tempb13(nd) = -(dy0*temp7*tempb30(nd)/temp8)
    tempb14(nd) = dy1*c1b(nd)/temp9
    tempb15(nd) = -((x0**2+2.0_dp*(x1*x0))*tempb14(nd)/temp9)
    tempb16(nd) = y1*6.0_dp*c1b(nd)/temp10
    tempb17(nd) = -(x0*x1*tempb16(nd)/temp10)
    tempb9(nd) = -(6.0_dp*c1b(nd)/temp11)
    tempb18(nd) = -(x0*x1*y0*tempb9(nd)/temp11)
    tempb11(nd) = 3.0_dp*c2b(nd)/temp3
    tempb29(nd) = -(y0*(x0+x1)*tempb11(nd)/temp3)
    tempb28(nd) = -(dy1*c2b(nd)/temp4)
    tempb27(nd) = -((2.0_dp*x0+x1)*tempb28(nd)/temp4)
    tempb26(nd) = -(c2b(nd)/temp5)
    dy0b(nd) = temp7*tempb30(nd) + c3b(nd)/temp1 + (x0+2.0_dp*x1)*&
&     tempb26(nd) + x0*tempb5(nd)
    tempb25(nd) = -(dy0*(x0+2.0_dp*x1)*tempb26(nd)/temp5)
    tempb24(nd) = -(y1*3.0_dp*c2b(nd)/temp6)
    tempb23(nd) = -((x0+x1)*tempb24(nd)/temp6)
    tempb19(nd) = -(y1*2.0_dp*c3b(nd)/temp**2)
    tempb10(nd) = -(2.0_dp*c3b(nd)/temp0)
    y0b(nd) = x0*x1*tempb9(nd) + tempb10(nd) + (x0+x1)*tempb11(nd) + &
&     temp12*tempb(nd)
    tempb20(nd) = -(y0*tempb10(nd)/temp0)
    tempb21(nd) = -(dy0*c3b(nd)/temp1**2)
    tempb22(nd) = -(dy1*c3b(nd)/temp2**2)
    x0b(nd) = 2.0_dp*x1*tempb12(nd) + (2*x0-2.0_dp*x1)*tempb13(nd) + (&
&     2.0_dp*x1+2*x0)*tempb14(nd) + (2*x0-2.0_dp*x1)*tempb15(nd) + x1*&
&     tempb16(nd) + (3.0_dp*x1**2-x1*3.0_dp*2*x0+3*x0**2)*tempb17(nd) + &
&     y0*x1*tempb9(nd) + (3.0_dp*x1**2-x1*3.0_dp*2*x0+3*x0**2)*tempb18(&
&     nd) + (3.0_dp*x1**2-x1*3.0_dp*2*x0+3*x0**2)*tempb19(nd) + (3.0_dp*&
&     x1**2-x1*3.0_dp*2*x0+3*x0**2)*tempb20(nd) + (2*x0-2.0_dp*x1)*&
&     tempb21(nd) + (2*x0-2.0_dp*x1)*tempb22(nd) + (3.0_dp*x1**2-x1*&
&     3.0_dp*2*x0+3*x0**2)*tempb23(nd) + tempb24(nd) + (2*x0-2.0_dp*x1)*&
&     tempb25(nd) + dy0*tempb26(nd) + (2*x0-2.0_dp*x1)*tempb27(nd) + &
&     2.0_dp*tempb28(nd) + (3.0_dp*x1**2-x1*3.0_dp*2*x0+3*x0**2)*tempb29&
&     (nd) + y0*tempb11(nd) + (2*x0-2.0_dp*x1)*tempb8(nd) + x1*2*x0*&
&     tempb7(nd) + (2*x0-2.0_dp*x1)*tempb6(nd) + dy0*tempb5(nd) + (&
&     3.0_dp*x1**2-x1*3.0_dp*2*x0+3*x0**2)*tempb3(nd) + (x1*3.0_dp*2*x0-&
&     3*x0**2)*tempb2(nd) + (3.0_dp*x1**2-x1*3.0_dp*2*x0+3*x0**2)*tempb1&
&     (nd) + 3.0_dp*x1**2*tempb0(nd)
    x1b(nd) = (2.0_dp*x0+2*x1)*tempb12(nd) + (2*x1-2.0_dp*x0)*tempb13(nd&
&     ) + 2.0_dp*x0*tempb14(nd) + (2*x1-2.0_dp*x0)*tempb15(nd) + x0*&
&     tempb16(nd) + (x0*3.0_dp*2*x1-3*x1**2-3.0_dp*x0**2)*tempb17(nd) + &
&     y0*x0*tempb9(nd) + (x0*3.0_dp*2*x1-3*x1**2-3.0_dp*x0**2)*tempb18(&
&     nd) + (x0*3.0_dp*2*x1-3*x1**2-3.0_dp*x0**2)*tempb19(nd) + (x0*&
&     3.0_dp*2*x1-3*x1**2-3.0_dp*x0**2)*tempb20(nd) + (2*x1-2.0_dp*x0)*&
&     tempb21(nd) + (2*x1-2.0_dp*x0)*tempb22(nd) + (x0*3.0_dp*2*x1-3*x1&
&     **2-3.0_dp*x0**2)*tempb23(nd) + tempb24(nd) + (2*x1-2.0_dp*x0)*&
&     tempb25(nd) + dy0*2.0_dp*tempb26(nd) + (2*x1-2.0_dp*x0)*tempb27(nd&
&     ) + tempb28(nd) + (x0*3.0_dp*2*x1-3*x1**2-3.0_dp*x0**2)*tempb29(nd&
&     ) + y0*tempb11(nd) + (2*x1-2.0_dp*x0)*tempb8(nd) + x0**2*tempb7(nd&
&     ) + (2*x1-2.0_dp*x0)*tempb6(nd) + x0*dy0*2*x1*tempb4(nd) + (x0*&
&     3.0_dp*2*x1-3*x1**2-3.0_dp*x0**2)*tempb3(nd) + 3.0_dp*x0**2*tempb2&
&     (nd) + (x0*3.0_dp*2*x1-3*x1**2-3.0_dp*x0**2)*tempb1(nd) + (x0*&
&     3.0_dp*2*x1-3*x1**2)*tempb0(nd)
  END DO
END SUBROUTINE HERMITE_SPLINE_BV

!  Differentiation of calcoverlapareas in reverse (adjoint) mode:
!   gradient     of useful results: rotordiameter turbiney wakeoverlaptrel_mat
!                wakediameters wakecenters
!   with respect to varying inputs: rotordiameter turbiney wakediameters
!                wakecenters
SUBROUTINE CALCOVERLAPAREAS_BV(nturbines, turbinex, turbiney, turbineyb&
& , rotordiameter, rotordiameterb, wakediameters, wakediametersb, &
& wakecenters, wakecentersb, wakeoverlaptrel_mat, wakeoverlaptrel_matb, &
& nbdirs)
  
!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
! do turbI = 1, nTurbines
!         do turb = 1, nTurbines
!             do zone = 1, 3
!                 print *, "wakeOverlapTRel_mat[", turbI, ", ", turb, ", ", zone, "] = ", wakeOverlapTRel_mat(turbI, turb, zone)
!             end do
!         end do
!     end do
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
! in
  INTEGER, INTENT(IN) :: nturbines
  REAL(dp), DIMENSION(nturbines), INTENT(IN) :: turbinex, turbiney, &
& rotordiameter
  REAL(dp), DIMENSION(nbdirs, nturbines) :: turbineyb, rotordiameterb
  REAL(dp), DIMENSION(nturbines, nturbines, 3), INTENT(IN) :: &
& wakediameters
  REAL(dp), DIMENSION(nbdirs, nturbines, nturbines, 3) :: &
& wakediametersb
  REAL(dp), DIMENSION(nturbines, nturbines), INTENT(IN) :: wakecenters
  REAL(dp), DIMENSION(nbdirs, nturbines, nturbines) :: wakecentersb
! out    
  REAL(dp), DIMENSION(nturbines, nturbines, 3) :: wakeoverlaptrel_mat
  REAL(dp), DIMENSION(nbdirs, nturbines, nturbines, 3) :: &
& wakeoverlaptrel_matb
! local
  INTEGER :: turb, turbi, zone
  REAL(dp), PARAMETER :: pi=3.141592653589793_dp, tol=0.000001_dp
  REAL(dp) :: ovdyd, ovr, ovrr, ovl, ovz
  REAL(dp), DIMENSION(nbdirs) :: ovdydb, ovrb, ovrrb, ovlb, ovzb
  REAL(dp), DIMENSION(nturbines, nturbines, 3) :: wakeoverlap
  REAL(dp), DIMENSION(nbdirs, nturbines, nturbines, 3) :: &
& wakeoverlapb
  INTRINSIC KIND
  INTRINSIC ABS
  INTRINSIC SQRT
  INTRINSIC DACOS
  INTEGER :: nd
  INTEGER :: branch
  INTEGER :: nbdirs
  REAL(dp) :: temp1
  REAL(dp) :: temp0
  REAL(dp) :: tempb2(nbdirs, nturbines, 3)
  REAL(dp) :: tempb1(nbdirs)
  REAL(dp) :: tempb0(nbdirs)
  REAL(dp) :: tempb(nbdirs)
  REAL(dp) :: temp
  wakeoverlap = 0.0_dp
  DO turb=1,nturbines
    DO turbi=1,nturbines
      IF (turbinex(turbi) .GT. turbinex(turb)) THEN
! distance between wake center and rotor center
        CALL PUSHREAL4ARRAY(ovdyd, dp/4)
        ovdyd = wakecenters(turbi, turb) - turbiney(turbi)
! rotor diameter
        CALL PUSHREAL4ARRAY(ovr, dp/4)
        ovr = rotordiameter(turbi)/2
        DO zone=1,3
! wake diameter
          ovrr = wakediameters(turbi, turb, zone)/2.0_dp
          IF (ovdyd .GE. 0.) THEN
            CALL PUSHREAL4ARRAY(ovdyd, dp/4)
            ovdyd = ovdyd
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHREAL4ARRAY(ovdyd, dp/4)
            ovdyd = -ovdyd
            CALL PUSHCONTROL1B(1)
          END IF
          IF (ovdyd .GE. 0.0_dp + tol) THEN
! calculate the distance from the wake center to the vertical line between
! the two circle intersection points
            CALL PUSHREAL4ARRAY(ovl, dp/4)
            ovl = (-(ovr*ovr)+ovrr*ovrr+ovdyd*ovdyd)/(2.0_dp*ovdyd)
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHREAL4ARRAY(ovl, dp/4)
            ovl = 0.0_dp
            CALL PUSHCONTROL1B(1)
          END IF
          CALL PUSHREAL4ARRAY(ovz, dp/4)
          ovz = ovrr*ovrr - ovl*ovl
! Finish calculating the distance from the intersection line to the outer edge of the wake zone
          IF (ovz .GT. 0.0_dp + tol) THEN
            CALL PUSHREAL4ARRAY(ovz, dp/4)
            ovz = SQRT(ovz)
            CALL PUSHCONTROL1B(0)
          ELSE
            ovz = 0.0_dp
            CALL PUSHCONTROL1B(1)
          END IF
          IF (ovdyd .LT. ovr + ovrr) THEN
! if the rotor overlaps the wake zone
            IF (ovl .LT. ovrr .AND. ovdyd - ovl .LT. ovr) THEN
              wakeoverlap(turbi, turb, zone) = ovrr*ovrr*DACOS(ovl/ovrr)&
&               + ovr*ovr*DACOS((ovdyd-ovl)/ovr) - ovdyd*ovz
              CALL PUSHCONTROL2B(3)
            ELSE IF (ovrr .GT. ovr) THEN
              wakeoverlap(turbi, turb, zone) = pi*ovr*ovr
              CALL PUSHCONTROL2B(2)
            ELSE
              wakeoverlap(turbi, turb, zone) = pi*ovrr*ovrr
              CALL PUSHCONTROL2B(1)
            END IF
          ELSE
            wakeoverlap(turbi, turb, zone) = 0.0_dp
            CALL PUSHCONTROL2B(0)
          END IF
        END DO
        CALL PUSHCONTROL1B(1)
      ELSE
        CALL PUSHCONTROL1B(0)
      END IF
    END DO
  END DO
  DO turb=1,nturbines
    DO turbi=1,nturbines
      wakeoverlap(turbi, turb, 3) = wakeoverlap(turbi, turb, 3) - &
&       wakeoverlap(turbi, turb, 2)
      wakeoverlap(turbi, turb, 2) = wakeoverlap(turbi, turb, 2) - &
&       wakeoverlap(turbi, turb, 1)
    END DO
  END DO
  wakeoverlaptrel_mat = wakeoverlap
  DO turbi=nturbines,1,-1
    temp1 = pi*rotordiameter(turbi)**2
    DO nd=1,nbdirs
      tempb2(nd, :, :) = 4.0_dp*wakeoverlaptrel_matb(nd, turbi, :, :)/&
&       temp1
      rotordiameterb(nd, turbi) = rotordiameterb(nd, turbi) + pi*2*&
&       rotordiameter(turbi)*SUM(-(wakeoverlaptrel_mat(turbi, :, :)*&
&       tempb2(nd, :, :)/temp1))
      wakeoverlaptrel_matb(nd, turbi, :, :) = tempb2(nd, :, :)
    END DO
  END DO
  DO nd=1,nbdirs
    wakeoverlapb(nd, :, :, :) = 0.0
    wakeoverlapb(nd, :, :, :) = wakeoverlaptrel_matb(nd, :, :, :)
  END DO
  DO turb=nturbines,1,-1
    DO turbi=nturbines,1,-1
      DO nd=1,nbdirs
        wakeoverlapb(nd, turbi, turb, 1) = wakeoverlapb(nd, turbi, turb&
&         , 1) - wakeoverlapb(nd, turbi, turb, 2)
        wakeoverlapb(nd, turbi, turb, 2) = wakeoverlapb(nd, turbi, turb&
&         , 2) - wakeoverlapb(nd, turbi, turb, 3)
      END DO
    END DO
  END DO
  DO turb=nturbines,1,-1
    DO turbi=nturbines,1,-1
      CALL POPCONTROL1B(branch)
      IF (branch .NE. 0) THEN
        DO nd=1,nbdirs
          ovdydb(nd) = 0.0
          ovrb(nd) = 0.0
        END DO
        DO zone=3,1,-1
          CALL POPCONTROL2B(branch)
          IF (branch .LT. 2) THEN
            IF (branch .EQ. 0) THEN
              DO nd=1,nbdirs
                wakeoverlapb(nd, turbi, turb, zone) = 0.0
              END DO
              ovrr = wakediameters(turbi, turb, zone)/2.0_dp
              DO nd=1,nbdirs
                ovlb(nd) = 0.0
                ovrrb(nd) = 0.0
                ovzb(nd) = 0.0
              END DO
              GOTO 100
            ELSE
              ovrr = wakediameters(turbi, turb, zone)/2.0_dp
              DO nd=1,nbdirs
                ovrrb(nd) = pi*2*ovrr*wakeoverlapb(nd, turbi, turb, zone&
&                 )
                wakeoverlapb(nd, turbi, turb, zone) = 0.0
              END DO
            END IF
          ELSE IF (branch .EQ. 2) THEN
            DO nd=1,nbdirs
              ovrb(nd) = ovrb(nd) + pi*2*ovr*wakeoverlapb(nd, turbi, &
&               turb, zone)
              wakeoverlapb(nd, turbi, turb, zone) = 0.0
            END DO
            ovrr = wakediameters(turbi, turb, zone)/2.0_dp
            DO nd=1,nbdirs
              ovrrb(nd) = 0.0
            END DO
          ELSE
            ovrr = wakediameters(turbi, turb, zone)/2.0_dp
            temp = ovl/ovrr
            temp0 = (ovdyd-ovl)/ovr
            DO nd=1,nbdirs
              IF (temp .EQ. 1.0 .OR. temp .EQ. (-1.0)) THEN
                tempb0(nd) = 0.0
              ELSE
                tempb0(nd) = -(ovrr*wakeoverlapb(nd, turbi, turb, zone)/&
&                 SQRT(1.D0-temp**2))
              END IF
              IF (temp0 .EQ. 1.0 .OR. temp0 .EQ. (-1.0)) THEN
                tempb1(nd) = 0.0
              ELSE
                tempb1(nd) = -(ovr*wakeoverlapb(nd, turbi, turb, zone)/&
&                 SQRT(1.D0-temp0**2))
              END IF
              ovrrb(nd) = DACOS(temp)*2*ovrr*wakeoverlapb(nd, turbi, &
&               turb, zone) - temp*tempb0(nd)
              ovlb(nd) = tempb0(nd) - tempb1(nd)
              ovrb(nd) = ovrb(nd) + DACOS(temp0)*2*ovr*wakeoverlapb(nd, &
&               turbi, turb, zone) - temp0*tempb1(nd)
              ovdydb(nd) = ovdydb(nd) + tempb1(nd) - ovz*wakeoverlapb(nd&
&               , turbi, turb, zone)
              ovzb(nd) = -(ovdyd*wakeoverlapb(nd, turbi, turb, zone))
              wakeoverlapb(nd, turbi, turb, zone) = 0.0
            END DO
            GOTO 100
          END IF
          DO nd=1,nbdirs
            ovlb(nd) = 0.0
            ovzb(nd) = 0.0
          END DO
 100      CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            CALL POPREAL4ARRAY(ovz, dp/4)
            DO nd=1,nbdirs
              IF (ovz .EQ. 0.0) THEN
                ovzb(nd) = 0.0
              ELSE
                ovzb(nd) = ovzb(nd)/(2.0*SQRT(ovz))
              END IF
            END DO
          ELSE
            DO nd=1,nbdirs
              ovzb(nd) = 0.0
            END DO
          END IF
          CALL POPREAL4ARRAY(ovz, dp/4)
          DO nd=1,nbdirs
            ovrrb(nd) = ovrrb(nd) + 2*ovrr*ovzb(nd)
            ovlb(nd) = ovlb(nd) - 2*ovl*ovzb(nd)
          END DO
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            CALL POPREAL4ARRAY(ovl, dp/4)
            DO nd=1,nbdirs
              tempb(nd) = ovlb(nd)/(2.0_dp*ovdyd)
              ovrrb(nd) = ovrrb(nd) + 2*ovrr*tempb(nd)
              ovrb(nd) = ovrb(nd) - 2*ovr*tempb(nd)
              ovdydb(nd) = ovdydb(nd) + (2*ovdyd-(ovrr**2-ovr**2+ovdyd**&
&               2)/ovdyd)*tempb(nd)
            END DO
          ELSE
            CALL POPREAL4ARRAY(ovl, dp/4)
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            CALL POPREAL4ARRAY(ovdyd, dp/4)
          ELSE
            CALL POPREAL4ARRAY(ovdyd, dp/4)
            DO nd=1,nbdirs
              ovdydb(nd) = -ovdydb(nd)
            END DO
          END IF
          DO nd=1,nbdirs
            wakediametersb(nd, turbi, turb, zone) = wakediametersb(nd, &
&             turbi, turb, zone) + ovrrb(nd)/2.0_dp
          END DO
        END DO
        CALL POPREAL4ARRAY(ovr, dp/4)
        CALL POPREAL4ARRAY(ovdyd, dp/4)
        DO nd=1,nbdirs
          rotordiameterb(nd, turbi) = rotordiameterb(nd, turbi) + ovrb(&
&           nd)/2
          wakecentersb(nd, turbi, turb) = wakecentersb(nd, turbi, turb) &
&           + ovdydb(nd)
          turbineyb(nd, turbi) = turbineyb(nd, turbi) - ovdydb(nd)
        END DO
      END IF
    END DO
  END DO
END SUBROUTINE CALCOVERLAPAREAS_BV

!  Differentiation of cttoaxialind in reverse (adjoint) mode:
!   gradient     of useful results: axial_induction ct
!   with respect to varying inputs: ct
SUBROUTINE CTTOAXIALIND_BV(ct, ctb, nturbines, axial_induction, &
& axial_inductionb, nbdirs)
  
!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
! in
  INTEGER, INTENT(IN) :: nturbines
  REAL(dp), DIMENSION(nturbines), INTENT(IN) :: ct
  REAL(dp), DIMENSION(nbdirs, nturbines) :: ctb
! local
  INTEGER :: i
! out
  REAL(dp), DIMENSION(nturbines) :: axial_induction
  REAL(dp), DIMENSION(nbdirs, nturbines) :: axial_inductionb
  INTRINSIC KIND
  INTRINSIC SQRT
  INTEGER :: nd
  INTEGER :: branch
  INTEGER :: nbdirs
! execute
  DO i=1,nturbines
    IF (ct(i) .GT. 0.96) THEN
      CALL PUSHCONTROL1B(1)
    ELSE
      CALL PUSHCONTROL1B(0)
    END IF
  END DO
  DO i=nturbines,1,-1
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      DO nd=1,nbdirs
        IF (.NOT.1.0_dp - ct(i) .EQ. 0.0) ctb(nd, i) = ctb(nd, i) + &
&           0.5_dp*axial_inductionb(nd, i)/(2.0*SQRT(1.0_dp-ct(i)))
        axial_inductionb(nd, i) = 0.0
      END DO
    ELSE
      DO nd=1,nbdirs
        IF (.NOT.0.0203_dp - 0.6427_dp*(0.889_dp-ct(i)) .EQ. 0.0) ctb(nd&
&         , i) = ctb(nd, i) + 0.6427_dp*axial_inductionb(nd, i)/(2.0*&
&           SQRT(0.0203_dp-0.6427_dp*(0.889_dp-ct(i))))
        axial_inductionb(nd, i) = 0.0
      END DO
    END IF
  END DO
END SUBROUTINE CTTOAXIALIND_BV




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.11 (r5902M) - 15 Dec 2015 09:00
!
!  Differentiation of floris_unified in forward (tangent) mode:
!   variations   of useful results: velocitiesturbines
!   with respect to varying inputs: rotordiameter turbinexw yaw_deg
!                turbineyw ct a_in
!   RW status of diff variables: rotordiameter:in turbinexw:in
!                yaw_deg:in velocitiesturbines:out turbineyw:in
!                ct:in a_in:in
SUBROUTINE FLORIS_UNIFIED_DV(nturbines, turbinexw, turbinexwd, turbineyw&
& , turbineywd, yaw_deg, yaw_degd, rotordiameter, rotordiameterd, vinf, &
& ct, ctd, a_in, a_ind, ke_in, kd, me, initialwakedisplacement, bd, mu, &
& au, bu, initialwakeangle, cos_spread, kecorrct, region2ct, kecorrarray&
& , usewakeangle, adjustinitialwakediamtoyaw, axialindprovided, useaubu&
& , velocitiesturbinesd, nbdirs)
  
!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
! in
  INTEGER, INTENT(IN) :: nturbines
  REAL(dp), INTENT(IN) :: kd, initialwakedisplacement, initialwakeangle&
& , ke_in
  REAL(dp), INTENT(IN) :: kecorrct, region2ct, bd, cos_spread, vinf, &
& kecorrarray
  REAL(dp), DIMENSION(nturbines), INTENT(IN) :: yaw_deg, ct, a_in, &
& turbinexw, turbineyw
  REAL(dp), DIMENSION(nbdirs, nturbines), INTENT(IN) :: yaw_degd, ctd&
& , a_ind, turbinexwd, turbineywd
  REAL(dp), DIMENSION(nturbines), INTENT(IN) :: rotordiameter
  REAL(dp), DIMENSION(nbdirs, nturbines), INTENT(IN) :: &
& rotordiameterd
  REAL(dp), DIMENSION(3), INTENT(IN) :: me, mu
  REAL(dp), INTENT(IN) :: au, bu
  LOGICAL, INTENT(IN) :: usewakeangle, adjustinitialwakediamtoyaw, &
& axialindprovided, useaubu
! local (General)
  REAL(dp), DIMENSION(nturbines) :: ke, yaw
  REAL(dp), DIMENSION(nbdirs, nturbines) :: ked, yawd
  REAL(dp) :: deltax
  REAL(dp), DIMENSION(nbdirs) :: deltaxd
  INTEGER :: turb, turbi, zone
  REAL(dp), PARAMETER :: pi=3.141592653589793_dp
! local (Wake centers and diameters)
! in rotor diameters    
  REAL(dp) :: spline_bound
  REAL(dp) :: wakeangleinit, zeroloc
  REAL(dp), DIMENSION(nbdirs) :: wakeangleinitd, zerolocd
  REAL(dp) :: factor, displacement, x, x1, x2, y1, y2, dy1, dy2
  REAL(dp), DIMENSION(nbdirs) :: factord, displacementd, xd, x1d, x2d&
& , y1d, dy1d
  REAL(dp) :: wakediameter0
  REAL(dp), DIMENSION(nbdirs) :: wakediameter0d
  REAL(dp), DIMENSION(nturbines, nturbines, 3) :: wakediameterst_mat
  REAL(dp), DIMENSION(nbdirs, nturbines, nturbines, 3) :: &
& wakediameterst_matd
  REAL(dp), DIMENSION(nturbines, nturbines) :: wakecentersyt_mat
  REAL(dp), DIMENSION(nbdirs, nturbines, nturbines) :: &
& wakecentersyt_matd
! local (Wake overlap)
  REAL(dp) :: rmax
  REAL(dp), DIMENSION(nbdirs) :: rmaxd
  REAL(dp), DIMENSION(nturbines, nturbines, 3) :: wakeoverlaptrel_mat
  REAL(dp), DIMENSION(nbdirs, nturbines, nturbines, 3) :: &
& wakeoverlaptrel_matd
  REAL(dp), DIMENSION(nturbines, nturbines, 3) :: cosfac_mat
  REAL(dp), DIMENSION(nbdirs, nturbines, nturbines, 3) :: cosfac_matd
! local (Velocity)
  REAL(dp), DIMENSION(nturbines) :: a, kearray
  REAL(dp), DIMENSION(nbdirs, nturbines) :: ad, kearrayd
  REAL(dp), DIMENSION(3) :: mmu
  REAL(dp), DIMENSION(nbdirs, 3) :: mmud
  REAL(dp) :: s, wakeeffcoeff, wakeeffcoeffperzone
  REAL(dp), DIMENSION(nbdirs) :: sd, wakeeffcoeffd, &
& wakeeffcoeffperzoned
! out
  REAL(dp), DIMENSION(nturbines) :: velocitiesturbines
  REAL(dp), DIMENSION(nbdirs, nturbines), INTENT(OUT) :: &
& velocitiesturbinesd
  REAL(dp), DIMENSION(nturbines*nturbines) :: &
& wakecentersyt_vec
  REAL(dp), DIMENSION(3*nturbines*nturbines) :: &
& wakediameterst_vec
  REAL(dp), DIMENSION(3*nturbines*nturbines) :: &
& wakeoverlaptrel_vec
  INTRINSIC COS
  INTRINSIC KIND
  INTRINSIC SIN
  INTRINSIC DABS
  INTRINSIC SUM
  INTRINSIC SQRT
  REAL(dp) :: arg1
  REAL(dp), DIMENSION(nbdirs) :: arg1d
  REAL(dp), DIMENSION(nturbines) :: arg10
  REAL(dp), DIMENSION(nbdirs, nturbines) :: arg10d
  REAL(dp) :: result1
  REAL(dp), DIMENSION(nbdirs) :: result1d
  INTEGER :: nd
  INTEGER :: nbdirs
  DOUBLE PRECISION :: dabs0d(nbdirs)
  DOUBLE PRECISION :: dabs0
!!!!!!!!!!!!!!!!!!!!!!!!!!!! Wake Centers and Diameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
  spline_bound = 1.0_dp
  DO nd=1,nbdirs
    yawd(nd, :) = pi*yaw_degd(nd, :)/180.0_dp
  END DO
  yaw = yaw_deg*pi/180.0_dp
! calculate y-locations of wake centers in wind ref. frame
  wakecentersyt_mat = 0.0_dp
  DO nd=1,nbdirs
    wakecentersyt_matd(nd, :, :) = 0.0
  END DO
  DO turb=1,nturbines
    DO nd=1,nbdirs
      wakeangleinitd(nd) = 0.5_dp*(yawd(nd, turb)*COS(yaw(turb))*ct(turb&
&       )+SIN(yaw(turb))*ctd(nd, turb))
    END DO
    wakeangleinit = 0.5_dp*SIN(yaw(turb))*ct(turb)
    IF (usewakeangle) wakeangleinit = wakeangleinit + initialwakeangle*&
&       pi/180.0_dp
    DO turbi=1,nturbines
      IF (turbinexw(turb) .LT. turbinexw(turbi)) THEN
        deltax = turbinexw(turbi) - turbinexw(turb)
        factor = 2.0_dp*kd*deltax/rotordiameter(turb) + 1.0_dp
        DO nd=1,nbdirs
          deltaxd(nd) = turbinexwd(nd, turbi) - turbinexwd(nd, turb)
          factord(nd) = (2.0_dp*kd*deltaxd(nd)*rotordiameter(turb)-&
&           2.0_dp*kd*deltax*rotordiameterd(nd, turb))/rotordiameter(&
&           turb)**2
          wakecentersyt_matd(nd, turbi, turb) = turbineywd(nd, turb)
          displacementd(nd) = ((wakeangleinitd(nd)*(wakeangleinit*&
&           wakeangleinit+15.0_dp*factor*factor*factor*factor)+&
&           wakeangleinit*(wakeangleinitd(nd)*wakeangleinit+&
&           wakeangleinit*wakeangleinitd(nd)+15.0_dp*((factord(nd)*&
&           factor+factor*factord(nd))*factor**2+factor**2*(factord(nd)*&
&           factor+factor*factord(nd)))))*30.0_dp*kd*factor**5/&
&           rotordiameter(turb)-wakeangleinit*(wakeangleinit*&
&           wakeangleinit+15.0_dp*factor*factor*factor*factor)*(30.0_dp*&
&           kd*(((factord(nd)*factor+factor*factord(nd))*factor+factor**&
&           2*factord(nd))*factor**2+factor**3*(factord(nd)*factor+&
&           factor*factord(nd)))/rotordiameter(turb)-30.0_dp*kd*&
&           rotordiameterd(nd, turb)*factor**5/rotordiameter(turb)**2))/&
&           (30.0_dp*kd/rotordiameter(turb)*(factor*factor*factor*factor&
&           *factor))**2
          displacementd(nd) = displacementd(nd) - ((wakeangleinitd(nd)*(&
&           wakeangleinit*wakeangleinit+15.0_dp)+wakeangleinit*(&
&           wakeangleinitd(nd)*wakeangleinit+wakeangleinit*&
&           wakeangleinitd(nd)))*30.0_dp*kd/rotordiameter(turb)+&
&           wakeangleinit*(wakeangleinit*wakeangleinit+15.0_dp)*30.0_dp*&
&           kd*rotordiameterd(nd, turb)/rotordiameter(turb)**2)/(30.0_dp&
&           *kd/rotordiameter(turb))**2
          wakecentersyt_matd(nd, turbi, turb) = wakecentersyt_matd(nd, &
&           turbi, turb) + displacementd(nd)
        END DO
        wakecentersyt_mat(turbi, turb) = turbineyw(turb)
        displacement = wakeangleinit*(wakeangleinit*wakeangleinit+&
&         15.0_dp*factor*factor*factor*factor)/(30.0_dp*kd/rotordiameter&
&         (turb)*(factor*factor*factor*factor*factor))
        displacement = displacement - wakeangleinit*(wakeangleinit*&
&         wakeangleinit+15.0_dp)/(30.0_dp*kd/rotordiameter(turb))
        wakecentersyt_mat(turbi, turb) = wakecentersyt_mat(turbi, turb) &
&         + initialwakedisplacement + displacement
        IF (usewakeangle .EQV. .false.) THEN
          DO nd=1,nbdirs
            wakecentersyt_matd(nd, turbi, turb) = wakecentersyt_matd(nd&
&             , turbi, turb) + bd*deltaxd(nd)
          END DO
          wakecentersyt_mat(turbi, turb) = wakecentersyt_mat(turbi, turb&
&           ) + bd*deltax
        END IF
      END IF
    END DO
  END DO
  DO nd=1,nbdirs
!adjust k_e to C_T, adjusted to yaw
    ked(nd, :) = kecorrct*ctd(nd, :)
  END DO
  ke = ke_in + kecorrct*(ct-region2ct)
  wakediameterst_mat = 0.0_dp
  DO nd=1,nbdirs
    wakediameterst_matd(nd, :, :, :) = 0.0
  END DO
  DO turb=1,nturbines
    IF (adjustinitialwakediamtoyaw) THEN
      DO nd=1,nbdirs
        wakediameter0d(nd) = rotordiameterd(nd, turb)*COS(yaw(turb)) - &
&         rotordiameter(turb)*yawd(nd, turb)*SIN(yaw(turb))
      END DO
      wakediameter0 = rotordiameter(turb)*COS(yaw(turb))
    ELSE
      DO nd=1,nbdirs
        wakediameter0d(nd) = rotordiameterd(nd, turb)
      END DO
      wakediameter0 = rotordiameter(turb)
    END IF
    DO turbi=1,nturbines
      zone = 1
      DO nd=1,nbdirs
! turbine separation
        deltaxd(nd) = turbinexwd(nd, turbi) - turbinexwd(nd, turb)
! x position of interest
        xd(nd) = turbinexwd(nd, turbi)
! define centerpoint of spline
        zerolocd(nd) = turbinexwd(nd, turb) - (wakediameter0d(nd)*2.0_dp&
&         *ke(turb)*me(zone)-wakediameter0*2.0_dp*me(zone)*ked(nd, turb)&
&         )/(2.0_dp*ke(turb)*me(zone))**2
      END DO
      deltax = turbinexw(turbi) - turbinexw(turb)
      x = turbinexw(turbi)
      zeroloc = turbinexw(turb) - wakediameter0/(2.0_dp*ke(turb)*me(zone&
&       ))
      IF (zeroloc + spline_bound*rotordiameter(turb) .LT. turbinexw(&
&         turbi)) THEN
        DO nd=1,nbdirs
! check this
          wakediameterst_matd(nd, turbi, turb, zone) = 0.0
        END DO
        wakediameterst_mat(turbi, turb, zone) = 0.0_dp
      ELSE IF (zeroloc - spline_bound*rotordiameter(turb) .LT. turbinexw&
&         (turbi)) THEN
        x1 = zeroloc - spline_bound*rotordiameter(turb)
        DO nd=1,nbdirs
!check this
!!!!!!!!!!!!!!!!!!!!!! calculate spline values !!!!!!!!!!!!!!!!!!!!!!!!!!
! position of upwind point
          x1d(nd) = zerolocd(nd) - spline_bound*rotordiameterd(nd, turb)
! diameter of upwind point
          y1d(nd) = wakediameter0d(nd) + 2.0_dp*me(zone)*(ked(nd, turb)*&
&           (x1-turbinexw(turb))+ke(turb)*(x1d(nd)-turbinexwd(nd, turb))&
&           )
! slope at upwind point
          dy1d(nd) = 2.0_dp*me(zone)*ked(nd, turb)
! position of downwind point
          x2d(nd) = zerolocd(nd) + spline_bound*rotordiameterd(nd, turb)
        END DO
        y1 = wakediameter0 + 2.0_dp*ke(turb)*me(zone)*(x1-turbinexw(turb&
&         ))
        dy1 = 2.0_dp*ke(turb)*me(zone)
        x2 = zeroloc + spline_bound*rotordiameter(turb)
! diameter at downwind point
        y2 = 0.0_dp
! slope at downwind point
        dy2 = 0.0_dp
! solve for the wake zone diameter and its derivative w.r.t. the downwind
! location at the point of interest
        CALL HERMITE_SPLINE_DV(x, xd, x1, x1d, x2, x2d, y1, y1d, dy1, &
&                        dy1d, y2, dy2, wakediameterst_mat(turbi, turb, &
&                        zone), wakediameterst_matd(:, turbi, turb, zone&
&                        ), nbdirs)
      ELSE IF (turbinexw(turb) .LT. turbinexw(turbi)) THEN
        DO nd=1,nbdirs
          wakediameterst_matd(nd, turbi, turb, zone) = wakediameter0d(nd&
&           ) + 2.0_dp*me(zone)*(ked(nd, turb)*deltax+ke(turb)*deltaxd(&
&           nd))
        END DO
        wakediameterst_mat(turbi, turb, zone) = wakediameter0 + 2.0_dp*&
&         ke(turb)*me(zone)*deltax
      END IF
      IF (turbinexw(turb) .LT. turbinexw(turbi)) THEN
        zone = 2
        DO nd=1,nbdirs
          wakediameterst_matd(nd, turbi, turb, zone) = wakediameter0d(nd&
&           ) + 2.0_dp*me(zone)*(ked(nd, turb)*deltax+ke(turb)*deltaxd(&
&           nd))
        END DO
        wakediameterst_mat(turbi, turb, zone) = wakediameter0 + 2.0_dp*&
&         ke(turb)*me(zone)*deltax
        zone = 3
        DO nd=1,nbdirs
          wakediameterst_matd(nd, turbi, turb, zone) = wakediameter0d(nd&
&           ) + 2.0_dp*me(zone)*(ked(nd, turb)*deltax+ke(turb)*deltaxd(&
&           nd))
        END DO
        wakediameterst_mat(turbi, turb, zone) = wakediameter0 + 2.0_dp*&
&         ke(turb)*me(zone)*deltax
      END IF
    END DO
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Wake Overlap !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate relative overlap
  CALL CALCOVERLAPAREAS_DV(nturbines, turbinexw, turbineyw, turbineywd, &
&                    rotordiameter, rotordiameterd, wakediameterst_mat, &
&                    wakediameterst_matd, wakecentersyt_mat, &
&                    wakecentersyt_matd, wakeoverlaptrel_mat, &
&                    wakeoverlaptrel_matd, nbdirs)
  DO nd=1,nbdirs
    cosfac_matd(nd, :, :, :) = 0.0
  END DO
! calculate cosine factor TODO: put this into the same loop as velocity
  DO turbi=1,nturbines
    DO turb=1,nturbines
      IF (turbinexw(turb) .LT. turbinexw(turbi)) THEN
        DO zone=1,3
          DO nd=1,nbdirs
            rmaxd(nd) = cos_spread*0.5_dp*(wakediameterst_matd(nd, turbi&
&             , turb, 3)+rotordiameterd(nd, turbi))
          END DO
          rmax = cos_spread*0.5_dp*(wakediameterst_mat(turbi, turb, 3)+&
&           rotordiameter(turbi))
          IF (wakecentersyt_mat(turbi, turb) - turbineyw(turbi) .GE. 0.&
&         ) THEN
            DO nd=1,nbdirs
              dabs0d(nd) = wakecentersyt_matd(nd, turbi, turb) - &
&               turbineywd(nd, turbi)
            END DO
            dabs0 = wakecentersyt_mat(turbi, turb) - turbineyw(turbi)
          ELSE
            DO nd=1,nbdirs
              dabs0d(nd) = -(wakecentersyt_matd(nd, turbi, turb)-&
&               turbineywd(nd, turbi))
            END DO
            dabs0 = -(wakecentersyt_mat(turbi, turb)-turbineyw(turbi))
          END IF
          arg1 = pi*dabs0/rmax
          DO nd=1,nbdirs
            arg1d(nd) = (pi*dabs0d(nd)*rmax-pi*dabs0*rmaxd(nd))/rmax**2
            cosfac_matd(nd, turbi, turb, zone) = -(0.5_dp*arg1d(nd)*SIN(&
&             arg1))
          END DO
          cosfac_mat(turbi, turb, zone) = 0.5_dp*(1.0_dp+COS(arg1))
        END DO
      ELSE
        DO nd=1,nbdirs
!cosFac_mat(turbI, turb, zone) = 1.0_dp
          cosfac_matd(nd, turbi, turb, :) = 0.0
        END DO
        cosfac_mat(turbi, turb, :) = 1.0_dp
      END IF
    END DO
  END DO
  DO nd=1,nbdirs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Velocity !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! convert yaw from degrees to radians
    yawd(nd, :) = pi*yaw_degd(nd, :)/180.0_dp
  END DO
  yaw = yaw_deg*pi/180.0_dp
  IF (axialindprovided) THEN
    DO nd=1,nbdirs
      ad(nd, :) = a_ind(nd, :)
    END DO
    a = a_in
  ELSE
    CALL CTTOAXIALIND_DV(ct, ctd, nturbines, a, ad, nbdirs)
  END IF
  DO nd=1,nbdirs
! adjust ke to Ct as adjusted to yaw
    ked(nd, :) = kecorrct*ctd(nd, :)
  END DO
  ke = ke_in + kecorrct*(ct-region2ct)
  DO nd=1,nbdirs
    kearrayd(nd, :) = 0.0
  END DO
  DO turb=1,nturbines
    arg10(:) = wakeoverlaptrel_mat(turb, :, 1) + wakeoverlaptrel_mat(&
&     turb, :, 2)
    s = SUM(arg10(:))
    DO nd=1,nbdirs
      arg10d(nd, :) = wakeoverlaptrel_matd(nd, turb, :, 1) + &
&       wakeoverlaptrel_matd(nd, turb, :, 2)
      sd(nd) = SUM(arg10d(nd, :))
      kearrayd(nd, turb) = ked(nd, turb)*(1+s*kecorrarray) + ke(turb)*&
&       kecorrarray*sd(nd)
    END DO
    kearray(turb) = ke(turb)*(1+s*kecorrarray)
  END DO
! find effective wind speeds at downstream turbines, then predict the power of 
! downstream turbine    
  velocitiesturbines = vinf
  DO nd=1,nbdirs
    velocitiesturbinesd(nd, :) = 0.0
    mmud(nd, :) = 0.0
  END DO
  DO turbi=1,nturbines
    wakeeffcoeff = 0.0_dp
    DO nd=1,nbdirs
      wakeeffcoeffd(nd) = 0.0
    END DO
! find overlap-area weighted effect of each wake zone
    DO turb=1,nturbines
      wakeeffcoeffperzone = 0.0_dp
      DO nd=1,nbdirs
        deltaxd(nd) = turbinexwd(nd, turbi) - turbinexwd(nd, turb)
      END DO
      deltax = turbinexw(turbi) - turbinexw(turb)
      IF (useaubu) THEN
        arg1 = au*pi/180.0_dp + bu*yaw(turb)
        DO nd=1,nbdirs
          arg1d(nd) = bu*yawd(nd, turb)
          mmud(nd, :) = -((-(mu*arg1d(nd)*SIN(arg1)))/COS(arg1)**2)
        END DO
        mmu = mu/COS(arg1)
      END IF
      IF (deltax .GT. 0 .AND. turbi .NE. turb) THEN
        DO nd=1,nbdirs
          wakeeffcoeffperzoned(nd) = 0.0
        END DO
        DO zone=1,3
          IF (useaubu) THEN
            DO nd=1,nbdirs
              wakeeffcoeffperzoned(nd) = wakeeffcoeffperzoned(nd) + 2*&
&               cosfac_mat(turbi, turb, zone)*rotordiameter(turb)*((&
&               cosfac_matd(nd, turbi, turb, zone)*rotordiameter(turb)+&
&               cosfac_mat(turbi, turb, zone)*rotordiameterd(nd, turb))*&
&               (rotordiameter(turb)+2.0_dp*kearray(turb)*mmu(zone)*&
&               deltax)-cosfac_mat(turbi, turb, zone)*rotordiameter(turb&
&               )*(rotordiameterd(nd, turb)+2.0_dp*((kearrayd(nd, turb)*&
&               deltax+kearray(turb)*deltaxd(nd))*mmu(zone)+kearray(turb&
&               )*deltax*mmud(nd, zone))))*wakeoverlaptrel_mat(turbi, &
&               turb, zone)/(rotordiameter(turb)+2.0_dp*kearray(turb)*&
&               mmu(zone)*deltax)**3 + cosfac_mat(turbi, turb, zone)**2*&
&               rotordiameter(turb)**2*wakeoverlaptrel_matd(nd, turbi, &
&               turb, zone)/(rotordiameter(turb)+2.0_dp*kearray(turb)*&
&               mmu(zone)*deltax)**2
            END DO
            wakeeffcoeffperzone = wakeeffcoeffperzone + (cosfac_mat(&
&             turbi, turb, zone)*rotordiameter(turb)/(rotordiameter(turb&
&             )+2.0_dp*kearray(turb)*mmu(zone)*deltax))**2*&
&             wakeoverlaptrel_mat(turbi, turb, zone)
          ELSE
            DO nd=1,nbdirs
              wakeeffcoeffperzoned(nd) = wakeeffcoeffperzoned(nd) + 2*&
&               cosfac_mat(turbi, turb, zone)*rotordiameter(turb)*((&
&               cosfac_matd(nd, turbi, turb, zone)*rotordiameter(turb)+&
&               cosfac_mat(turbi, turb, zone)*rotordiameterd(nd, turb))*&
&               (rotordiameter(turb)+2.0_dp*kearray(turb)*mu(zone)*&
&               deltax)-cosfac_mat(turbi, turb, zone)*rotordiameter(turb&
&               )*(rotordiameterd(nd, turb)+2.0_dp*mu(zone)*(kearrayd(nd&
&               , turb)*deltax+kearray(turb)*deltaxd(nd))))*&
&               wakeoverlaptrel_mat(turbi, turb, zone)/(rotordiameter(&
&               turb)+2.0_dp*kearray(turb)*mu(zone)*deltax)**3 + &
&               cosfac_mat(turbi, turb, zone)**2*rotordiameter(turb)**2*&
&               wakeoverlaptrel_matd(nd, turbi, turb, zone)/(&
&               rotordiameter(turb)+2.0_dp*kearray(turb)*mu(zone)*deltax&
&               )**2
            END DO
            wakeeffcoeffperzone = wakeeffcoeffperzone + (cosfac_mat(&
&             turbi, turb, zone)*rotordiameter(turb)/(rotordiameter(turb&
&             )+2.0_dp*kearray(turb)*mu(zone)*deltax))**2*&
&             wakeoverlaptrel_mat(turbi, turb, zone)
          END IF
        END DO
        DO nd=1,nbdirs
          wakeeffcoeffd(nd) = wakeeffcoeffd(nd) + 2*a(turb)*&
&           wakeeffcoeffperzone*(ad(nd, turb)*wakeeffcoeffperzone+a(turb&
&           )*wakeeffcoeffperzoned(nd))
        END DO
        wakeeffcoeff = wakeeffcoeff + (a(turb)*wakeeffcoeffperzone)**2
      END IF
    END DO
    DO nd=1,nbdirs
      IF (wakeeffcoeff .EQ. 0.0) THEN
        result1d(nd) = 0.0
      ELSE
        result1d(nd) = wakeeffcoeffd(nd)/(2.0*SQRT(wakeeffcoeff))
      END IF
      wakeeffcoeffd(nd) = -(2.0_dp*result1d(nd))
    END DO
    result1 = SQRT(wakeeffcoeff)
    wakeeffcoeff = 1.0_dp - 2.0_dp*result1
    DO nd=1,nbdirs
! multiply the inflow speed with the wake coefficients to find effective wind 
! speed at turbine
      velocitiesturbinesd(nd, turbi) = velocitiesturbinesd(nd, turbi)*&
&       wakeeffcoeff + velocitiesturbines(turbi)*wakeeffcoeffd(nd)
    END DO
    velocitiesturbines(turbi) = velocitiesturbines(turbi)*wakeeffcoeff
  END DO
! pack desired matrices into vectors for output
  DO turbi=1,nturbines
! wake centers
    wakecentersyt_vec(nturbines*(turbi-1)+1:nturbines*(turbi-1)+&
&   nturbines) = wakecentersyt_mat(turbi, :)
! wake diameters
    wakediameterst_vec(3*nturbines*(turbi-1)+1:3*nturbines*(turbi-1)+&
&   nturbines) = wakediameterst_mat(turbi, :, 1)
    wakediameterst_vec(3*nturbines*(turbi-1)+nturbines+1:3*nturbines*(&
&   turbi-1)+2*nturbines) = wakediameterst_mat(turbi, :, 2)
    wakediameterst_vec(3*nturbines*(turbi-1)+2*nturbines+1:nturbines*(&
&   turbi-1)+3*nturbines) = wakediameterst_mat(turbi, :, 3)
! relative wake overlap
    wakeoverlaptrel_vec(3*nturbines*(turbi-1)+1:3*nturbines*(turbi-1)+&
&   nturbines) = wakeoverlaptrel_mat(turbi, :, 1)
    wakeoverlaptrel_vec(3*nturbines*(turbi-1)+nturbines+1:3*nturbines*(&
&   turbi-1)+2*nturbines) = wakeoverlaptrel_mat(turbi, :, 2)
    wakeoverlaptrel_vec(3*nturbines*(turbi-1)+2*nturbines+1:3*nturbines*&
&   (turbi-1)+3*nturbines) = wakeoverlaptrel_mat(turbi, :, 3)
  END DO
END SUBROUTINE FLORIS_UNIFIED_DV

!  Differentiation of hermite_spline in forward (tangent) mode:
!   variations   of useful results: y
!   with respect to varying inputs: x x0 x1 dy0 y0
! Flow field calculations have been intentionally left out to save development time.
! The flow field can be calculated using the pure python version of floris 
! This implementation is fully smooth and differentiable with the exception of a 
! discontinuity at the hub of each turbine. The discontinuity only presents issues if
! turbines are place within 1E-15 * rotor diameter of one another, which is extremely 
! unlikely during optimization if the user does not explicitly place them there.
SUBROUTINE HERMITE_SPLINE_DV(x, xd, x0, x0d, x1, x1d, y0, y0d, dy0, dy0d&
& , y1, dy1, y, yd, nbdirs)
  
!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
!dy_dx = c3*3*x**2 + c2*2*x + c1
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
! in
  REAL(dp), INTENT(IN) :: x, x0, x1, y0, dy0, y1, dy1
  REAL(dp), DIMENSION(nbdirs), INTENT(IN) :: xd, x0d, x1d, y0d, dy0d
! out
!, dy_dx
  REAL(dp), INTENT(OUT) :: y
  REAL(dp), DIMENSION(nbdirs), INTENT(OUT) :: yd
! local
  REAL(dp) :: c3, c2, c1, c0
  REAL(dp), DIMENSION(nbdirs) :: c3d, c2d, c1d, c0d
  INTRINSIC KIND
  INTEGER :: nd
  INTEGER :: nbdirs
  c3 = 2.0_dp*y1/(x0**3-3.0_dp*x0**2*x1+3.0_dp*x0*x1**2-x1**3) - 2.0_dp*&
&   y0/(x0**3-3.0_dp*x0**2*x1+3.0_dp*x0*x1**2-x1**3) + dy0/(x0**2-2.0_dp&
&   *x0*x1+x1**2) + dy1/(x0**2-2.0_dp*x0*x1+x1**2)
  c2 = 3.0_dp*y0*(x0+x1)/(x0**3-3.0_dp*x0**2*x1+3.0_dp*x0*x1**2-x1**3) -&
&   dy1*(2.0_dp*x0+x1)/(x0**2-2.0_dp*x0*x1+x1**2) - dy0*(x0+2.0_dp*x1)/(&
&   x0**2-2.0_dp*x0*x1+x1**2) - 3.0_dp*y1*(x0+x1)/(x0**3-3.0_dp*x0**2*x1&
&   +3.0_dp*x0*x1**2-x1**3)
  c1 = dy0*(x1**2+2.0_dp*x0*x1)/(x0**2-2.0_dp*x0*x1+x1**2) + dy1*(x0**2+&
&   2.0_dp*x1*x0)/(x0**2-2.0_dp*x0*x1+x1**2) - 6.0_dp*x0*x1*y0/(x0**3-&
&   3.0_dp*x0**2*x1+3.0_dp*x0*x1**2-x1**3) + 6.0_dp*x0*x1*y1/(x0**3-&
&   3.0_dp*x0**2*x1+3.0_dp*x0*x1**2-x1**3)
  DO nd=1,nbdirs
! initialize coefficients for parametric cubic spline
    c3d(nd) = (dy0d(nd)*(x0**2-2.0_dp*x0*x1+x1**2)-dy0*(2*x0*x0d(nd)-&
&     2.0_dp*(x0d(nd)*x1+x0*x1d(nd))+2*x1*x1d(nd)))/(x0**2-2.0_dp*x0*x1+&
&     x1**2)**2 - (2.0_dp*y0d(nd)*(x0**3-3.0_dp*x0**2*x1+3.0_dp*x0*x1**2&
&     -x1**3)-2.0_dp*y0*(3*x0**2*x0d(nd)-3.0_dp*(2*x0*x0d(nd)*x1+x0**2*&
&     x1d(nd))+3.0_dp*(x0d(nd)*x1**2+x0*2*x1*x1d(nd))-3*x1**2*x1d(nd)))/&
&     (x0**3-3.0_dp*x0**2*x1+3.0_dp*x0*x1**2-x1**3)**2 - 2.0_dp*y1*(3*x0&
&     **2*x0d(nd)-3.0_dp*(2*x0*x0d(nd)*x1+x0**2*x1d(nd))+3.0_dp*(x0d(nd)&
&     *x1**2+x0*2*x1*x1d(nd))-3*x1**2*x1d(nd))/(x0**3-3.0_dp*x0**2*x1+&
&     3.0_dp*x0*x1**2-x1**3)**2 - dy1*(2*x0*x0d(nd)-2.0_dp*(x0d(nd)*x1+&
&     x0*x1d(nd))+2*x1*x1d(nd))/(x0**2-2.0_dp*x0*x1+x1**2)**2
    c2d(nd) = (3.0_dp*(y0d(nd)*(x0+x1)+y0*(x0d(nd)+x1d(nd)))*(x0**3-&
&     3.0_dp*x0**2*x1+3.0_dp*x0*x1**2-x1**3)-3.0_dp*y0*(x0+x1)*(3*x0**2*&
&     x0d(nd)-3.0_dp*(2*x0*x0d(nd)*x1+x0**2*x1d(nd))+3.0_dp*(x0d(nd)*x1&
&     **2+x0*2*x1*x1d(nd))-3*x1**2*x1d(nd)))/(x0**3-3.0_dp*x0**2*x1+&
&     3.0_dp*x0*x1**2-x1**3)**2 - (dy1*(2.0_dp*x0d(nd)+x1d(nd))*(x0**2-&
&     2.0_dp*x0*x1+x1**2)-dy1*(2.0_dp*x0+x1)*(2*x0*x0d(nd)-2.0_dp*(x0d(&
&     nd)*x1+x0*x1d(nd))+2*x1*x1d(nd)))/(x0**2-2.0_dp*x0*x1+x1**2)**2 - &
&     ((dy0d(nd)*(x0+2.0_dp*x1)+dy0*(x0d(nd)+2.0_dp*x1d(nd)))*(x0**2-&
&     2.0_dp*x0*x1+x1**2)-dy0*(x0+2.0_dp*x1)*(2*x0*x0d(nd)-2.0_dp*(x0d(&
&     nd)*x1+x0*x1d(nd))+2*x1*x1d(nd)))/(x0**2-2.0_dp*x0*x1+x1**2)**2 - &
&     (3.0_dp*y1*(x0d(nd)+x1d(nd))*(x0**3-3.0_dp*x0**2*x1+3.0_dp*x0*x1**&
&     2-x1**3)-3.0_dp*y1*(x0+x1)*(3*x0**2*x0d(nd)-3.0_dp*(2*x0*x0d(nd)*&
&     x1+x0**2*x1d(nd))+3.0_dp*(x0d(nd)*x1**2+x0*2*x1*x1d(nd))-3*x1**2*&
&     x1d(nd)))/(x0**3-3.0_dp*x0**2*x1+3.0_dp*x0*x1**2-x1**3)**2
    c1d(nd) = ((dy0d(nd)*(x1**2+2.0_dp*x0*x1)+dy0*(2*x1*x1d(nd)+2.0_dp*(&
&     x0d(nd)*x1+x0*x1d(nd))))*(x0**2-2.0_dp*x0*x1+x1**2)-dy0*(x1**2+&
&     2.0_dp*x0*x1)*(2*x0*x0d(nd)-2.0_dp*(x0d(nd)*x1+x0*x1d(nd))+2*x1*&
&     x1d(nd)))/(x0**2-2.0_dp*x0*x1+x1**2)**2 + (dy1*(2*x0*x0d(nd)+&
&     2.0_dp*(x1d(nd)*x0+x1*x0d(nd)))*(x0**2-2.0_dp*x0*x1+x1**2)-dy1*(x0&
&     **2+2.0_dp*x1*x0)*(2*x0*x0d(nd)-2.0_dp*(x0d(nd)*x1+x0*x1d(nd))+2*&
&     x1*x1d(nd)))/(x0**2-2.0_dp*x0*x1+x1**2)**2 - (6.0_dp*((x0d(nd)*x1+&
&     x0*x1d(nd))*y0+x0*x1*y0d(nd))*(x0**3-3.0_dp*x0**2*x1+3.0_dp*x0*x1&
&     **2-x1**3)-6.0_dp*x0*x1*y0*(3*x0**2*x0d(nd)-3.0_dp*(2*x0*x0d(nd)*&
&     x1+x0**2*x1d(nd))+3.0_dp*(x0d(nd)*x1**2+x0*2*x1*x1d(nd))-3*x1**2*&
&     x1d(nd)))/(x0**3-3.0_dp*x0**2*x1+3.0_dp*x0*x1**2-x1**3)**2 + (&
&     6.0_dp*y1*(x0d(nd)*x1+x0*x1d(nd))*(x0**3-3.0_dp*x0**2*x1+3.0_dp*x0&
&     *x1**2-x1**3)-6.0_dp*x0*x1*y1*(3*x0**2*x0d(nd)-3.0_dp*(2*x0*x0d(nd&
&     )*x1+x0**2*x1d(nd))+3.0_dp*(x0d(nd)*x1**2+x0*2*x1*x1d(nd))-3*x1**2&
&     *x1d(nd)))/(x0**3-3.0_dp*x0**2*x1+3.0_dp*x0*x1**2-x1**3)**2
    c0d(nd) = ((y0d(nd)*(-(x1**3)+3.0_dp*x0*x1**2)+y0*(3.0_dp*(x0d(nd)*&
&     x1**2+x0*2*x1*x1d(nd))-3*x1**2*x1d(nd)))*(x0**3-3.0_dp*x0**2*x1+&
&     3.0_dp*x0*x1**2-x1**3)-y0*(-(x1**3)+3.0_dp*x0*x1**2)*(3*x0**2*x0d(&
&     nd)-3.0_dp*(2*x0*x0d(nd)*x1+x0**2*x1d(nd))+3.0_dp*(x0d(nd)*x1**2+&
&     x0*2*x1*x1d(nd))-3*x1**2*x1d(nd)))/(x0**3-3.0_dp*x0**2*x1+3.0_dp*&
&     x0*x1**2-x1**3)**2 - (y1*(3.0_dp*(x1d(nd)*x0**2+x1*2*x0*x0d(nd))-3&
&     *x0**2*x0d(nd))*(x0**3-3.0_dp*x0**2*x1+3.0_dp*x0*x1**2-x1**3)-y1*(&
&     -(x0**3)+3.0_dp*x1*x0**2)*(3*x0**2*x0d(nd)-3.0_dp*(2*x0*x0d(nd)*x1&
&     +x0**2*x1d(nd))+3.0_dp*(x0d(nd)*x1**2+x0*2*x1*x1d(nd))-3*x1**2*x1d&
&     (nd)))/(x0**3-3.0_dp*x0**2*x1+3.0_dp*x0*x1**2-x1**3)**2 - (((x0d(&
&     nd)*dy0+x0*dy0d(nd))*x1**2+x0*dy0*2*x1*x1d(nd))*(x0**2-2.0_dp*x0*&
&     x1+x1**2)-x0*x1**2*dy0*(2*x0*x0d(nd)-2.0_dp*(x0d(nd)*x1+x0*x1d(nd)&
&     )+2*x1*x1d(nd)))/(x0**2-2.0_dp*x0*x1+x1**2)**2 - (dy1*(2*x0*x0d(nd&
&     )*x1+x0**2*x1d(nd))*(x0**2-2.0_dp*x0*x1+x1**2)-x0**2*x1*dy1*(2*x0*&
&     x0d(nd)-2.0_dp*(x0d(nd)*x1+x0*x1d(nd))+2*x1*x1d(nd)))/(x0**2-&
&     2.0_dp*x0*x1+x1**2)**2
!    print *, 'c3 = ', c3
!    print *, 'c2 = ', c2
!    print *, 'c1 = ', c1
!    print *, 'c0 = ', c0
! Solve for y and dy values at the given point
    yd(nd) = c3d(nd)*x**3 + c3*3*x**2*xd(nd) + c2d(nd)*x**2 + c2*2*x*xd(&
&     nd) + c1d(nd)*x + c1*xd(nd) + c0d(nd)
  END DO
  c0 = y0*(-(x1**3)+3.0_dp*x0*x1**2)/(x0**3-3.0_dp*x0**2*x1+3.0_dp*x0*x1&
&   **2-x1**3) - y1*(-(x0**3)+3.0_dp*x1*x0**2)/(x0**3-3.0_dp*x0**2*x1+&
&   3.0_dp*x0*x1**2-x1**3) - x0*x1**2*dy0/(x0**2-2.0_dp*x0*x1+x1**2) - &
&   x0**2*x1*dy1/(x0**2-2.0_dp*x0*x1+x1**2)
  y = c3*x**3 + c2*x**2 + c1*x + c0
END SUBROUTINE HERMITE_SPLINE_DV

!  Differentiation of calcoverlapareas in forward (tangent) mode:
!   variations   of useful results: wakeoverlaptrel_mat
!   with respect to varying inputs: rotordiameter turbiney wakediameters
!                wakecenters
SUBROUTINE CALCOVERLAPAREAS_DV(nturbines, turbinex, turbiney, turbineyd&
& , rotordiameter, rotordiameterd, wakediameters, wakediametersd, &
& wakecenters, wakecentersd, wakeoverlaptrel_mat, wakeoverlaptrel_matd, &
& nbdirs)
  
!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
! do turbI = 1, nTurbines
!         do turb = 1, nTurbines
!             do zone = 1, 3
!                 print *, "wakeOverlapTRel_mat[", turbI, ", ", turb, ", ", zone, "] = ", wakeOverlapTRel_mat(turbI, turb, zone)
!             end do
!         end do
!     end do
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
! in
  INTEGER, INTENT(IN) :: nturbines
  REAL(dp), DIMENSION(nturbines), INTENT(IN) :: turbinex, turbiney, &
& rotordiameter
  REAL(dp), DIMENSION(nbdirs, nturbines), INTENT(IN) :: turbineyd, &
& rotordiameterd
  REAL(dp), DIMENSION(nturbines, nturbines, 3), INTENT(IN) :: &
& wakediameters
  REAL(dp), DIMENSION(nbdirs, nturbines, nturbines, 3), INTENT(IN) ::&
& wakediametersd
  REAL(dp), DIMENSION(nturbines, nturbines), INTENT(IN) :: wakecenters
  REAL(dp), DIMENSION(nbdirs, nturbines, nturbines), INTENT(IN) :: &
& wakecentersd
! out    
  REAL(dp), DIMENSION(nturbines, nturbines, 3), INTENT(OUT) :: &
& wakeoverlaptrel_mat
  REAL(dp), DIMENSION(nbdirs, nturbines, nturbines, 3), INTENT(OUT) &
& :: wakeoverlaptrel_matd
! local
  INTEGER :: turb, turbi, zone
  REAL(dp), PARAMETER :: pi=3.141592653589793_dp, tol=0.000001_dp
  REAL(dp) :: ovdyd, ovr, ovrr, ovl, ovz
  REAL(dp), DIMENSION(nbdirs) :: ovdydd, ovrd, ovrrd, ovld, ovzd
  REAL(dp), DIMENSION(nturbines, nturbines, 3) :: wakeoverlap
  REAL(dp), DIMENSION(nbdirs, nturbines, nturbines, 3) :: &
& wakeoverlapd
  INTRINSIC KIND
  INTRINSIC ABS
  INTRINSIC SQRT
  INTRINSIC DACOS
  REAL(dp) :: arg1
  REAL(dp), DIMENSION(nbdirs) :: arg1d
  DOUBLE PRECISION :: result1
  DOUBLE PRECISION, DIMENSION(nbdirs) :: result1d
  REAL(dp) :: arg2
  REAL(dp), DIMENSION(nbdirs) :: arg2d
  DOUBLE PRECISION :: result2
  DOUBLE PRECISION, DIMENSION(nbdirs) :: result2d
  INTEGER :: nd
  INTEGER :: nbdirs
  wakeoverlaptrel_mat = 0.0_dp
  wakeoverlap = 0.0_dp
  DO nd=1,nbdirs
    wakeoverlapd(nd, :, :, :) = 0.0
  END DO
  DO turb=1,nturbines
    DO turbi=1,nturbines
      IF (turbinex(turbi) .GT. turbinex(turb)) THEN
        DO nd=1,nbdirs
! distance between wake center and rotor center
          ovdydd(nd) = wakecentersd(nd, turbi, turb) - turbineyd(nd, &
&           turbi)
! rotor diameter
          ovrd(nd) = rotordiameterd(nd, turbi)/2
        END DO
        ovdyd = wakecenters(turbi, turb) - turbiney(turbi)
        ovr = rotordiameter(turbi)/2
        DO zone=1,3
          DO nd=1,nbdirs
! wake diameter
            ovrrd(nd) = wakediametersd(nd, turbi, turb, zone)/2.0_dp
          END DO
          ovrr = wakediameters(turbi, turb, zone)/2.0_dp
          IF (ovdyd .GE. 0.) THEN
            ovdyd = ovdyd
          ELSE
            DO nd=1,nbdirs
              ovdydd(nd) = -ovdydd(nd)
            END DO
            ovdyd = -ovdyd
          END IF
          IF (ovdyd .GE. 0.0_dp + tol) THEN
            DO nd=1,nbdirs
! calculate the distance from the wake center to the vertical line between
! the two circle intersection points
              ovld(nd) = ((ovrrd(nd)*ovrr-ovr*ovrd(nd)-ovrd(nd)*ovr+ovrr&
&               *ovrrd(nd)+ovdydd(nd)*ovdyd+ovdyd*ovdydd(nd))*2.0_dp*&
&               ovdyd-(-(ovr*ovr)+ovrr*ovrr+ovdyd*ovdyd)*2.0_dp*ovdydd(&
&               nd))/(2.0_dp*ovdyd)**2
            END DO
            ovl = (-(ovr*ovr)+ovrr*ovrr+ovdyd*ovdyd)/(2.0_dp*ovdyd)
          ELSE
            ovl = 0.0_dp
            DO nd=1,nbdirs
              ovld(nd) = 0.0
            END DO
          END IF
          DO nd=1,nbdirs
            ovzd(nd) = ovrrd(nd)*ovrr + ovrr*ovrrd(nd) - ovld(nd)*ovl - &
&             ovl*ovld(nd)
          END DO
          ovz = ovrr*ovrr - ovl*ovl
! Finish calculating the distance from the intersection line to the outer edge of the wake zone
          IF (ovz .GT. 0.0_dp + tol) THEN
            DO nd=1,nbdirs
              IF (ovz .EQ. 0.0) THEN
                ovzd(nd) = 0.0
              ELSE
                ovzd(nd) = ovzd(nd)/(2.0*SQRT(ovz))
              END IF
            END DO
            ovz = SQRT(ovz)
          ELSE
            ovz = 0.0_dp
            DO nd=1,nbdirs
              ovzd(nd) = 0.0
            END DO
          END IF
          IF (ovdyd .LT. ovr + ovrr) THEN
! if the rotor overlaps the wake zone
            IF (ovl .LT. ovrr .AND. ovdyd - ovl .LT. ovr) THEN
              arg1 = ovl/ovrr
              result1 = DACOS(arg1)
              arg2 = (ovdyd-ovl)/ovr
              result2 = DACOS(arg2)
              DO nd=1,nbdirs
                arg1d(nd) = (ovld(nd)*ovrr-ovl*ovrrd(nd))/ovrr**2
                IF (arg1 .EQ. 1.0 .OR. arg1 .EQ. (-1.0)) THEN
                  result1d(nd) = 0.D0
                ELSE
                  result1d(nd) = -(arg1d(nd)/SQRT(1.D0-arg1**2))
                END IF
                arg2d(nd) = ((ovdydd(nd)-ovld(nd))*ovr-(ovdyd-ovl)*ovrd(&
&                 nd))/ovr**2
                IF (arg2 .EQ. 1.0 .OR. arg2 .EQ. (-1.0)) THEN
                  result2d(nd) = 0.D0
                ELSE
                  result2d(nd) = -(arg2d(nd)/SQRT(1.D0-arg2**2))
                END IF
                wakeoverlapd(nd, turbi, turb, zone) = (ovrrd(nd)*ovrr+&
&                 ovrr*ovrrd(nd))*result1 + ovrr**2*result1d(nd) + (ovrd&
&                 (nd)*ovr+ovr*ovrd(nd))*result2 + ovr**2*result2d(nd) -&
&                 ovdydd(nd)*ovz - ovdyd*ovzd(nd)
              END DO
              wakeoverlap(turbi, turb, zone) = ovrr*ovrr*result1 + ovr*&
&               ovr*result2 - ovdyd*ovz
            ELSE IF (ovrr .GT. ovr) THEN
              DO nd=1,nbdirs
                wakeoverlapd(nd, turbi, turb, zone) = pi*(ovrd(nd)*ovr+&
&                 ovr*ovrd(nd))
              END DO
              wakeoverlap(turbi, turb, zone) = pi*ovr*ovr
            ELSE
              DO nd=1,nbdirs
                wakeoverlapd(nd, turbi, turb, zone) = pi*(ovrrd(nd)*ovrr&
&                 +ovrr*ovrrd(nd))
              END DO
              wakeoverlap(turbi, turb, zone) = pi*ovrr*ovrr
            END IF
          ELSE
            DO nd=1,nbdirs
              wakeoverlapd(nd, turbi, turb, zone) = 0.0
            END DO
            wakeoverlap(turbi, turb, zone) = 0.0_dp
          END IF
        END DO
      END IF
    END DO
  END DO
  DO turb=1,nturbines
    DO turbi=1,nturbines
      DO nd=1,nbdirs
        wakeoverlapd(nd, turbi, turb, 3) = wakeoverlapd(nd, turbi, turb&
&         , 3) - wakeoverlapd(nd, turbi, turb, 2)
        wakeoverlapd(nd, turbi, turb, 2) = wakeoverlapd(nd, turbi, turb&
&         , 2) - wakeoverlapd(nd, turbi, turb, 1)
      END DO
      wakeoverlap(turbi, turb, 3) = wakeoverlap(turbi, turb, 3) - &
&       wakeoverlap(turbi, turb, 2)
      wakeoverlap(turbi, turb, 2) = wakeoverlap(turbi, turb, 2) - &
&       wakeoverlap(turbi, turb, 1)
    END DO
  END DO
  DO nd=1,nbdirs
    wakeoverlaptrel_matd(nd, :, :, :) = wakeoverlapd(nd, :, :, :)
  END DO
  wakeoverlaptrel_mat = wakeoverlap
  DO turbi=1,nturbines
    DO nd=1,nbdirs
      wakeoverlaptrel_matd(nd, turbi, :, :) = (wakeoverlaptrel_matd(nd, &
&       turbi, :, :)*pi*rotordiameter(turbi)**2/4.0_dp-&
&       wakeoverlaptrel_mat(turbi, :, :)*pi*(rotordiameterd(nd, turbi)*&
&       rotordiameter(turbi)+rotordiameter(turbi)*rotordiameterd(nd, &
&       turbi))/4.0_dp)/(pi*rotordiameter(turbi)*rotordiameter(turbi)/&
&       4.0_dp)**2
    END DO
    wakeoverlaptrel_mat(turbi, :, :) = wakeoverlaptrel_mat(turbi, :, :)/&
&     (pi*rotordiameter(turbi)*rotordiameter(turbi)/4.0_dp)
  END DO
END SUBROUTINE CALCOVERLAPAREAS_DV

!  Differentiation of cttoaxialind in forward (tangent) mode:
!   variations   of useful results: axial_induction
!   with respect to varying inputs: ct
SUBROUTINE CTTOAXIALIND_DV(ct, ctd, nturbines, axial_induction, &
& axial_inductiond, nbdirs)
  
!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
! in
  INTEGER, INTENT(IN) :: nturbines
  REAL(dp), DIMENSION(nturbines), INTENT(IN) :: ct
  REAL(dp), DIMENSION(nbdirs, nturbines), INTENT(IN) :: ctd
! local
  INTEGER :: i
! out
  REAL(dp), DIMENSION(nturbines), INTENT(OUT) :: axial_induction
  REAL(dp), DIMENSION(nbdirs, nturbines), INTENT(OUT) :: &
& axial_inductiond
  INTRINSIC KIND
  INTRINSIC SQRT
  REAL(dp) :: arg1
  REAL(dp), DIMENSION(nbdirs) :: arg1d
  REAL(dp) :: result1
  REAL(dp), DIMENSION(nbdirs) :: result1d
  INTEGER :: nd
  INTEGER :: nbdirs
  axial_induction = 0.0_dp
  DO nd=1,nbdirs
    axial_inductiond(nd, :) = 0.0
  END DO
! execute
  DO i=1,nturbines
    IF (ct(i) .GT. 0.96) THEN
      arg1 = 0.0203_dp - 0.6427_dp*(0.889_dp-ct(i))
      DO nd=1,nbdirs
! Glauert condition
        arg1d(nd) = 0.6427_dp*ctd(nd, i)
        IF (arg1 .EQ. 0.0) THEN
          result1d(nd) = 0.0
        ELSE
          result1d(nd) = arg1d(nd)/(2.0*SQRT(arg1))
        END IF
        axial_inductiond(nd, i) = result1d(nd)
      END DO
      result1 = SQRT(arg1)
      axial_induction(i) = 0.143_dp + result1
    ELSE
      arg1 = 1.0_dp - ct(i)
      DO nd=1,nbdirs
        arg1d(nd) = -ctd(nd, i)
        IF (arg1 .EQ. 0.0) THEN
          result1d(nd) = 0.0
        ELSE
          result1d(nd) = arg1d(nd)/(2.0*SQRT(arg1))
        END IF
        axial_inductiond(nd, i) = -(0.5_dp*result1d(nd))
      END DO
      result1 = SQRT(arg1)
      axial_induction(i) = 0.5_dp*(1.0_dp-result1)
    END IF
  END DO
END SUBROUTINE CTTOAXIALIND_DV


