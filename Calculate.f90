subroutine Calculate
use dflib
use dfwin
use dflogm
use global

integer :: n, nend, jj, ticks, index
complex :: temp, tempin, tempout, tm1, tp1, phasef
real :: freq, omega
real*8 :: time
COMPLEX, ALLOCATABLE, Dimension (nmax) :: tm(:), tp(:)


! print *, 'Calculation started. PLEASE WAIT ...'

ALLOCATE (tm(nmax))
ALLOCATE (tp(nmax))

	   jj = 1
       freq = -probef
       do WHILE (freq < 1.0e12)

			freq = freq+probef

			omega = ABS(freq-pumpf)*pi2
			call gausscalc(omega,tm1)

			omega = (freq+pumpf)*pi2
			call gausscalc(omega,tp1)

!	add a gaussian roll-off to reduce some of the ringing
	
		    tm(jj) = tm1*EXP(-2.0*(freq/1.0e12)**2)				
			tp(jj) = tp1*EXP(-2.0*(freq/1.0e12)**2)

			jj = jj+1
       END do

       tm(1) = CONJG(tm(1))
       nend = jj-1

!	CALCULATE the positive time points

        time = 80.0e-12
		output=0; index=0;
	do WHILE (time < etime)
			index=index+1
			time = (80.0e-12)*((3500.0/80.0)**(index/64.0))
	
			tempin  = (tp(1)+tm(1))/2.0
			tempout = (tp(1)-tm(1))/2.0
			do n = 2, nend
				freq = (n-1)*probef
				phasef = EXP(i*time*freq*pi*2)
				tempin = tempin+(tp(n)+tm(n))*phasef
				tempout = tempout+(tp(n)-tm(n))*phasef
			END do

			temp = REAL(tempin) + i*AIMAG(tempout)
			temp = temp*EXP(i*time*pumpf*pi*2)			
			
			output(index,1) = time*1.0e12;
			output(index,2) = REAL(temp);
			output(index,3) = AIMAG(temp);
			output(index,4) = (-REAL(temp)/AIMAG(temp));
			
!			if (abs(time-100.0e-12)<1e-12) then 
!				print *, 'ratio at 100 ps =', output(index,4) 
!			end if
	end do

DEALLOCATE (tm)
DEALLOCATE (tp)



end subroutine Calculate

!*************************************************************************************
subroutine gausscalc(omega,tt)

USE global
USE nr, ONLY : qromo, midpnt, midinf, qromb
implicit none
!
!              subroutine to calculate the complex temperature of a multilayer
!              slab at frequency omega, full calculation or radial heat flow
!
	interface
		function radialr(m)
		USE nrtype
        implicit none
        REAL(SP), DIMENSION(:), INTENT(IN) :: m
        REAL(SP), DIMENSION(size(m)) :: radialr
        END function radialr
    end interface

    interface
        function radiali(m)
		USE nrtype
        implicit none
        REAL(SP), DIMENSION(:), INTENT(IN) :: m
        REAL(SP), DIMENSION(size(m)) :: radiali
        END function radiali
    end interface

    interface
          function radial(m1)
          implicit none
          REAL, INTENT(IN) :: m1
          COMPLEX :: radial
          END function radial
    end interface


	REAL, INTENT(IN) :: omega
	COMPLEX, INTENT(out) :: tt
	INTEGER :: jj
	REAL :: ttr, tti, cross, k2min

!       Calculate iw/D for each layer at this frequency w

	do jj=1, nlayers
		qn(jj) = i*omega/diff(jj)
	end do
!		check if full integration is needed

    k2min = omega/diffmax
    if (k2min*w0**2 .lt. 1e5 ) then

!		Do the real and imaginary parts of the integration

		cross = 1.5/w0

        ttr = qromb(radialr, 0.0, cross)
        tti = qromb(radiali, 0.0, cross)

        tt = ttr+i*tti
    else
         tt = radial(0.0)*0.5/pisq/w0**2
    END if
        tt = tt*pi2

    return
END subroutine gausscalc

!***************************************************************************


complex function radial(m)
USE global
implicit none

	interface
        function u(j,m2)
        implicit none
        COMPLEX :: u
        integer, INTENT(IN) :: j
        REAL, INTENT(IN) :: m2
		END function u                                               
    end interface

    real, INTENT(IN) :: m
    real::m2

    COMPLEX, DIMENSION(2) :: B, Btemp, Btemp2
    COMPLEX, DIMENSION(2,2) :: alpha, beta=0.0
    INTEGER :: j
    COMPLEX :: gj, gjp1, uj, ujp1

!              initialize the B vector

	B(1) = CMPLX(0.0, 0.0)
    B(2) = CMPLX(1.0, 0.0)

    m2 = m**2

    j = nlayers-1

    do WHILE(j >= 1)

		uj = u(j, m2)
	    ujp1 = u(j+1, m2)

        gjp1 = cond(j+1)*ujp1
        gj = cond(j)*uj
        alpha(1,1) = gj+gjp1
        alpha(1,2) = gj-gjp1
	    alpha(2,1) = alpha(1,2)
        alpha(2,2) = alpha(1,1)

        Btemp = matmul(alpha,B)

        beta(1,1) = EXP(-uj*h(j))
        beta(2,2) = 1.0/beta(1,1)

        Btemp2 = MATMUL(beta,Btemp)

        B = Btemp2/gj
	       
		j = j-1

	END do

	radial = (B(1)+B(2))/(B(2)-B(1))/u(1,m2)/cond(1)

	return
END function radial

!****************************************************************

function u(j,m2)
USE global
implicit none
	INTEGER, INTENT(IN) :: j
	REAL, INTENT(IN) :: m2
	COMPLEX :: u

	u = (qn(j)+4.0*pisq*m2)**0.5		
	return
end function u

!******************************************************************
function radiali(m)
USE global
use nrtype
implicit none

	interface
		function radial(m1)
		implicit none
		REAL, INTENT(IN) :: m1
		COMPLEX :: radial
		END function radial
	end interface

     REAL(SP), DIMENSION(:), INTENT(IN) :: m
     REAL(SP), DIMENSION(size(m)) :: radiali
     REAL :: msc
     INTEGER :: j

     do j = 1, SIZE(m)
		msc = m(j)
        radiali(j) = aimag(radial(msc))*EXP(-pisq*msc**2*w0**2)*msc
     end do

     return

END function radiali
!
!******************************************************************
!

function radialr(m)
use global
use nrtype
implicit none

	interface
		function radial(m1)
		implicit none
		REAL, INTENT(IN) :: m1
		COMPLEX :: radial
		END function radial
	end interface

    REAL(SP), DIMENSION(:), INTENT(IN) :: m
    REAL(SP), DIMENSION(size(m)) :: radialr
    INTEGER :: j
    REAL :: msc

    do j = 1, SIZE(m)
		msc = m(j)
        radialr(j) = real(radial(msc))*EXP(-pisq*msc**2*w0**2)*msc
    end do

    return

END function radialr
      

!*************************************************************************
