SUBROUTINE trapzd(func,a,b,s,n)
	USE nrtype; USE nrutil, ONLY : arth
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: a,b
	REAL(SP), INTENT(INOUT) :: s
	INTEGER(I4B), INTENT(IN) :: n
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: func
		END FUNCTION func
	END INTERFACE
	REAL(SP) :: del,fsum
	INTEGER(I4B) :: it
	if (n == 1) then
		s=0.5_sp*(b-a)*sum(func( (/ a,b /) ))
	else
		it=2**(n-2)
		del=(b-a)/it
		fsum=sum(func(arth(a+0.5_sp*del,del,it)))
		s=0.5_sp*(s+del*fsum)
	end if
	END SUBROUTINE trapzd
!**********************************************************************************
	SUBROUTINE polint(xa,ya,x,y,dy)
	USE nrtype; USE nrutil, ONLY : assert_eq,iminloc,nrerror
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: xa,ya
	REAL(SP), INTENT(IN) :: x
	REAL(SP), INTENT(OUT) :: y,dy
	INTEGER(I4B) :: m,n,ns
	REAL(SP), DIMENSION(size(xa)) :: c,d,den,ho
	n=assert_eq(size(xa),size(ya),'polint')
	c=ya
	d=ya
	ho=xa-x
	ns=iminloc(abs(x-xa))
	y=ya(ns)
	ns=ns-1
	do m=1,n-1
		den(1:n-m)=ho(1:n-m)-ho(1+m:n)
		if (any(den(1:n-m) == 0.0)) &
			call nrerror('polint: calculation failure')
		den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
		d(1:n-m)=ho(1+m:n)*den(1:n-m)
		c(1:n-m)=ho(1:n-m)*den(1:n-m)
		if (2*ns < n-m) then
			dy=c(ns+1)
		else
			dy=d(ns)
			ns=ns-1
		end if
		y=y+dy
	end do
	END SUBROUTINE polint

	
!**********************************************************************************	
	FUNCTION qromb(func,a,b)
	USE nrtype; USE nrutil, ONLY : nrerror
	USE nr, ONLY : polint,trapzd
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: a,b
	REAL(SP) :: qromb
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: JMAX=20,JMAXP=JMAX+1,K=5,KM=K-1
	REAL(SP), PARAMETER :: EPS=1.0e-6_sp
	REAL(SP), DIMENSION(JMAXP) :: h,s
	REAL(SP) :: dqromb
	INTEGER(I4B) :: j
	h(1)=1.0
	do j=1,JMAX
		call trapzd(func,a,b,s(j),j)
		if (j >= K) then
			call polint(h(j-KM:j),s(j-KM:j),0.0_sp,qromb,dqromb)
			if (abs(dqromb) <= EPS*abs(qromb)) RETURN
		end if
		s(j+1)=s(j)
		h(j+1)=0.25_sp*h(j)
	end do
	call nrerror('qromb: too many steps')
	END FUNCTION qromb
!****************************************************************
