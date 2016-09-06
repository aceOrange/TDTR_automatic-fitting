!*************************************************************************************************
module global
implicit none

	real :: w0
	real, dimension(7) :: cond,hcap,h,diff
	real :: diffmax
	real :: pumpf, probef
	integer::nlayers
	integer::nmax=2e4
	real,parameter :: pi=3.14159, pisq=3.14159**2, pi2=3.14159*2
	real, dimension(300,5) :: output
	REAL, PARAMETER :: etime = 3.3e-9
	
	complex, parameter :: i=(0.0,1.0)
	complex*16, dimension(10) :: qn
	
	real,dimension(300,5)::expm



end module global

!*************************************************************************************************

program Thermo_Calc

use dflib
use dfwin
use global
	real:: x0, y0
	open(3, file='Result.txt', form='formatted', status='old')
		read(3,*) x0, y0
	close(3)

	call load_exp_data
	call fitting(x0,y0);

	
end program Thermo_Calc

!*************************************************************************************************
subroutine load_exp_data
use global
implicit none
	
	real ::sum, v0, vin, vout, thelta
	real, dimension(25):: th
	real,dimension(300,5)::expm2
	integer :: jj

	open(2, file='data', form="formatted", status='old', BLANK='NULL');
	jj=0
	do while (.true.)
		jj=jj+1
		read(2,*,end=911) expm(jj,1), expm(jj,2), expm(jj,3), expm(jj,4), expm(jj,5)
	end do
911 	close(2);

	sum=0;
	do jj=1,19
		sum=sum+expm(jj,4);
	end do
	v0=sum/19;

	do jj=21,40
		vin=expm(jj,3);
		vout=expm(jj,4);
		th(jj-20)=(vin-sqrt(vin**2.0-2.0*vout*(v0-vout)))/vout;
	end do

	sum=0;
	do jj=1,20
		sum=sum+th(jj);
	end do
	thelta=sum/20;

	expm2=expm;
	jj=1;
	do while(expm(jj,5)>0)
		expm2(jj,3)=expm(jj,3)*cos(thelta)-expm(jj,4)*sin(thelta);
		expm2(jj,4)=expm(jj,3)*sin(thelta)+expm(jj,4)*cos(thelta);
		jj=jj+1;
	end do
	expm=expm2;
	return
end subroutine load_exp_data


!*************************************************************************************************
