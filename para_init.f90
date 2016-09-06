
!*************************************************************************************
subroutine para_init(a,b)
	
use global
implicit none
real,intent(in)::a
real,INTENT(IN)::b
integer:: jj
		
	probef=8.0000e+7; pumpf=9.3e6;

	w0=0.65000e-5;

	nlayers=4;

	cond(1)=19.0e2;			hcap(1)=24.200e6;		h(1)=1.0e-9;
	cond(2)=1.90e2;			hcap(2)=2.4200e6;		h(2)=91e-9;
	cond(3)=b;				hcap(3)=0.1000e6;	    h(3)=1.0e-9;
	cond(4)=a;				hcap(4)=1.64e6;			h(4)=1e-3;
	

	do jj=1,nlayers
		diff(jj)=cond(jj)/hcap(jj);
	end do

	diffmax=0
	do jj=1,nlayers
		if (diffmax < diff(jj)) then
			diffmax=diff(jj)
		endif
	end do
	return
end subroutine para_init