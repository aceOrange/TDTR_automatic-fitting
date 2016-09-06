subroutine fitting(x0,y0)
use global
use dflib
use dfwin
implicit none
real,intent(in)::x0,y0

	real:: path(-200:200,-200:200)
	real:: dmin,xstep,ystep,tempx,tempy
	real, dimension(300,5):: temparr, temparr1,temparr2,temparr3,temparr4
	integer:: m,n,jj
	
	
	m=0;n=0;path=0;
	tempx=x0; tempy=y0;

	if (tempy < 0.05) then 
		tempy=0.05;
	endif
	call para_init(tempx,tempy);
	call Calculate;
	call distance(path(0,0));
	temparr=output;

		do while (.true.)

		write(*,'(A,F6.2,A1,F7.3,A,F8.3)') 'Now we are checking: ',tempx,', ',tempy,':  ', path(m,n)
!		Calculate the stepsize according to the difference
		if ((abs(m) >= 100) .or. (abs(n) >= 100)) then 
			goto 911
		endif

		if (path(m,n) > 128) then
			xstep=tempx/2.0;
			ystep=tempy/2.0;
			elseif ((path(m,n) <= 128).and.(path(m,n) > 64)) then
				xstep=tempx/4.0;
				ystep=tempy/4.0;
				elseif ((path(m,n) <= 64).and.(path(m,n) > 63)) then
					xstep=tempx/8.0;
					ystep=tempy/8.0;
					elseif ((path(m,n) <= 63).and.(path(m,n) > 8)) then
						xstep=tempx/16.0;
						ystep=tempy/16.0;
						elseif ((path(m,n) <= 8).and.(path(m,n) > 4)) then
							xstep=tempx/32.0;
							ystep=tempy/32.0;
							elseif ((path(m,n) <= 4).and.(path(m,n) > 2)) then
								xstep=tempx/60.0;
								ystep=tempy/60.0;
								elseif ((path(m,n) <= 2).and.(path(m,n) > 0.2)) then
									xstep=tempx/100.0;
									ystep=tempy/100.0;
									elseif ((path(m,n) <= 0.2)) then
										xstep=tempx/150.0;
										ystep=tempy/150.0;

		endif
		if (abs(path(m,n+1)-0) < 1e-6) then
			call para_init(tempx,tempy+ystep);
			call Calculate;
			call distance(path(m,n+1));
			temparr1=output;
		endif
		if ((abs(path(m,n-1)-0) < 1e-8).and.((tempy-ystep) > 0)) then
			call para_init(tempx,tempy-ystep);
			call Calculate
			call distance(path(m,n-1));
			temparr2=output;
		endif
		if (abs(path(m+1,n)-0) < 1e-8) then
			call para_init(tempx+xstep,tempy);
			call Calculate;
			call distance(path(m+1,n));
			temparr3=output;
		endif
		if ((abs(path(m-1,n)-0) < 1e-8) .and. ((tempx-xstep) > 0)) then
			call para_init(tempx-xstep,tempy);
			call Calculate;
			call distance(path(m-1,n));
			temparr4=output;
		endif
		
		dmin=path(m,n);
		do jj=-1,1,2
			if ((dmin > path(m+jj,n)).and.(path(m+jj,n) > 0)) then
				dmin=path(m+jj,n);
			endif 
			if ((dmin > path(m,n+jj)).AND.(path(m,n+jj) > 0)) then
				dmin=path(m,n+jj);
			endif
		end do
		
		if (abs(dmin-path(m,n)) < 1e-8) then
!			write(*,'(A)') '*******************************************************************'
!			write(*,'(A,F8.2,F8.3,F8.3)') 'The final results are: ',tempx,tempy,path(m,n)
			call save_data(tempx,tempy,path(m,n))
			output=temparr;
			call save_generated_data
			goto 911
		else
			do jj=-1 , 1 , 2
				if (abs(dmin - path(m,n+jj)) < 1e-8) then
					n=n+jj;
					tempy=tempy+jj*ystep;
					if (abs(jj+1)<1e-8) then
						temparr=temparr2;
					elseif (abs(jj-1)<1e-8) then
						temparr=temparr1;
					endif
				elseif (abs(dmin - path(m+jj,n)) < 1e-8) then
					m=m+jj;
					tempx=tempx+jj*xstep;
					if (abs(jj+1)<1e-8) then
						temparr=temparr4;
					elseif (abs(jj-1)<1e-8) then
						temparr=temparr3;
					endif
				endif
			end do
		endif
	end do
911 return
end subroutine fitting

!********************************************************************
subroutine distance(tt)
use global
use dflib
use dfwin
implicit none

real,intent(out)::tt

	real::sum
	integer :: index1, index2
	
	index1=1;
	index2=1;
	do while (expm(index1,2) < 100.0)
		index1=index1+1
	end do

	do while (output(index2,1) < 100.0)
		index2=index2+1
	end do
	
	sum=0.0
	do while (output(index2,1) > 0)
		if (abs(expm(index1,2)-output(index2,1)) < 1e-1) then
			sum=sum+(expm(index1,5)-output(index2,4))**2.0;
		endif
		index1=index1+1;
		index2=index2+1;
	end do

	tt=sum
	return
end subroutine distance

!***********************************************************************************************

subroutine save_data(x,y,d)
use dfwin
use dflib
use global
implicit none
real,intent(in)::x,y,d

		open(7, file = 'Result.txt',form='formatted',status='replace')
		write(7, '(f8.3,f8.3,f8.3)') x, y, d
		close(7)


end subroutine save_data
!************************************************************************************************