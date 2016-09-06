subroutine save_generated_data
use dfwin
use dflib
use global
implicit none
		integer :: jj
		jj=1;
		open(6, file = 'TDTRdat.dat',form='formatted',status='replace')
		do while (output(jj,4) > 0)
			write(6, '(f10.3,f10.3,f10.3,f8.3)') output(jj,1),output(jj,2),output(jj,3),output(jj,4)
			jj=jj+1;
		end do
		close(6)


end subroutine save_generated_data
!************************************************************************************************