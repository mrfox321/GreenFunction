!======================================
program ifc_output

implicit none

integer :: i,j,k,l,m,n,o,xi,xj,ai,aj,m1,m2,m3
real(8) :: ifc
real(8) :: ifc_r(8,8,1,4,4,3,3)

ifc_r = 0d0


open(9998,file='phosphorene.ifc2',status='old',action='read')
do i=1,8
    read(9998,*)
enddo

do i=1,9360
    if (mod(i,8**2+1).eq.1) then
        read(9998,*) xi,xj,ai,aj
    else
        read(9998,*) m1,m2,m3,ifc
        ifc_r(m1,m2,m3,ai,aj,xi,xj) = ifc
    endif
enddo
close(9998)

open(9997,file='phosphorene.ifc',status='unknown')
do i=1,8
    do j=1,8
        do k=1,1
            do l=1,4
                do m=1,4
                    do n=1,3
                        do o=1,3
                            write(9997,*)i,j,k,l,m,n,o,ifc_r(i,j,k,l,m,n,o)
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo
enddo
close(9997)


end program ifc_output
