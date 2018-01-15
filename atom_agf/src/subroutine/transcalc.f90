!==========================================================
program gunit



implicit none

integer :: i,j
integer :: k_grid,k_mesh,cnt
integer, allocatable :: sym_factor(:)
character(len=2), dimension(21) :: k_int
real(8) :: trans(21,2,100),avg_tran(2,100)



do i=1,21
    if (i.lt.10) then
        write(k_int(i),'(i1)') i
        !print*,'transmission'//k_int(i)//'.dat'
    else
        write(k_int(i),'(i2)') i
        !print*,'transmission'//k_int(i)//'.dat'
    endif
enddo


k_grid = 10

k_mesh = (k_grid/2+2)*(k_grid/2+1)/2

allocate(sym_factor(k_mesh))

cnt = 0
do i=1,k_grid/2+1
    do j=1,k_grid/2+1
        if (j.gt.i) then
            exit
        else
            cnt = cnt+1

            if (i.eq.j) then
                if (i.eq.1) then
                    sym_factor(cnt) = 1
                elseif (i.eq.(k_grid/2+1)) then
                    sym_factor(cnt) = 1
                else
                    sym_factor(cnt) = 4
                endif
            elseif (j.eq.1) then
                if (i.eq.1) then
                    sym_factor(cnt) = 1
                elseif (i.eq.(k_grid/2+1)) then
                    sym_factor(cnt) = 2
                else
                    sym_factor(cnt) = 4
                endif
            elseif (i.eq.(k_grid/2+1)) then
                if (j.eq.i) then
                    sym_factor(cnt) = 1
                elseif (j.eq.1) then
                    sym_factor(cnt) = 2
                else
                    sym_factor(cnt) = 4
                endif
            else
                sym_factor(cnt) = 8
            endif
        endif
    enddo
enddo


avg_tran = 0d0

do i=1,21
    open(i,file='transmission'//k_int(i)//'.dat',status='old',action='read')
    do j=1,100
        read(i,*) trans(i,1,j),trans(i,2,j)
        if (isnan(trans(i,2,j))) then
            trans(i,2,j) = 0d0
        endif
        avg_tran(2,j) = avg_tran(2,j) + sym_factor(i)*trans(i,2,j)
    enddo
    close(i)
enddo



avg_tran(1,:) = trans(1,1,:)

open(22,file='trans.dat',status='unknown')
do i=1,100
    write(22,*) avg_tran(1,i),avg_tran(2,i)
enddo
close(22)



end program gunit
!==========================================================





