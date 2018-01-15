program scat

use my_functions

implicit none

integer :: i,j,chain_length
complex(8) :: eig_v(60,60,12,12),green,t_mat(33,33),dyn(12,12)
real(8) :: eig_f(60,60,12),k_pt(60,60,2),k(2)
real(8) :: R_i(2),R_j(2),energy
real(8) :: phi_odd(3,3),phi_even(3,3),D_chain(30,30)
real(8) :: mass
real(8) :: ifc_r(8,8,1,4,4,3,3)
real :: start, finish
character(len=2) :: k_input
integer :: k_inp

call getarg(1,k_input)

read(k_input,'(i10)') k_inp

!chain_length = 10
!
!call read_spectrum(k_pt,eig_f,eig_v)
!
!open(2333,file='scat'//trim(k_input)//'.dat',status='unknown')
!do i=1,60
!    do j=1,60
!        write(2333,*)i,j,scat_rate(k_pt,eig_f,eig_v,i,j,k_inp,chain_length)
!    enddo
!enddo
!close(2333)

mass = 30.973762d0
k(1) = 0.1243d0
k(2) = 0.3412d0
dyn = dyn_mat(k,mass)

do i=1,12
    print*,dyn(i,i)
enddo

end program scat


