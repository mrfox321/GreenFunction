!==========================================================
program gunit

use my_functions
use atomic
use superlattice

implicit none

integer :: size_hp,num_hp,num_hp1,is_random,i,j,k_mesh
integer :: nx,ny,om_mesh,k_grid,cnt,num_avg,nn
real :: finish,start
complex(8) :: transmit
type(layers), allocatable :: layer_2(:),layer_3(:)
real(8) :: omega_max,om_imag,om_real,pi
complex(8), allocatable :: omeg(:),self_r(:,:,:),self_l(:,:,:)
real(8), allocatable :: k_pt(:,:)
real(8), allocatable :: transmission(:)
character(len=2) :: k_input
integer :: k_inp



call getarg(1,k_input)

read(k_input,'(i10)') k_inp

print*,'the k-point number =', k_inp

pi=3.14159265358979323846d0

num_avg = 750

nx = 3
ny = 3
size_hp = 5
num_hp = 4
is_random = 0

call layered(layer_2,nx,ny,size_hp,num_hp,is_random)

nn = 3*size(layer_2(1)%a)

k_grid = 10

k_mesh = (k_grid/2+2)*(k_grid/2+1)/2


allocate(k_pt(k_mesh,3))

cnt = 0
do i=1,k_grid/2+1
    do j=1,k_grid/2+1
        if (j.gt.i) then
            exit
        else
            cnt = cnt+1
            k_pt(cnt,1) = (1d0)*(pi)*(real(i-1))/(real(k_grid/2))/real(nx)
            k_pt(cnt,2) = (1d0)*(pi)*(real(j-1))/(real(k_grid/2))/real(ny)
            k_pt(cnt,3) = 0d0
        endif
    enddo
enddo

omega_max = 0.0040d0
om_mesh = 100
om_imag = omega_max*real(1)/real(3000*om_mesh)
allocate(omeg(om_mesh))
omeg = 0d0
do i=1,om_mesh
    om_real = sqrt(omega_max)*real(i)/real(om_mesh)
    omeg(i) = om_real**2 + dcmplx(0,1)*om_imag
enddo

allocate(self_l(om_mesh,nn,nn))
allocate(self_r(om_mesh,nn,nn))

self_l = 0d0
self_r = 0d0

!do i=1,om_mesh
call self_energy(layer_2,k_pt(1,:),omeg(22),size(layer_2),nn,size_hp,self_l(22,:,:),self_r(22,:,:))
!enddo

num_hp1 = 400
is_random = 1

allocate(transmission(num_avg))
transmission = 0d0

do i=1,num_avg
    call layered(layer_3,nx,ny,size_hp,num_hp1,is_random)
call setup_surface(k_pt(1,:),layer_3,size(layer_3),nn,omeg(22),self_l(22,:,:),self_r(22,:,:),transmit)
    transmission(i) = real(transmit)
    print*,real(omeg(22)),transmission(i)
    deallocate(layer_3)
enddo

!do i=1,num_avg
!    transmission(i) = transmission(i)/real(num_avg)
!enddo

open(2333,file='transmission'//k_input//'.dat',status='unknown')
do i=1,num_avg
    write(2333,*)real(omeg(22)),transmission(i)
enddo
close(2333)



end program gunit
!==========================================================





