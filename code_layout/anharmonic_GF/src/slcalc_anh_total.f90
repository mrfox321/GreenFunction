!==========================================================
program gunit

use my_functions
use atomic
use superlattice

implicit none

integer :: size_hp,num_hp,num_hp1,is_random,i,j,k,k_mesh
integer :: nx,ny,om_mesh,k_grid,cnt,num_avg,nn
real :: finish,start
complex(8) :: transmit,trans_a,trans_e
type(layers), allocatable :: layer_2(:),layer_3(:)
real(8) :: omega_max,om_imag,om_real,pi,temperature
complex(8), allocatable :: omeg(:),self_r(:,:,:),self_l(:,:,:),self_a(:,:,:)
real(8), allocatable :: k_pt(:,:)
real(8), allocatable :: transmission(:),transmission_a(:),transmission_e(:)
character(len=2) :: k_input
integer :: k_inp
real(8) :: self_real,self_imag


call getarg(1,k_input)

read(k_input,'(i10)') k_inp

print*,'new'
print*,'the k-point number =', k_inp

pi=3.14159265358979323846d0

num_avg = 1  !number of configations to average over

nx = 3  !unit cell size in x-direction
ny = 3  !unit cell size in y-direction
size_hp = 5 !# of unit cells per layer
num_hp = 4 !# of layers in total system
is_random = 0 !binary variable specifying roughness + nanoparticles

call layered(layer_2,nx,ny,size_hp,num_hp,is_random)  !create "layer" data structure

nn = 3*size(layer_2(1)%a)  !degrees of freedom in a single "layer"

k_grid = 10  !mesh size

k_mesh = (k_grid/2+2)*(k_grid/2+1)/2  !# of k-points in BZ


allocate(k_pt(k_mesh,3))

!constructing the k-mesh
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

!phonon frequencies  (DFT units)
omega_max = 0.0040d0  !upper bound for phonon frequency (need to look at dispersion)
om_mesh = 100  !phonon frequency mesh size
om_imag = omega_max*real(1)/real(3000*om_mesh)  !imaginary phonon frequency (retarded green's function)
allocate(omeg(om_mesh))
omeg = 0d0
do i=1,om_mesh
    om_real = sqrt(omega_max)*real(i)/real(om_mesh)
    omeg(i) = om_real**2 + dcmplx(0,1)*om_imag
enddo

allocate(self_l(om_mesh,nn,nn))
allocate(self_r(om_mesh,nn,nn))
allocate(self_a(om_mesh,3,3))

self_l = 0d0
self_r = 0d0
!self_a has negative imaginary part to preserve sign of retarded green's function

do i=k_inp*5-4,k_inp*5
    call self_energy_sl(layer_2,k_pt(21,:),omeg(i),size(layer_2),nn,size_hp,self_l(i,:,:),self_r(i,:,:))
enddo
open(5000,file='self_en_l21.dat',status='old',action='write')
open(5001,file='self_en_r21.dat',status='old',action='write')
do i=k_inp*5-4,k_inp*5
    do j=1,nn
        do k=1,nn
            write(5000,*) i,j,k,real(self_l(i,j,k)),aimag(self_l(i,j,k))
            write(5001,*) i,j,k,real(self_r(i,j,k)),aimag(self_r(i,j,k))
!             read(5000,*) self_real,self_imag
!             self_l(i,j,k) = self_real + dcmplx(0,1)*self_imag
!             read(5001,*) self_real,self_imag
!             self_r(i,j,k) = self_real + dcmplx(0,1)*self_imag
        enddo
    enddo
enddo
close(5000)
close(5001)


allocate(transmission(om_mesh))
allocate(transmission_a(om_mesh))
allocate(transmission_e(om_mesh))
transmission = 0d0
transmission_a = 0d0
transmission_e = 0d0




open(2445,file='total'//trim(k_input)//'.dat',status='unknown')
do i=1,9  !iteration over length
    if (i.le.5) then
        num_hp1 = 2*i
    elseif (i.gt.5) then
        num_hp1 = (i-4)*10
    endif
    is_random = 1
    call layered(layer_3,nx,ny,size_hp,num_hp1,is_random)
    do k=1,10 !iteration over temperature
        temperature = 10d0*real(k)
        self_a = 0d0
        call self_energy_anh(omeg,om_mesh,temperature,self_a)
        do j=k_inp*5-4,k_inp*5!om_mesh
            call setup_surface(k_pt(21,:),layer_3,size(layer_3),nn,omeg(j),self_l(j,:,:),self_r(j,:,:),self_a(j,:,:)&
            &,transmit,trans_a,k_inp,trans_e)

            write(2445,*),real(omeg(j)),num_hp1,temperature,real(transmit),real(trans_a),real(transmit+trans_a),real(trans_e)
        enddo
    enddo
    deallocate(layer_3)
enddo
close(2445)

end program gunit
!==========================================================





