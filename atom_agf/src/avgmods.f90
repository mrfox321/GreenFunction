module my_functions
implicit none
contains

!==========================================================
subroutine create_full_dyn(d,a,n,nl)

implicit none

integer :: i,j,k,n,nl
complex(8), intent(in) :: d(nl,n,n),a(nl,n,n)
complex(8) :: dyn(nl*n,nl*n)

do i=1,nl*n 
    do j=1,nl*n
        dyn(i,j) = 0d0
    enddo
enddo

!setup diagonal blocks

do i=1,nl
    do j=1,n
        do k=1,n
            dyn((i-1)*n+j,(i-1)*n+k) = d(i,j,k)
        enddo
    enddo
enddo

!setup off-diagonal blocks 

do i=1,nl-1
    do j=1,n
        do k=1,n
            dyn(i*n+j,(i-1)*n+k) = a(i,j,k) !(superdiagonal)
            dyn((i-1)*n+j,i*n+k) = a(i,j,k) !(subdiagonal)
        enddo
    enddo
enddo

open(9001,file='real_g.dat',status='unknown')
open(9002,file='imag_g.dat',status='unknown')
do i=1,nl*n
    do j=1,nl*n
        write(9001,*)i,j,real(dyn(i,j))
        write(9002,*)i,j,imag(dyn(i,j))
    enddo
enddo
close(9001)
close(9002)

end subroutine create_full_dyn
!==========================================================
subroutine setup_superlattice_mass(m,xx,yy,nl,pp)

implicit none

integer :: i,j,k
integer, intent(in) :: xx,yy,nl,pp
real(8) :: m1,m2,m(nl,xx,yy)

m1 = 1d0
m2 = 2d0

do i=1,nl
    do j=1,xx
        do k=1,yy
            if (mod(i-1,2*pp).lt.pp) then
                m(i,j,k) = m1
            else
                m(i,j,k) = m2
            endif
        enddo
    enddo
enddo

end subroutine  setup_superlattice_mass
!==========================================================
subroutine setup_layers(xx,yy,kk,dyn_d,nl,np,dyn_a)

implicit none

integer :: xx,yy,cnt,i,j,k,nl,np
integer :: loc(xx,yy)
complex(8) :: kmat(nl,xx*yy,xx*yy),coupl(nl,xx*yy,xx*yy)
complex(8), intent(out) :: dyn_d(nl,xx*yy,xx*yy), dyn_a(nl,xx*yy,xx*yy)
real(8) :: kk,m(nl,xx,yy)
complex(8) :: sum(nl,xx*yy)

cnt = 0
do i=1,xx
    do j=1,yy
        cnt = cnt + 1
        loc(i,j) = cnt
    enddo
enddo

call setup_superlattice_mass(m,xx,yy,nl,np)

do i=1,nl
    do j=1,xx*yy
        do k=1,xx*yy
            dyn_d(i,j,k) = 0d0
        enddo
    enddo
enddo

do k=1,nl
    do i=1,xx
        do j=1,yy
            if ((i.eq.xx).and.(j.eq.yy)) then
                cycle
            elseif (j.eq.yy) then
                kmat(k,loc(i,j),loc(i+1,j)) = -kk
                kmat(k,loc(i+1,j),loc(i,j)) = -kk

                dyn_d(k,loc(i,j),loc(i+1,j)) = -kk/sqrt(m(k,i,j)*m(k,i+1,j))
                dyn_d(k,loc(i+1,j),loc(i,j)) = -kk/sqrt(m(k,i,j)*m(k,i+1,j))
            elseif (i.eq.xx) then
                kmat(k,loc(i,j),loc(i,j+1)) = -kk
                kmat(k,loc(i,j+1),loc(i,j)) = -kk

                dyn_d(k,loc(i,j),loc(i,j+1)) = -kk/sqrt(m(k,i,j)*m(k,i,j+1))
                dyn_d(k,loc(i,j+1),loc(i,j)) = -kk/sqrt(m(k,i,j)*m(k,i,j+1))
            else
                kmat(k,loc(i,j),loc(i+1,j)) = -kk
                kmat(k,loc(i+1,j),loc(i,j)) = -kk

                dyn_d(k,loc(i,j),loc(i+1,j)) = -kk/sqrt(m(k,i,j)*m(k,i+1,j))
                dyn_d(k,loc(i+1,j),loc(i,j)) = -kk/sqrt(m(k,i,j)*m(k,i+1,j))

                kmat(k,loc(i,j),loc(i,j+1)) = -kk
                kmat(k,loc(i,j+1),loc(i,j)) = -kk

                dyn_d(k,loc(i,j),loc(i,j+1)) = -kk/sqrt(m(k,i,j)*m(k,i,j+1))
                dyn_d(k,loc(i,j+1),loc(i,j)) = -kk/sqrt(m(k,i,j)*m(k,i,j+1))
            endif
        enddo
    enddo
enddo

do i=1,nl
    do j=1,xx*yy
        sum(i,j) = 0
    enddo
enddo



do i=1,nl-1
    do j=1,xx*yy
        coupl(i,j,j) = -kk
    enddo
enddo


do i=1,nl
    do j=1,xx*yy
        do k=1,xx*yy
            sum(i,j) = sum(i,j) + kmat(i,j,k)
        enddo
        
        if (i.eq.1) then
            sum(i,j) = sum(i,j) + coupl(i,j,j)
        elseif (i.eq.nl) then
            sum(i,j) = sum(i,j) + coupl(i-1,j,j)
        else 
            sum(i,j) = sum(i,j) + 2*coupl(i,j,j)
        endif

        kmat(i,j,j) = -sum(i,j)
    enddo
enddo

do k=1,nl
    cnt = 0
    do i=1,xx
        do j=1,yy
            cnt = cnt+1
            dyn_d(k,cnt,cnt) = kmat(k,cnt,cnt)/sqrt(m(k,i,j)*m(k,i,j))
        enddo
    enddo
enddo

do i=1,nl-1
    do j=1,xx*yy
        do k=1,xx*yy
            dyn_a(i,j,k) = 0d0
        enddo
    enddo
enddo

do k=1,nl-1
    cnt = 0
    do i=1,xx
        do j=1,yy
            cnt = cnt+1
            dyn_a(k,cnt,cnt) = coupl(k,cnt,cnt)/sqrt(m(k,i,j)*m(k+1,i,j))
        enddo
    enddo
enddo


end subroutine setup_layers
!==========================================================
!==========================================================
subroutine calculate_delta_sigma(d,a,nl,n,del,sig)

implicit none   

integer :: i,nl,n,j,k
complex(8), intent(in) :: d(nl,n,n),a(nl,n,n)
complex(8), intent(out) :: del(nl,n,n),sig(nl,n,n)
complex(8) :: sigi(n,n),deli(n,n)

del(1,:,:) = d(1,:,:)
sig(nl,:,:) = d(nl,:,:)

do i=2,nl
call z_inv(del(i-1,:,:),deli,n)
call z_inv(sig(nl-i+2,:,:),sigi,n)
del(i,:,:) = d(i,:,:) - matmul_z(a(i-1,:,:),matmul_z(deli,transpose(conjg(a(i-1,:,:))),n),n)
sig(nl-i+1,:,:) = d(nl-i+1,:,:) - matmul_z(transpose(conjg(a(nl-i+1,:,:))),matmul_z(sigi,a(nl-i+1,:,:),n),n)
enddo

!do i=1,nl
!print*, sig(i,40,40)
!enddo

end subroutine calculate_delta_sigma
!==========================================================
subroutine calculate_u_v(d,a,nl,n,uv)

implicit none

integer :: i,j,k
integer, intent(in) :: nl,n
complex(8), intent(in) :: d(nl,n,n), a(nl,n,n)
complex(8), intent(out) :: uv(nl,n,n)
complex(8) :: del(nl,n,n), sig(nl,n,n),u(n,n),v(nl,n,n)
complex(8) :: v_c(nl,n,n), sigi_c(n,n)

call calculate_delta_sigma(d,a,nl,n,del,sig)

u = del(1,:,:)

call z_inv(sig(1,:,:),sigi_c,n)

v(1,:,:) = sigi_c

do i=2,nl
call z_inv(sig(i,:,:),sigi_c,n)
v(i,:,:) = matmul_z(v(i-1,:,:),matmul_z(transpose(conjg(a(i-1,:,:))),sigi_c,n),n)
enddo

do i=1,nl
    uv(i,:,:) = matmul_z(u,v(i,:,:),n)
enddo

end subroutine calculate_u_v    
!==========================================================
subroutine z_inv(mat_a,mat_b,nc)

implicit none

integer :: info, nc,i
integer :: ipiv(nc)
complex(8) :: work(nc)
complex(8), intent(in) :: mat_a(nc,nc)
complex(8), intent(out) :: mat_b(nc,nc)
!complex(8) :: mat_c(nc,nc)

mat_b = mat_a

call zgetrf(nc,nc,mat_b,nc,ipiv,info)
call zgetri(nc,mat_b,nc,ipiv,work,nc,info)

end subroutine z_inv
!==========================================================

function matmul_z(a,b,n) result(c)


integer, intent(in) :: n
complex(8) :: a(n,n),b(n,n)
complex(8) :: c(n,n), alpha, beta

alpha = cmplx(1d0,0d0)
beta = 0d0

call zgemm('N','N',n,n,n,alpha,a,n,b,n,beta,c,n)

end function matmul_z
!==========================================================
subroutine green_iter(H_d,H_a,omega,n,G_s)

implicit none

integer :: itt,i,n
complex(8), intent(in) :: H_d(n,n),H_a(n,n),omega(n,n)
complex(8), intent(out) :: G_s(n,n)
complex(8) :: T(n,n),inv(n,n)

itt = 10000

!initialize

call z_inv(omega-H_d,inv,n)

T = matmul_z(inv,transpose(conjg(H_a)),n)

do i=1,itt
    call z_inv(omega-H_d-matmul_z(H_a,T,n),inv,n)
    T = matmul_z(inv,transpose(conjg(H_a)),n)
enddo

call z_inv(omega-H_d-matmul_z(H_a,T,n),G_s,n)
    

end subroutine green_iter
!==========================================================
subroutine lopez_green(H_d,H_a,omega,n,G_s)

implicit none

integer :: itt,i,n
complex(8) :: t_a(n,n),t_b(n,n),s_a(n,n),s_b(n,n)
complex(8) :: T(n,n),r(n,n),inv(n,n),id(n,n),G_inv(n,n),prod(n,n)
complex(8), intent(in) :: H_d(n,n),H_a(n,n),omega(n,n)
complex(8), intent(out) :: G_s(n,n)


itt = 300

id = 0d0

do i=1,n
    id(i,i) = 1d0
enddo

call z_inv(omega-H_d,inv,n)

t_a = matmul_z(inv,transpose(conjg(H_a)),n)
t_b = matmul_z(inv,H_a,n)

T = t_a
r = t_b

do i=1,itt
    prod = id - matmul_z(t_a,t_b,n) - matmul_z(t_b,t_a,n)
    
    call z_inv(prod,inv,n)

    s_a = matmul_z(inv,matmul_z(t_a,t_a,n),n)
    s_b = matmul_z(inv,matmul_z(t_b,t_b,n),n)

    T = T + matmul_z(r,s_a,n)

    r = matmul_z(r,s_b,n)

    t_a = s_a
    t_b = s_b
enddo

G_inv = omega - H_d - matmul_z(H_a,T,n)

call z_inv(G_inv,G_s,n)



end subroutine lopez_green
!==========================================================
subroutine green_surface(omega,H_lr,green_s,tau,g_surface,n)

implicit none

integer, intent(in) :: n
complex(8) :: prod(n,n)
complex(8), intent(in) :: H_lr(n,n),green_s(n,n),tau(n,n),omega(n,n)
complex(8), intent(out) :: g_surface(n,n)


prod = omega - H_lr - matmul_z(tau,matmul_z(green_s,transpose(conjg(tau)),n),n)

call z_inv(prod,g_surface,n)


end subroutine green_surface
!==========================================================
end module my_functions

!======================================
module atomic

implicit none
!======================================
type atoms
    integer :: c(3),n(46,2),t !c is coordinate of atom, n are the neighbors (layer #, atom #)
    real(8) :: m
end type atoms
!======================================
type layers 
    type(atoms), allocatable :: a(:)
end type layers
!======================================
end module atomic
!======================================
module superlattice
implicit none
contains
!======================================
subroutine layered(layer,xx,yy,size_hper,num_hper,is_random)

use atomic
implicit none

integer, intent(in) :: xx,yy,size_hper,num_hper,is_random
integer :: rl(8,3),blocks,nk,nl,is_periodic
integer :: i,j,k,l,m,a,b,c,aa,bb,cc,n(2,46,3)
integer, allocatable :: find_atom(:,:,:,:)
type(layers), intent(inout), allocatable :: layer(:)
!complex(8), intent(out), allocatable :: dynmat_d(:,:,:),dynmat_a(:,:,:)
!real(8), intent(in) :: kpt(3)

rl(1,1) = 0 
rl(1,2) = 0
rl(1,3) = 0
rl(2,1) = 0
rl(2,2) = 2
rl(2,3) = 2
rl(3,1) = 2
rl(3,2) = 0
rl(3,3) = 2
rl(4,1) = 2
rl(4,2) = 2
rl(4,3) = 0
rl(5,1) = 3
rl(5,2) = 1
rl(5,3) = 1
rl(6,1) = 1
rl(6,2) = 3
rl(6,3) = 1
rl(7,1) = 1
rl(7,2) = 1
rl(7,3) = 3
rl(8,1) = 3
rl(8,2) = 3
rl(8,3) = 3

blocks = size_hper*num_hper 

allocate(layer(blocks))
allocate(find_atom(0:xx*4-1,0:yy*4-1,0:blocks*4-1,2))


find_atom = 0

do m=1,2
nk = 0
nl = 0
    do i=1,blocks
        nk=0
        nl=nl+1
        do j=1,xx
            do k=1,yy
                do l=1,8
                    nk=nk+1
                    if (m.eq.2) then
                        layer(nl)%a(nk)%c(1)= rl(l,1)+4*(j-1)
                        layer(nl)%a(nk)%c(2)= rl(l,2)+4*(k-1)
                        layer(nl)%a(nk)%c(3)= rl(l,3)+4*(i-1)

                        a = rl(l,1)+4*(j-1)
                        b = rl(l,2)+4*(k-1)
                        c = rl(l,3)+4*(i-1)

                        layer(nl)%a(nk)%t = mod(a+b+c+1,2)+1

                        find_atom(a,b,c,1) = nl
                        find_atom(a,b,c,2) = nk
                    endif
                enddo
            enddo
        enddo
        if (m.eq.1) then
            allocate(layer(nl)%a(nk))
        endif
    enddo
enddo


!open(2001,file='coords.dat',status='unknown')
!do i=1,size(layer)
!    do j=1,size(layer(i)%a)
!        write(2001,*)layer(i)%a(j)%c(1),layer(i)%a(j)%c(2),layer(i)%a(j)%c(3)
!    enddo
!enddo
!close(2001)

call read_neighbors(n)

is_periodic = 0

do i=1,size(layer)
    do j=1,size(layer(i)%a)
        do k=1,46
            a = layer(i)%a(j)%c(1)
            b = layer(i)%a(j)%c(2)
            c = layer(i)%a(j)%c(3)
                        
            aa = modulo(a+n(layer(i)%a(j)%t,k,1),xx*4)
            bb = modulo(b+n(layer(i)%a(j)%t,k,2),yy*4)
            if (is_periodic.eq.0) then
                cc = c+n(layer(i)%a(j)%t,k,3)
            elseif (is_periodic.eq.1) then
                cc = modulo(c+n(layer(i)%a(j)%t,k,3),4*blocks)
            endif
    
            if ((0.le.cc) .and. (cc.le.(4*blocks-1))) then
                if (find_atom(aa,bb,cc,2).eq.0) then
                   print*,'fucked',aa,bb,cc
                endif
                layer(i)%a(j)%n(k,1) = find_atom(aa,bb,cc,1)
                layer(i)%a(j)%n(k,2) = find_atom(aa,bb,cc,2)
            else
                layer(i)%a(j)%n(k,1) = 0
                layer(i)%a(j)%n(k,2) = 0
            endif
        enddo
    enddo
enddo

call setup_mass(layer,size_hper,blocks,is_random)


end subroutine layered
!======================================
subroutine read_neighbors(n)

implicit none

integer :: i,j
integer, intent(out) :: n(2,46,3)

open(9999,file='neighbors.dat',status='old',action='read')

do i=1,46
    read(9999,*)n(1,i,1),n(1,i,2),n(1,i,3)
enddo
close(9999)

do i=1,46
    do j=1,3
        n(2,i,j) = -n(1,i,j)
    enddo
enddo

end subroutine read_neighbors
!======================================
subroutine setup_mass(layer,hper,bloks,is_random)

use atomic
implicit none

real(8) :: m1,m2,m3,m4
integer :: i,j,k,cnt
integer, intent(in) :: hper,bloks,is_random
type(layers),intent(inout) :: layer(bloks)
integer :: mm,clock
integer, allocatable :: seed(:)
real :: y
real(8) :: p



m1 = 69.723d0  !(gallium)
m2 = 26.9827d0 !(aluminum)
m3 = 74.9216d0  !(arsenic)
m4 = 167.259d0 !(erbium)

!m1 = 1d0 !(silicon)
!m2 = 2.4825d0 !(germanium)
!m3 = 0.4285714285d0 !(carbon)

do i=1,bloks
    do j=1,size(layer(i)%a)
        if (mod(i-1,2*hper).lt.hper) then
            if (layer(i)%a(j)%t.eq.1) then
                layer(i)%a(j)%m = m1  !gallium
            elseif (layer(i)%a(j)%t.eq.2) then
                layer(i)%a(j)%m = m3 !arsenic
            endif
        else
            if (layer(i)%a(j)%t.eq.1) then
                layer(i)%a(j)%m = m2 !aluminum
            elseif (layer(i)%a(j)%t.eq.2) then
                layer(i)%a(j)%m = m3 !arsenic
            endif
        endif
    enddo
enddo

call random_seed(size=mm)
allocate(seed(mm))
CALL SYSTEM_CLOCK(COUNT=clock)
seed = clock + 37 * (/ (i - 1, i = 1, mm) /)
!seed  = 1367986020 + 37 * (/ (i - 1, i = 1, mm) /)
!print*,'the CLOCK VALUE is', clock
!print*,'the SEED IS',seed
CALL RANDOM_SEED(PUT = seed)
deallocate(seed)


p = 0.5d0
if (is_random.eq.1) then
    do i=2,bloks-1
        do j=1,size(layer(i)%a)
            if (layer(i)%a(j)%t.eq.1) then
                if (mod(i-1,2*hper).eq.0)then
                    call random_number(y)
                    if (y.lt.p) then
                        layer(i)%a(j)%m = m2
                    endif
                elseif (mod(i-1,2*hper).eq.(hper-1)) then
                    call random_number(y)
                    if (y.lt.p) then   
                        layer(i)%a(j)%m = m2
                    endif
                elseif (mod(i-1,2*hper).eq.(hper)) then
                    call random_number(y)
                    if (y.lt.p) then
                        layer(i)%a(j)%m = m1
                    endif
                elseif (mod(i-1,2*hper).eq.(2*hper-1)) then
                    call random_number(y)
                    if (y.lt.p) then
                        layer(i)%a(j)%m = m1
                    endif
                endif
            endif
        enddo
    enddo
endif



call random_seed(size=mm)
allocate(seed(mm))
CALL SYSTEM_CLOCK(COUNT=clock)
seed = clock + 37 * (/ (i - 1, i = 1, mm) /)
!seed  = 1729137069 + 37 * (/ (i - 1, i = 1, mm) /)
!print*,'the CLOCK VALUE is', clock
!print*,'the SEED IS',seed
CALL RANDOM_SEED(PUT = seed)
deallocate(seed)

!cnt = 0
!p = 0.0016d0
!if (is_random.eq.1) then
!    do i=3,bloks-2
!        do j=1,size(layer(i)%a)
!            call random_number(y)
!            if (y.lt.p) then
!                cnt = cnt+1
!                print*,'nanoparticle',i,j,cnt
!                if (layer(i)%a(j)%t.eq.1) then
!                    layer(i)%a(j)%m = m4
!                endif
!                do k=1,46
!                    if (layer(layer(i)%a(j)%n(k,1))%a(layer(i)%a(j)%n(k,2))%t.eq.1) then
!                        layer(layer(i)%a(j)%n(k,1))%a(layer(i)%a(j)%n(k,2))%m = m4
!                    endif
!                enddo
!            endif
!        enddo
!    enddo
!endif


!p = 0.004d0
p = 0.0004d0
cnt = 0
if (is_random.eq.1) then
    do i=3,bloks-2
        do j=1,size(layer(i)%a)
            if (mod(i-1,2*hper).eq.0)then
                call random_number(y)
                if (y.lt.p) then
                    cnt = cnt+1
                    print*,'nanoparticle',i,j,cnt
                    if (layer(i)%a(j)%t.eq.1) then
                        layer(i)%a(j)%m = m4
                    endif
                    do k=1,46
                        if (layer(layer(i)%a(j)%n(k,1))%a(layer(i)%a(j)%n(k,2))%t.eq.1) then
                            layer(layer(i)%a(j)%n(k,1))%a(layer(i)%a(j)%n(k,2))%m = m4
                        endif
                    enddo
                endif
            elseif (mod(i-1,2*hper).eq.(hper-1)) then
                call random_number(y)
                if (y.lt.p) then
                    cnt = cnt+1
                    print*,'nanoparticle',i,j,cnt
                    if (layer(i)%a(j)%t.eq.1) then
                        layer(i)%a(j)%m = m4
                    endif
                    do k=1,46
                        if (layer(layer(i)%a(j)%n(k,1))%a(layer(i)%a(j)%n(k,2))%t.eq.1) then
                            layer(layer(i)%a(j)%n(k,1))%a(layer(i)%a(j)%n(k,2))%m = m4
                        endif
                    enddo
                endif
            elseif (mod(i-1,2*hper).eq.(hper)) then
                call random_number(y)
                if (y.lt.p) then
                    cnt = cnt+1
                    print*,'nanoparticle',i,j,cnt
                    if (layer(i)%a(j)%t.eq.1) then
                        layer(i)%a(j)%m = m4
                    endif
                    do k=1,46
                        if (layer(layer(i)%a(j)%n(k,1))%a(layer(i)%a(j)%n(k,2))%t.eq.1) then
                            layer(layer(i)%a(j)%n(k,1))%a(layer(i)%a(j)%n(k,2))%m = m4
                        endif
                    enddo
                endif
            elseif (mod(i-1,2*hper).eq.(2*hper-1)) then
                call random_number(y)
                if (y.lt.p) then
                    cnt = cnt+1
                    print*,'nanoparticle',i,j,cnt
                    if (layer(i)%a(j)%t.eq.1) then
                        layer(i)%a(j)%m = m4
                    endif
                    do k=1,46
                        if (layer(layer(i)%a(j)%n(k,1))%a(layer(i)%a(j)%n(k,2))%t.eq.1) then
                            layer(layer(i)%a(j)%n(k,1))%a(layer(i)%a(j)%n(k,2))%m = m4
                        endif
                    enddo
                endif
            endif
        enddo
    enddo
endif



end subroutine setup_mass
!======================================
subroutine setup_ifc(kpt,layer,bloks,num_atom,dyn_d,dyn_a)

use atomic
implicit none

integer, intent(in) :: bloks,num_atom
type(layers), intent(in) :: layer(bloks)
integer :: i,j,k,l,m,n,ref,num_1,num_2,cnt,n_vec(2,46,3)
real(8) :: ifc_d(3*num_atom,3*num_atom),ifc_a(3*num_atom,3*num_atom)
real(8) :: ifc(2,46,3,3),sum_f(3*num_atom)
complex(8), intent(inout) :: dyn_d(bloks,3*num_atom,3*num_atom), dyn_a(bloks-1,3*num_atom,3*num_atom)
complex(8) :: ph
real(8) :: pi,rk1,rk2,rk3,m_i,m_j
real(8), intent(in) :: kpt(3)


pi=3.14159265358979323846d0

ifc_d = 0d0
ifc_a = 0d0
dyn_d = 0d0
dyn_a = 0d0

ref = 3


!call ifc_si(ifc)
call ifc_gaas(ifc)

call read_neighbors(n_vec)


do i=1,num_atom
    do j=1,46
        if (layer(ref)%a(i)%n(j,1).eq.ref) then
            l=layer(ref)%a(i)%n(j,2)

            !if (i.eq.l) then
            !    print*, 'overlap',i,l
            !endif

            do m=1,3
                do n=1,3
                    ifc_d(3*(i-1)+m,3*(l-1)+n) = ifc(layer(ref)%a(i)%t,j,m,n)+ifc_d(3*(i-1)+m,3*(l-1)+n)
                enddo
            enddo
        elseif (layer(ref)%a(i)%n(j,1).eq.(ref+1)) then
            l = layer(ref)%a(i)%n(j,2)
            do m=1,3
                do n=1,3
                    ifc_a(3*(i-1)+m,3*(l-1)+n) = ifc(layer(ref)%a(i)%t,j,m,n)+ifc_a(3*(i-1)+m,3*(l-1)+n)
                enddo
            enddo
        endif
    enddo
enddo


!diagonal components of ifc_d
do i=1,3*num_atom
    sum_f(i) = 0
    do j=1,3*num_atom
        sum_f(i) = sum_f(i) + ifc_d(i,j) + ifc_a(i,j) + ifc_a(j,i)
    enddo
    ifc_d(i,i) = ifc_d(i,i)-sum_f(i)
enddo

!do i=1,3*num_atom
!    print*,ifc_d(i,i),-sum_f(i)
!enddo

do i=1,bloks
    do j=1,num_atom
        do k=1,46
            if (layer(i)%a(j)%n(k,1).ne.0) then
                num_1 = layer(i)%a(j)%n(k,2)
                rk1 = kpt(1)*n_vec(layer(i)%a(j)%t,k,1)/4d0  
                rk2 = kpt(2)*n_vec(layer(i)%a(j)%t,k,2)/4d0 
                rk3 = kpt(3)*n_vec(layer(i)%a(j)%t,k,3)/4d0 
                ph = dcmplx(0,1)*(rk1+rk2+rk3)
                m_i = layer(i)%a(j)%m
                m_j = layer(layer(i)%a(j)%n(k,1))%a(layer(i)%a(j)%n(k,2))%m
                do m=1,3
                    do n=1,3
                        if (i.eq.layer(i)%a(j)%n(k,1)) then
dyn_d(i,3*(j-1)+m,3*(num_1-1)+n) = ifc(layer(i)%a(j)%t,k,m,n)*zexp(ph)/sqrt(m_i*m_j)+dyn_d(i,3*(j-1)+m,3*(num_1-1)+n)
                        elseif ((i+1).eq.layer(i)%a(j)%n(k,1).and.(i.ne.bloks)) then
dyn_a(i,3*(j-1)+m,3*(num_1-1)+n) = ifc(layer(i)%a(j)%t,k,m,n)*zexp(ph)/sqrt(m_i*m_j)+dyn_a(i,3*(j-1)+m,3*(num_1-1)+n)
                        endif
                    enddo
                enddo
            endif
        enddo
    enddo
enddo
do i=1,bloks
    cnt = 0
    do j=1,num_atom
        m_i = layer(i)%a(j)%m
        do m=1,3
            cnt = cnt+1
            dyn_d(i,cnt,cnt) = dyn_d(i,cnt,cnt)-sum_f(cnt)/m_i
        enddo
    enddo
enddo

do i=1,bloks-1
    do j=1,num_atom
        !print*,i,dyn_a(i,3*(j-1)+1,3*(j-1)+1)
    enddo
enddo



end subroutine setup_ifc
!======================================
subroutine self_energy(layer_1,kpt,omeg,bloks,n,size_hper,self_l,self_r)

use my_functions
use atomic

implicit none

integer :: i
integer, intent(in) :: size_hper,n,bloks
type(layers), intent(in) :: layer_1(bloks)
complex(8) :: dyn_d(bloks,n,n), dyn_a(bloks-1,n,n),omega(n,n)
complex(8) :: g_r(n,n),g_l(n,n)
real(8), intent(in) :: kpt(3)
complex(8), intent(in) :: omeg
complex(8), intent(out) :: self_l(n,n),self_r(n,n)

dyn_d = 0d0
dyn_a = 0d0

call setup_ifc(kpt,layer_1,size(layer_1),size(layer_1(1)%a),dyn_d,dyn_a)

g_r = 0d0
g_l = 0d0

self_l = 0d0
self_r = 0d0

omega = 0d0

do i=1,n
    omega(i,i) = omeg
enddo


call lopez_green(dyn_d(1,:,:),dyn_a(1,:,:),omega,n,g_r)
call lopez_green(dyn_d(1,:,:),transpose(conjg(dyn_a(1,:,:))),omega,n,g_l)

!comment for homogenous system
call z_inv(omega-dyn_d(2,:,:)-matmul_z(dyn_a(2,:,:),matmul_z(g_r,transpose(conjg(dyn_a(2,:,:))),n),n),g_r,n)
call z_inv(omega-dyn_d(2,:,:)-matmul_z(transpose(conjg(dyn_a(2,:,:))),matmul_z(g_l,dyn_a(2,:,:),n),n),g_l,n)

!GaAs-GaAs
self_l = matmul_z(transpose(conjg(dyn_a(2,:,:))),matmul_z(g_l,dyn_a(2,:,:),n),n)
!AlAs-GaAs
self_r = matmul_z(dyn_a(2*size_hper,:,:),matmul_z(g_r,transpose(conjg(dyn_a(2*size_hper,:,:))),n),n)

end subroutine self_energy
!======================================
subroutine setup_surface(kpt,layer_2,bloks1,n,omegas,self_l,self_r,transmission)

use my_functions
use atomic
implicit none

integer, intent(in) :: bloks1,n
complex(8), intent(in) :: self_l(n,n),self_r(n,n),omegas
real(8), intent(in) :: kpt(3)
type(layers), intent(in) :: layer_2(bloks1)
integer :: i,j,k,l
complex(8) :: omega_1(n,n)
complex(8) :: dyn_d_1(bloks1,n,n),dyn_a_1(bloks1-1,n,n)
complex(8) :: g_1n(n,n),gam_l(n,n),gam_r(n,n),trans(n,n)
complex(8), intent(out) :: transmission

dyn_d_1 = 0d0
dyn_a_1 = 0d0

call setup_ifc(kpt,layer_2,size(layer_2),size(layer_2(1)%a),dyn_d_1,dyn_a_1)


omega_1 = 0d0

do i=1,n
    omega_1(i,i) = omegas
enddo

g_1n = 0d0

call green_n(real(omega_1),dyn_d_1,dyn_a_1,self_l,self_r,g_1n,size(dyn_d_1(:,1,1)),n)

gam_l = 0d0
gam_r = 0d0

gam_l = self_l - transpose(conjg(self_l))
gam_r = self_r - transpose(conjg(self_r))

do i=1,n
    do j=1,n
        gam_l(i,j) = dcmplx(0,1)*gam_l(i,j)
        gam_r(i,j) = dcmplx(0,1)*gam_r(i,j)
    enddo
enddo

trans = 0d0
trans = matmul_z(gam_l,matmul_z(g_1n,matmul_z(gam_r,transpose(conjg(g_1n)),n),n),n)

transmission = 0d0

do i=1,n
    transmission = transmission + trans(i,i)
enddo

print*, transmission
 
end subroutine setup_surface
!======================================
subroutine green_inverse(omega,dyn_d,dyn_a,self_l,self_r,g_1n,bloks,n)

use atomic
use my_functions
implicit none

integer, intent(in) :: bloks,n
integer :: i,j,k,l
real(8), intent(in) :: omega
!real(8) :: omega_1(n*bloks,n*bloks)
complex(8), intent(in) :: dyn_d(bloks,n,n),dyn_a(bloks-1,n,n)
complex(8) :: g_gg(n*bloks,n*bloks)
complex(8), intent(in) :: self_l(n,n),self_r(n,n)
complex(8), intent(out) :: g_1n(n,n)
complex(8) :: dyn(bloks*n,bloks*n)


dyn = 0d0
g_1n = 0d0
g_gg = 0d0

do j=1,bloks-1
    do k=1,n  
        do l=1,n            
            dyn(n*(j-1)+k,n*(j-1)+n+l) = dyn_a(j,k,l)
        enddo
    enddo
enddo
dyn = dyn+transpose(conjg(dyn))
do j=1,bloks
    do k=1,n
        do l=1,n
            dyn(n*(j-1)+k,n*(j-1)+l) = dyn_d(j,k,l)
            if (j.eq.1) then
                dyn(n*(j-1)+k,n*(j-1)+l) = dyn(n*(j-1)+k,n*(j-1)+l) + self_l(k,l)
            elseif (j.eq.bloks) then
                dyn(n*(j-1)+k,n*(j-1)+l) = dyn(n*(j-1)+k,n*(j-1)+l) + self_r(k,l)
            endif
        enddo
    enddo
enddo

!omega_1 = 0d0

!do j=1,n*bloks
!    omega_1(j,j) = omega(1,1)
!enddo

!call z_inv(omega_1-dyn,g_gg,n*bloks)

do j=1,n*bloks
    do k=1,n*bloks
        if (j.eq.k) then
            dyn(j,k) = omega - dyn(j,k)
        else
            dyn(j,k) = -dyn(j,k)
        endif
    enddo
enddo

call z_inv(dyn,g_gg,n*bloks)

!open(5500,file='dyn_mat.dat',status='unknown')
!do j=1,n*bloks
!    do k=1,n*bloks
!        write(5500,*)j,k,real(dyn(j,k)),imag(dyn(j,k))
!    enddo
!enddo
!close(5500)
!
!open(5501,file='g_gg.dat',status='unknown')
!do j=1,n*bloks
!    do k=1,n*bloks
!        write(5501,*)j,k,real(g_gg(j,k)),imag(g_gg(j,k))
!    enddo
!enddo
!close(5501)


do j=1,n
    do k=1,n
        g_1n(j,k) = g_gg(j,n*(bloks-1)+k)
    enddo
enddo

!open(5555,file='g_f.dat',status='unknown')
!do j=1,n
!    do k=1,n
!        write(5555,*)j,k,real(g_1n(j,k)),imag(g_1n(j,k))
!    enddo
!enddo
!close(5555)


end subroutine green_inverse
!======================================
subroutine green_n(omega,dyn_d,dyn_a,self_l,self_r,g_1n,bloks,n)

use my_functions
implicit none

integer, intent(in) :: bloks,n
integer :: i,j,k
real(8), intent(in) :: omega(n,n)
complex(8), intent(in) :: dyn_d(bloks,n,n),dyn_a(bloks-1,n,n)
complex(8) :: g_nn(n,n),g_gg(n,n),prod(n,n),prod1(n,n)
complex(8), intent(in) :: self_l(n,n),self_r(n,n)
complex(8), intent(out) :: g_1n(n,n)

g_nn = 0d0
g_gg = 0d0

!open(2002,file='dyn_a.dat',status='unknown')
!do i=1,n
!    do j=1,n
!        do k=1,bloks-1
!            write(2002,*)i,j,k,real(dyn_a(k,i,j))
!        enddo
!    enddo
!enddo
!close(2002)
!
!open(2003,file='dyn_d.dat',status='unknown')
!do i=1,n
!    do j=1,n
!        do k=1,bloks
!            write(2003,*)i,j,k,real(dyn_d(k,i,j))
!        enddo
!    enddo
!enddo
!close(2003)
!
!open(5502,file='self_l.dat',status='unknown')
!do j=1,n
!    do k=1,n
!        write(5502,*)j,k,real(self_l(j,k)),imag(self_l(j,k))
!    enddo
!enddo
!close(5502)
!
!open(5503,file='self_r.dat',status='unknown')
!do j=1,n
!    do k=1,n
!        write(5503,*)j,k,real(self_r(j,k)),imag(self_r(j,k))
!    enddo
!enddo
!close(5503)

call z_inv(omega-dyn_d(1,:,:)-self_l,g_nn,n)
g_1n = g_nn

do i=2,bloks-1
    call z_inv(omega-dyn_d(i,:,:)-matmul_z(transpose(conjg(dyn_a(i-1,:,:))),matmul_z(g_nn,dyn_a(i-1,:,:),n),n),g_nn,n)
    g_1n = matmul_z(g_1n,matmul_z(dyn_a(i-1,:,:),g_nn,n),n)
enddo

g_gg = omega-dyn_d(bloks,:,:)-self_r-matmul_z(transpose(conjg(dyn_a(bloks-1,:,:))),matmul_z(g_nn,dyn_a(bloks-1,:,:),n),n)
call z_inv(g_gg,g_nn,n)
g_1n = matmul_z(g_1n,matmul_z(dyn_a(bloks-1,:,:),g_nn,n),n)

!open(5556,file='g_n.dat',status='unknown')
!do j=1,n
!    do k=1,n
!        write(5556,*)j,k,real(g_1n(j,k)),imag(g_1n(j,k))
!    enddo
!enddo
!close(5556)


end subroutine green_n
!======================================
subroutine ifc_si(mat_dof)

implicit none

integer :: x,y
real(8) :: al,be,mu(2:5),nu(2:5),de(2:5),la(2:5)
real(8) :: a(3,3),b(3,3),c(3,3),d(3,3)
real(8) :: e(2:5,3,3),f(2:5,3,3),g(2:5,3,3),h(2:5,3,3),i(2:5,3,3),j(2:5,3,3)
real(8) :: k(2:5,3,3),l(2:5,3,3),m(2:5,3,3),n(2:5,3,3),o(2:5,3,3),p(2:5,3,3)
real(8), intent(out) :: mat_dof(2,46,3,3)

al = -33.14d0
be = -23.17d0

mu(2) = -1.90d0
mu(3) = 0.37d0
mu(4) = -0.22d0
mu(5) = -0.17d0

nu(2) = -1.84d0
nu(3) = -0.37d0
nu(4) = 0d0
nu(5) = -0.26d0

de(2) = 1.10d0
de(3) = 0.29d0
de(4) = 0d0 
de(5) = 0.56d0

la(2) = 4.22d0
la(3) = 0.09d0
la(4) = -0.05d0
la(5) = -1.79d0

do x=1,3
    do y=1,3
        if (x.eq.y) then
            a(x,y) = al
        else 
            a(x,y) = be
        endif
    enddo
enddo

do x=1,3
    do y=1,3
        if (x.eq.y) then
            b(x,y) = al
        elseif ((x+y).eq.5) then !picks out i,j = 2,3 or 3,2
            b(x,y) = be
        else    
            b(x,y) = -be
        endif
    enddo
enddo

do x=1,3
    do y=1,3
        if (x.eq.y) then
            c(x,y) = al
        elseif (((x+y).eq.4).and.(x.ne.y)) then !picks out 1,3 or 3,1
            c(x,y) = be
        else
            c(x,y) = -be
        endif
    enddo
enddo

do x=1,3
    do y=1,3
        if (x.eq.y) then
            d(x,y) = al
        elseif ((x+y).eq.3) then !picks out 1,2 or 2,1
            d(x,y) = be
        else
            d(x,y) = -be
        endif
    enddo
enddo


do x=2,5
    e(x,1,1) = mu(x)
    e(x,1,2) = nu(x)
    e(x,1,3) = de(x)
    e(x,2,1) = nu(x)
    e(x,2,2) = mu(x)
    e(x,2,3) = de(x)
    if (x.eq.2) then
        e(x,3,1) = -de(x) 
        e(x,3,2) = -de(x)
    elseif (x.eq.3) then
        e(x,3,1) = de(x)
        e(x,3,2) = de(x)
    elseif (x.eq.4) then
        e(x,3,1) = 0d0
        e(x,3,2) = 0d0
    elseif (x.eq.5) then
        e(x,3,1) = de(x)
        e(x,3,2) = de(x)
    endif
    e(x,3,3) = la(x)
enddo

do x=2,5
    f(x,1,1) = mu(x)
    f(x,1,2) = de(x)
    f(x,1,3) = nu(x)
    f(x,2,2) = la(x)
    if (x.eq.2) then
        f(x,2,1) = -de(x) 
        f(x,2,3) = -de(x)
    elseif (x.eq.3) then
        f(x,2,1) = de(x)
        f(x,2,3) = de(x)
    elseif (x.eq.4) then
        f(x,2,1) = 0d0
        f(x,2,3) = 0d0
    elseif (x.eq.5) then
        f(x,2,1) = de(x)
        f(x,2,3) = de(x)
    endif
    f(x,3,1) = nu(x)
    f(x,3,2) = de(x)
    f(x,3,3) = mu(x)
enddo

do x=2,5
    g(x,1,1) = la(x)
    if (x.eq.2) then
        g(x,1,2) = -de(x) 
        g(x,1,3) = -de(x)
    elseif (x.eq.3) then
        g(x,1,2) = de(x)
        g(x,1,3) = de(x)
    elseif (x.eq.4) then
        g(x,1,2) = 0d0
        g(x,1,3) = 0d0
    elseif (x.eq.5) then
        g(x,1,2) = de(x)
        g(x,1,3) = de(x)
    endif
    g(x,2,1) = de(x)
    g(x,2,2) = mu(x)
    g(x,2,3) = nu(x)
    g(x,3,1) = de(x)
    g(x,3,2) = nu(x)
    g(x,3,3) = mu(x)
enddo

do x=2,5
    h(x,1,1) = mu(x)
    h(x,1,2) = -nu(x)
    h(x,1,3) = de(x)
    h(x,2,1) = -nu(x)
    h(x,2,2) = mu(x)
    h(x,2,3) = -de(x)
    if (x.eq.2) then
        h(x,3,1) = -de(x) 
        h(x,3,2) = de(x)
    elseif (x.eq.3) then
        h(x,3,1) = de(x)
        h(x,3,2) = -de(x)
    elseif (x.eq.4) then
        h(x,3,1) = 0d0
        h(x,3,2) = 0d0
    elseif (x.eq.5) then
        h(x,3,1) = de(x)
        h(x,3,2) = -de(x)
    endif
    h(x,3,3) = la(x)
enddo

do x=2,5
    i(x,1,1) = mu(x)
    i(x,1,2) = de(x)
    i(x,1,3) = -nu(x)
    i(x,2,2) = la(x)
    if (x.eq.2) then
        i(x,2,1) = -de(x) 
        i(x,2,3) = de(x)
    elseif (x.eq.3) then
        i(x,2,1) = de(x)
        i(x,2,3) = -de(x)
    elseif (x.eq.4) then
        i(x,2,1) = 0d0
        i(x,2,3) = 0d0
    elseif (x.eq.5) then
        i(x,2,1) = de(x)
        i(x,2,3) = -de(x)
    endif
    i(x,3,1) = -nu(x)
    i(x,3,2) = -de(x)
    i(x,3,3) = mu(x)
enddo

do x=2,5
    j(x,1,1) = la(x)
    if (x.eq.2) then
        j(x,1,2) = -de(x) 
        j(x,1,3) = de(x)
    elseif (x.eq.3) then
        j(x,1,2) = de(x)
        j(x,1,3) = -de(x)
    elseif (x.eq.4) then
        j(x,1,2) = 0d0
        j(x,1,3) = 0d0
    elseif (x.eq.5) then
        j(x,1,2) = de(x)
        j(x,1,3) = -de(x)
    endif
    j(x,2,1) = de(x)
    j(x,2,2) = mu(x)
    j(x,2,3) = -nu(x)
    j(x,3,1) = -de(x)
    j(x,3,2) = -nu(x)
    j(x,3,3) = mu(x)
enddo

do x=2,5
    k(x,1,1) = mu(x)
    k(x,1,2) = -nu(x)
    k(x,1,3) = -de(x)
    k(x,2,1) = -nu(x)
    k(x,2,2) = mu(x)
    k(x,2,3) = de(x)
    if (x.eq.2) then
        k(x,3,1) = de(x) 
        k(x,3,2) = -de(x)
    elseif (x.eq.3) then
        k(x,3,1) = -de(x)
        k(x,3,2) = de(x)
    elseif (x.eq.4) then
        k(x,3,1) = 0d0
        k(x,3,2) = 0d0
    elseif (x.eq.5) then
        k(x,3,1) = -de(x)
        k(x,3,2) = de(x)
    endif
    k(x,3,3) = la(x)
enddo

do x=2,5
    l(x,1,1) = mu(x)
    l(x,1,2) = -de(x)
    l(x,1,3) = -nu(x)
    l(x,2,2) = la(x)
    if (x.eq.2) then
        l(x,2,1) = de(x) 
        l(x,2,3) = -de(x)
    elseif (x.eq.3) then
        l(x,2,1) = -de(x)
        l(x,2,3) = de(x)
    elseif (x.eq.4) then
        l(x,2,1) = 0d0
        l(x,2,3) = 0d0
    elseif (x.eq.5) then
        l(x,2,1) = -de(x)
        l(x,2,3) = de(x)
    endif
    l(x,3,1) = -nu(x)
    l(x,3,2) = de(x)
    l(x,3,3) = mu(x)
enddo

do x=2,5
    m(x,1,1) = la(x)
    if (x.eq.2) then
        m(x,1,2) = de(x) 
        m(x,1,3) = -de(x)
    elseif (x.eq.3) then
        m(x,1,2) = -de(x)
        m(x,1,3) = de(x)
    elseif (x.eq.4) then
        m(x,1,2) = 0d0
        m(x,1,3) = 0d0
    elseif (x.eq.5) then
        m(x,1,2) = -de(x)
        m(x,1,3) = de(x)
    endif
    m(x,2,1) = -de(x)
    m(x,2,2) = mu(x)
    m(x,2,3) = -nu(x)
    m(x,3,1) = de(x)
    m(x,3,2) = -nu(x)
    m(x,3,3) = mu(x)
enddo

do x=2,5
    n(x,1,1) = mu(x)
    n(x,1,2) = nu(x)
    n(x,1,3) = -de(x)
    n(x,2,1) = nu(x)
    n(x,2,2) = mu(x)
    n(x,2,3) = -de(x)
    if (x.eq.2) then
        n(x,3,1) = de(x) 
        n(x,3,2) = de(x)
    elseif (x.eq.3) then
        n(x,3,1) = -de(x)
        n(x,3,2) = -de(x)
    elseif (x.eq.4) then
        n(x,3,1) = 0d0
        n(x,3,2) = 0d0
    elseif (x.eq.5) then
        n(x,3,1) = -de(x)
        n(x,3,2) = -de(x)
    endif
    n(x,3,3) = la(x)
enddo

do x=2,5
    o(x,1,1) = mu(x)
    o(x,1,2) = -de(x)
    o(x,1,3) = nu(x)
    o(x,2,2) = la(x)
    if (x.eq.2) then
        o(x,2,1) = de(x) 
        o(x,2,3) = de(x)
    elseif (x.eq.3) then
        o(x,2,1) = -de(x)
        o(x,2,3) = -de(x)
    elseif (x.eq.4) then
        o(x,2,1) = 0d0
        o(x,2,3) = 0d0
    elseif (x.eq.5) then
        o(x,2,1) = -de(x)
        o(x,2,3) = -de(x)
    endif
    o(x,3,1) = nu(x)
    o(x,3,2) = -de(x)
    o(x,3,3) = mu(x)
enddo

do x=2,5
    p(x,1,1) = la(x)
    if (x.eq.2) then
        p(x,1,2) = de(x) 
        p(x,1,3) = de(x)
    elseif (x.eq.3) then
        p(x,1,2) = -de(x)
        p(x,1,3) = -de(x)
    elseif (x.eq.4) then
        p(x,1,2) = 0d0
        p(x,1,3) = 0d0
    elseif (x.eq.5) then
        p(x,1,2) = -de(x)
        p(x,1,3) = -de(x)
    endif
    p(x,2,1) = -de(x)
    p(x,2,2) = mu(x)
    p(x,2,3) = nu(x)
    p(x,3,1) = -de(x)
    p(x,3,2) = nu(x)
    p(x,3,3) = mu(x)
enddo

mat_dof(1,1,:,:) = a
mat_dof(1,2,:,:) = b
mat_dof(1,3,:,:) = c
mat_dof(1,4,:,:) = d
mat_dof(1,5,:,:) = e(2,:,:)
mat_dof(1,6,:,:) = f(2,:,:)
mat_dof(1,7,:,:) = g(2,:,:)
mat_dof(1,8,:,:) = h(2,:,:)
mat_dof(1,9,:,:) = i(2,:,:)
mat_dof(1,10,:,:) = j(2,:,:)
mat_dof(1,11,:,:) = k(2,:,:)
mat_dof(1,12,:,:) = l(2,:,:)
mat_dof(1,13,:,:) = m(2,:,:)
mat_dof(1,14,:,:) = n(2,:,:)
mat_dof(1,15,:,:) = o(2,:,:)
mat_dof(1,16,:,:) = p(2,:,:)
mat_dof(1,17,:,:) = e(3,:,:)
mat_dof(1,18,:,:) = f(3,:,:)
mat_dof(1,19,:,:) = g(3,:,:)
mat_dof(1,20,:,:) = h(3,:,:)
mat_dof(1,21,:,:) = i(3,:,:)
mat_dof(1,22,:,:) = j(3,:,:)
mat_dof(1,23,:,:) = k(3,:,:)
mat_dof(1,24,:,:) = l(3,:,:)
mat_dof(1,25,:,:) = m(3,:,:)
mat_dof(1,26,:,:) = n(3,:,:)
mat_dof(1,27,:,:) = o(3,:,:)
mat_dof(1,28,:,:) = p(3,:,:)
mat_dof(1,29,:,:) = e(4,:,:)
mat_dof(1,30,:,:) = f(4,:,:)
mat_dof(1,31,:,:) = g(4,:,:)
mat_dof(1,32,:,:) = h(4,:,:)
mat_dof(1,33,:,:) = i(4,:,:)
mat_dof(1,34,:,:) = j(4,:,:)
mat_dof(1,35,:,:) = e(5,:,:)
mat_dof(1,36,:,:) = f(5,:,:)
mat_dof(1,37,:,:) = g(5,:,:)
mat_dof(1,38,:,:) = h(5,:,:)
mat_dof(1,39,:,:) = i(5,:,:)
mat_dof(1,40,:,:) = j(5,:,:)
mat_dof(1,41,:,:) = k(5,:,:)
mat_dof(1,42,:,:) = l(5,:,:)
mat_dof(1,43,:,:) = m(5,:,:)
mat_dof(1,44,:,:) = n(5,:,:)
mat_dof(1,45,:,:) = o(5,:,:)
mat_dof(1,46,:,:) = p(5,:,:)

do x=1,46
    mat_dof(2,x,:,:) = transpose(mat_dof(1,x,:,:))
enddo

end subroutine ifc_si
!======================================
subroutine ifc_gaas(ifc_ga)

implicit none

integer :: basis(3),lvec(3),cnt
integer :: i,j,k,l,xi,xj,ai,aj,m1,m2,m3,mm1,mm2,mm3,a_num
real(8) :: ifc
real(8) :: ifc_r(6,6,6,2,2,3,3)
real(8), intent(out) :: ifc_ga(2,46,3,3)

ifc_r = 0d0
ifc_ga = 0d0

open(9998,file='GaAs.ifc',status='old',action='read')
do i=1,7
    read(9998,*)
enddo

basis(1) = 1
basis(2) = 1
basis(3) = 1

do i=1,7812
    if (mod(i,6**3+1).eq.1) then
        read(9998,*) xi,xj,ai,aj
    else
        read(9998,*) m1,m2,m3,ifc
        ifc_r(m1,m2,m3,ai,aj,xi,xj) = ifc

        if (m1.gt.4) then
            mm1 = m1-6
        else   
            mm1 = m1
        endif
        
        if (m2.gt.4) then
            mm2 = m2-6
        else   
            mm2 = m2
        endif
        
        if (m3.gt.4) then
            mm3 = m3-6
        else   
            mm3 = m3
        endif
    
        lvec(1) = 2*(mm1-1)+0*(mm2-1)+2*(mm3-1) - (aj-ai)*basis(1)
        lvec(2) = 0*(mm1-1)+2*(mm2-1)+2*(mm3-1) - (aj-ai)*basis(2)
        lvec(3) = 2*(mm1-1)+2*(mm2-1)+0*(mm3-1) - (aj-ai)*basis(3)

        call vec_neighbor(lvec,aj,a_num)
        

        if (a_num.ne.0) then
            if ((xi.eq.1).or.(xj.eq.1)) then
                if (xi.eq.xj) then
                    ifc_ga(aj,a_num,xi,xj) = ifc
                else
                    ifc_ga(aj,a_num,xi,xj) = -ifc
                endif
            else
                ifc_ga(aj,a_num,xi,xj) = ifc
            endif            
        endif
    endif
enddo
close(9998)

!open(9997,file='Gifc.ifc',status='unknown')
!do i=1,2
!    do j=1,46
!        do k=1,3
!            do l=1,3
!                write(9997,*)i,j,k,l,ifc_ga(i,j,k,l)
!            enddo
!        enddo
!    enddo
!enddo
!close(9997)


end subroutine ifc_gaas
!======================================
subroutine vec_neighbor(n_imp,atom_type,atom_number)

implicit none

integer, intent(in) :: n_imp(3),atom_type
integer :: n_vec(2,46,3),i
integer, intent(out) :: atom_number

call read_neighbors(n_vec)

atom_number = 0

do i=1,46
    if (n_imp(1).eq.n_vec(atom_type,i,1)) then
        if (n_imp(2).eq.n_vec(atom_type,i,2)) then 
            if (n_imp(3).eq.n_vec(atom_type,i,3)) then
                atom_number = i
            endif
        endif
    endif
enddo

end subroutine vec_neighbor
!======================================
subroutine super_dispersion(layer,bloks,num_atom)

use atomic
implicit none

integer, intent(in) :: bloks,num_atom
type(layers), intent(in) :: layer(bloks)
integer :: num_1,num_2,n_vec(2,46,3)
integer :: i,j,k,l,m,n,cnt
real(8) :: dyn(3*bloks*num_atom,3*bloks*num_atom),ifc(2,46,3,3),sum_f(3*bloks*num_atom)
real(8) :: pi,kpt(100,3),rk1,rk2,rk3,m_i,m_j
complex(8) :: dyn_k(100,3*bloks*num_atom,3*bloks*num_atom),ph

pi=3.14159265358979323846d0

!call ifc_si(ifc)

call ifc_gaas(ifc)

call read_neighbors(n_vec)

dyn = 0d0
dyn_k = 0d0

do i=1,bloks
    do j=1,num_atom
        do k=1,46
            num_1 = (i-1)*num_atom*3+(j-1)*3 
            num_2 = (layer(i)%a(j)%n(k,1)-1)*num_atom*3+(layer(i)%a(j)%n(k,2)-1)*3
            do m=1,3
                do n=1,3
                    dyn(num_1+m,num_2+n) = ifc(layer(i)%a(j)%t,k,m,n)+dyn(num_1+m,num_2+n)
                enddo
            enddo
        enddo
    enddo
enddo


do i=1,(3*bloks*num_atom)
    sum_f(i) = 0d0
    do j=1,(3*bloks*num_atom)
        sum_f(i) = sum_f(i) + dyn(i,j)
    enddo
    print*,sum_f(i)
    dyn(i,i) = dyn(i,i)-sum_f(i)
enddo

do i=1,100
    kpt(i,1) = (0d0)*(pi)*(real(i-1))/(real(100-1))/real(bloks) !need to change to xx
    kpt(i,2) = (0d0)*(pi)*(real(i-1))/(real(100-1))/real(bloks) !need to change to yy
    kpt(i,3) = (1d0)*(pi)*(real(i-1))/(real(100-1))/real(bloks)
enddo

do l=1,100
    do i=1,bloks
        do j=1,num_atom
            do k=1,46
                num_1 = (i-1)*num_atom*3+(j-1)*3 
                num_2 = (layer(i)%a(j)%n(k,1)-1)*num_atom*3+(layer(i)%a(j)%n(k,2)-1)*3
                rk1 = kpt(l,1)*n_vec(layer(i)%a(j)%t,k,1)/4d0  
                rk2 = kpt(l,2)*n_vec(layer(i)%a(j)%t,k,2)/4d0 
                rk3 = kpt(l,3)*n_vec(layer(i)%a(j)%t,k,3)/4d0 
                ph = dcmplx(0,1)*(rk1+rk2+rk3)
                m_i = layer(i)%a(j)%m
                m_j = layer(layer(i)%a(j)%n(k,1))%a(layer(i)%a(j)%n(k,2))%m
                do m=1,3
                    do n=1,3
dyn_k(l,num_1+m,num_2+n) = ifc(layer(i)%a(j)%t,k,m,n)*zexp(ph)/sqrt(m_i*m_j)+dyn_k(l,num_1+m,num_2+n) 
                    enddo
                enddo
            enddo
        enddo
    enddo
    cnt = 0
    do i=1,bloks
        do j=1,num_atom
            do m=1,3
                cnt = cnt+1
                m_i = layer(i)%a(j)%m
                dyn_k(l,cnt,cnt) = dyn_k(l,cnt,cnt)-sum_f(cnt)/m_i
            enddo
        enddo
    enddo
enddo

open(2004,file='superifc.dat',status='unknown')
do k=1,100
    do i=1,(3*bloks*num_atom)
        do j=1,(3*bloks*num_atom)
            write(2004,*)k,i,j,real(dyn_k(k,i,j)),imag(dyn_k(k,i,j))
        enddo
    enddo
enddo
close(2004)

end subroutine super_dispersion
!======================================
subroutine bulk_dispersion()

implicit none

integer :: i,j,k,l,m,n,nk,nn
integer :: n_vec(2,46,3)
real(8) :: dyn(6,6)
real(8) :: ifc(2,46,3,3),sum_f(6),pi
real(8) :: kpt(100,3),rk1,rk2,rk3
complex(8) :: dyn_k(100,6,6),ph
real(8) :: b1(3),b2(3),b3(3)
real(8) :: mass(2)

pi=3.14159265358979323846d0

!call ifc_si(ifc)

call ifc_gaas(ifc)

call read_neighbors(n_vec)

!m1 = 69.723d0  !(gallium)
!m2 = 26.9827d0 !(aluminum)
!m3 = 74.9216d0  !(arsenic)

mass(1) = 69.723d0 !(gallium)
!mass(1) = 26.9827d0 !(aluminum)
mass(2) = 74.9216d0 !(arsenic)

dyn = 0d0

do i=1,2
    do j=1,46
        l = modulo((i-1)*3+n_vec(i,j,1)+n_vec(i,j,2)+n_vec(i,j,3),2)+1
        do m=1,3
            do n=1,3
                dyn(3*(i-1)+m,3*(l-1)+n) = ifc(i,j,m,n)+dyn(3*(i-1)+m,3*(l-1)+n)
                !print*,ifc(i,j,m,n)
            enddo
        enddo
    enddo
enddo


!diagonal components

do i=1,6
    sum_f(i) = 0d0
    do j=1,6
        sum_f(i) = sum_f(i) + dyn(i,j) 
    enddo
    dyn(i,i) = dyn(i,i)-sum_f(i)    
enddo     

b1(1) = (2d0*pi)*(-1d0)
b1(2) = (2d0*pi)*(1d0)
b1(3) = (2d0*pi)*(1d0) 

b2(1) = (2d0*pi)*(1d0)
b2(2) = (2d0*pi)*(-1d0)
b2(3) = (2d0*pi)*(1d0)

b3(1) = (2d0*pi)*(1d0)
b3(2) = (2d0*pi)*(1d0)
b3(3) = (2d0*pi)*(-1d0)


nk = 100
!nn = 40
!do i=1,nn
!    do j=1,nn
!        do k=1,nn
!            nk = nk+1
!kpt(nk,1) = real(i-1)*b1(1)/real(nn+1)+real(j-1)*b2(1)/real(nn+1)+real(k-1)*b3(1)/real(nn+1)
!kpt(nk,2) = real(i-1)*b1(2)/real(nn+1)+real(j-1)*b2(2)/real(nn+1)+real(k-1)*b3(2)/real(nn+1)
!kpt(nk,3) = real(i-1)*b1(3)/real(nn+1)+real(j-1)*b2(3)/real(nn+1)+real(k-1)*b3(3)/real(nn+1)
!        enddo
!    enddo
!enddo

do i=1,nk
    kpt(i,1) = (0d0/4d0)*(2d0*pi)*(real(i-1))/(real(100-1))
    kpt(i,2) = (0d0/4d0)*(2d0*pi)*(real(i-1))/(real(100-1))
    kpt(i,3) = (4d0/4d0)*(2d0*pi)*(real(i-1))/(real(100-1))
enddo

do k=1,nk
    do i=1,2
        do j=1,46
            l = modulo((i-1)*3+n_vec(i,j,1)+n_vec(i,j,2)+n_vec(i,j,3),2)+1


            rk1 = kpt(k,1)*n_vec(i,j,1)/4d0  
            rk2 = kpt(k,2)*n_vec(i,j,2)/4d0 
            rk3 = kpt(k,3)*n_vec(i,j,3)/4d0 
            ph = dcmplx(0,1)*(rk1+rk2+rk3)
            do m=1,3
                do n=1,3
dyn_k(k,3*(i-1)+m,3*(l-1)+n) = ifc(i,j,m,n)*zexp(ph)/sqrt(mass(i)*mass(l))+dyn_k(k,3*(i-1)+m,3*(l-1)+n)
                enddo
            enddo
        enddo
        do j=1,3
            dyn_k(k,3*(i-1)+j,3*(i-1)+j) = dyn_k(k,3*(i-1)+j,3*(i-1)+j)-sum_f(3*(i-1)+j)/mass(i)
        enddo
    enddo
    !do i=1,6
    !    dyn_k(k,i,i) = dyn_k(k,i,i)-sum_f(i)
    !enddo
enddo


open(2004,file='bulkifc.dat',status='unknown')
do k=1,nk
    do i=1,6
        do j=1,6
            write(2004,*)k,i,j,real(dyn_k(k,i,j)),imag(dyn_k(k,i,j))
        enddo
    enddo
enddo
close(2004)


end subroutine bulk_dispersion
!======================================
!======================================
subroutine setup_dyn(layer,bloks,num_atom,ifc_d,ifc_a,dyn_d,dyn_a)

use atomic
implicit none

integer :: i,j,k,m,n
integer, intent(in) :: bloks,num_atom
real(8) :: m_i,m_j
real(8), intent(in) :: ifc_d(3*num_atom,3*num_atom),ifc_a(3*num_atom,3*num_atom)
complex(8), intent(out) :: dyn_d(bloks,3*num_atom,3*num_atom),dyn_a(bloks-1,3*num_atom,3*num_atom)
type(layers), intent(in) :: layer(bloks)

do i=1,bloks
    do j=1,num_atom
        do k=1,num_atom
            m_i = layer(i)%a(j)%m
            m_j = layer(i)%a(k)%m
            do m=1,3
                do n=1,3
                    dyn_d(i,3*(j-1)+m,3*(k-1)+n) = ifc_d(3*(j-1)+m,3*(k-1)+n)/sqrt(m_i*m_j)
                enddo
            enddo
        enddo
    enddo
enddo

do i=1,bloks-1
    do j=1,num_atom
        do k=1,num_atom
            m_i = layer(i)%a(j)%m
            m_j = layer(i+1)%a(k)%m
            do m=1,3
                do n=1,3
                    dyn_a(i,3*(j-1)+m,3*(k-1)+n) = ifc_a(3*(j-1)+m,3*(k-1)+n)/sqrt(m_i*m_j)
                enddo
            enddo
        enddo  
    enddo
enddo

end subroutine setup_dyn
!======================================
end module superlattice





















