module my_functions
implicit none
contains

!==========================================================
function matvec_z(a,x) result(y)

integer:: m,n,j
complex(8) :: a(:,:),x(:)
complex(8) :: alpha, beta
complex(8), allocatable :: y(:)

alpha = 1d0
beta = 0d0

m = size(a(:,1))
n = size(a(1,:))

allocate(y(m))

call zgemv('n',m,n,alpha,a,m,x,1,beta,y,1)

end function matvec_z
!==========================================================
function matmul_d(a,b) result(c)


integer :: n,m,j

real(8) :: a(:,:),b(:,:)
real(8) :: alpha, beta
real(8), allocatable :: c(:,:)

!complex(8) :: a(n,j),b(j,m)
!complex(8) :: c(n,m), alpha, beta

n = size(a(:,1))
m = size(b(1,:))
j = size(a(1,:))

allocate(c(n,m))

alpha = 1d0
beta = 0d0

call dgemm('N','N',n,m,j,alpha,a,n,b,j,beta,c,n)

end function matmul_d
!==========================================================
function matmul_g(a,b) result(c)


integer :: n,m,j

complex(8) :: a(:,:),b(:,:)
complex(8) :: alpha, beta
complex(8), allocatable :: c(:,:)

!complex(8) :: a(n,j),b(j,m)
!complex(8) :: c(n,m), alpha, beta

n = size(a(:,1))
m = size(b(1,:))
j = size(a(1,:))

allocate(c(n,m))

alpha = cmplx(1d0,0d0)
beta = 0d0

call zgemm('N','N',n,m,j,alpha,a,n,b,j,beta,c,n)

end function matmul_g
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
function green_calc_sum(E,k_pt,eig_f,eig_v,R_i,R_j,atom_i,atom_j,alpha,beta) result(green)

implicit none

integer :: cnt,i,j,k
complex(8) :: green,eig_v(:,:,:,:),f_num
real(8) :: E,k_pt(:,:,:),eig_f(:,:,:),R_i(:),R_j(:),phase,smear
integer :: atom_i,atom_j,alpha,beta

green = 0d0
cnt = 0

do i=1,size(k_pt,1)
    do j=1,size(k_pt,2)
        do k=1,size(eig_f,3)
            phase = dot_product(k_pt(i,j,:),R_i-R_j)
            f_num = zexp(dcmplx(0,1)*phase)
            f_num = f_num*eig_v(i,j,k,3*(atom_i-1)+alpha)*dconjg(eig_v(i,j,k,3*(atom_j-1)+beta))

            if ((i.eq.1).and.(j.eq.1)) then
                cnt = cnt+1
                smear = (abs(eig_f(i,j+1,k)-eig_f(i,j,k))+abs(eig_f(i+1,j,k)-eig_f(i,j,k)))/4d0
                green = green + f_num/(E-eig_f(i,j,k)+dcmplx(0,1)*smear)

            elseif (((i.eq.1).and.((j.gt.1)).and.(j.lt.size(k_pt,2)))) then
                cnt = cnt+2
                smear = (abs(eig_f(i,j+1,k)-eig_f(i,j,k))+abs(eig_f(i+1,j,k)-eig_f(i,j,k))&
                        +abs(eig_f(i,j,k)-eig_f(i,j-1,k)))/6d0
                green = green + 2d0*f_num/(E-eig_f(i,j,k)+dcmplx(0,1)*smear)

            elseif ((i.eq.1).and.(j.eq.size(k_pt,2))) then
                cnt = cnt+1
                smear = (abs(eig_f(i+1,j,k)-eig_f(i,j,k))+abs(eig_f(i,j,k)-eig_f(i,j-1,k)))/4d0
                green = green + f_num/(E-eig_f(i,j,k)+dcmplx(0,1)*smear)

            elseif (((i.gt.1).and.(i.lt.size(k_pt,1))).and.(j.eq.1)) then
                cnt = cnt+2
                smear = (abs(eig_f(i+1,j,k)-eig_f(i,j,k))+abs(eig_f(i,j,k)-eig_f(i-1,j,k))&
                        +abs(eig_f(i,j+1,k)-eig_f(i,j,k)))/6d0
                green = green + 2d0*f_num/(E-eig_f(i,j,k)+dcmplx(0,1)*smear)

            elseif ((j.eq.1).and.(i.eq.size(k_pt,1))) then
                cnt = cnt+1
                smear = (abs(eig_f(i,j,k)-eig_f(i-1,j,k))+abs(eig_f(i,j+1,k)-eig_f(i,j,k)))/4d0
                green = green + f_num/(E-eig_f(i,j,k)+dcmplx(0,1)*smear)

            elseif ((i.eq.size(k_pt,1)).and.((j.gt.1).and.(j.lt.size(k_pt,2)))) then
                cnt = cnt+2
                smear = (abs(eig_f(i,j+1,k)-eig_f(i,j,k))+abs(eig_f(i,j,k)-eig_f(i-1,j,k)) &
                         +abs(eig_f(i,j,k)-eig_f(i,j-1,k)))/6d0
                green = green + 2d0*f_num/(E-eig_f(i,j,k)+dcmplx(0,1)*smear)

            elseif ((j.eq.size(k_pt,2)).and.((i.gt.1).and.(i.lt.size(k_pt,1)))) then
                cnt = cnt+2
                smear = (abs(eig_f(i+1,j,k)-eig_f(i,j,k))+abs(eig_f(i,j,k)-eig_f(i-1,j,k)) &
                        +abs(eig_f(i,j,k)-eig_f(i,j-1,k)))/6d0
                green = green + 2d0*f_num/(E-eig_f(i,j,k)+dcmplx(0,1)*smear)

            elseif ((i.eq.size(k_pt,1)).and.(j.eq.size(k_pt,2))) then
                cnt = cnt+1
                smear = (abs(eig_f(i,j,k)-eig_f(i-1,j,k))+abs(eig_f(i,j,k)-eig_f(i,j-1,k)))/4d0
                green = green + f_num/(E-eig_f(i,j,k)+dcmplx(0,1)*smear)

            elseif (((i.gt.1).and.(i.lt.size(k_pt,1))).and.((j.gt.1).and.(j.lt.size(k_pt,2)))) then
                cnt = cnt+4
                smear = (abs(eig_f(i+1,j,k)-eig_f(i,j,k))+abs(eig_f(i,j,k)-eig_f(i-1,j,k))&
                        +abs(eig_f(i,j+1,k)-eig_f(i,j,k))+abs(eig_f(i,j,k)-eig_f(i,j-1,k)))/8d0
                green = green + 4d0*f_num/(E-eig_f(i,j,k)+dcmplx(0,1)*smear)

            endif
        enddo
    enddo
enddo


cnt = cnt/size(eig_f,3)
green = green/real(cnt)

end function green_calc_sum
!==========================================================
subroutine ifc_phosphorene(ifc_r)

implicit none

integer :: i,xi,xj,ai,aj,m1,m2,m3
real(8) :: ifc
real(8), intent(out) :: ifc_r(8,8,1,4,4,3,3)

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

end subroutine ifc_phosphorene
!==========================================================
subroutine generate_bond(phi_odd,phi_even,D_chain,chain_length,mass)

implicit none

integer :: i,j,m,n
real(8), intent(in) :: mass
integer, intent(in) :: chain_length
real(8), intent(out) :: phi_odd(3,3),phi_even(3,3),D_chain(3*chain_length,3*chain_length)
real(8) :: theta,R_z(3,3),R_y(3,3),R_1(3,3),R_2(3,3),D_bond(3,3),D1(3,3),D2(3,3)
real(8) :: ifc_r(8,8,1,4,4,3,3)

!grab the bond

D_chain = 0d0
phi_odd = 0d0
phi_even = 0d0


call ifc_phosphorene(ifc_r)

do i=1,3
    do j=1,3
        D_bond(i,j) = ifc_r(1,1,1,2,3,i,j)
    enddo
enddo


!odd bond
theta = -0.7343d0
R_z = rotate_z(theta)
theta = -2.25147d0
theta = -1.57079632679d0
R_y = rotate_y(theta)
R_1 = matmul_d(R_y,R_z)
theta = 0.7343d0
R_z = rotate_z(theta)
R_1 = matmul_d(R_z,R_1)

!even bond
theta = -0.7343d0
R_z = rotate_z(theta)
theta = -0.890118d0
theta = -1.57079632679d0
R_y = rotate_y(theta)
R_2 = matmul_d(R_y,R_z)
theta = 0.7343d0
R_z = rotate_z(theta)
R_2 = matmul_d(R_z,R_2)

!bond setup
D1(1:3,1:3) = 0d0
D2(1:3,1:3) = 0d0

do i=1,3
    do j=1,3
        do m=1,3
            do n=1,3
                !tensor contraction
                D1(i,j) = D1(i,j) + D_bond(m,n)*R_1(m,i)*R_1(n,j) !odd bond
                D2(i,j) = D2(i,j) + D_bond(m,n)*R_2(m,i)*R_2(n,j) !even bond
            enddo
        enddo
    enddo
enddo

phi_odd = D1
phi_even = D2

do i=1,(chain_length-1)
    do m=1,3
        do n=1,3
            if (mod(i,2).eq.1) then
                D_chain(3*(i-1)+m,3*i+n) = D2(m,n)
            elseif (mod(i,2).eq.0) then !even
                D_chain(3*(i-1)+m,3*i+n) = D1(m,n)
            endif
        enddo
    enddo
enddo

D_chain = D_chain+transpose(D_chain)

do i=1,3*chain_length
    D_chain(i,i) = -sum(D_chain(i,:))
enddo

if (chain_length.eq.1) then
    D_chain = 0d0
endif

do i=1,3*chain_length
    do j=1,3*chain_length
        D_chain(i,j) = D_chain(i,j)/mass
    enddo
enddo

D_chain(3*(chain_length-1)-2:3*(chain_length-1),3*(chain_length)-2:3*(chain_length)) = &
D_chain(3*(chain_length-1)-2:3*(chain_length-1),3*(chain_length)-2:3*(chain_length))/dsqrt(1000d0)
D_chain(3*(chain_length)-2:3*(chain_length),3*(chain_length-1)-2:3*(chain_length-1)) = &
D_chain(3*(chain_length)-2:3*(chain_length),3*(chain_length-1)-2:3*(chain_length-1))/dsqrt(1000d0)
D_chain(3*(chain_length)-2:3*(chain_length),3*(chain_length)-2:3*(chain_length)) = &
D_chain(3*(chain_length)-2:3*(chain_length),3*(chain_length)-2:3*(chain_length))/1000d0

end subroutine generate_bond
!==========================================================
function rotate_z(theta) result(R_z)

implicit none

real(8) :: theta,R_z(3,3)

R_z(1,1) = cos(theta)
R_z(1,2) = -sin(theta)
R_z(1,3) = 0d0

R_z(2,1) = sin(theta)
R_z(2,2) = cos(theta)
R_z(2,3) = 0d0

R_z(3,1) = 0d0
R_z(3,2) = 0d0
R_z(3,3) = 1d0

end function rotate_z
!==========================================================
function rotate_y(theta) result(R_y)

implicit none

real(8) :: theta,R_y(3,3)

R_y(1,1) = cos(theta)
R_y(1,2) = 0d0
R_y(1,3) = sin(theta)

R_y(2,1) = 0d0
R_y(2,2) = 1d0
R_y(2,3) = 0d0

R_y(3,1) = -sin(theta)
R_y(3,2) = 0d0
R_y(3,3) = cos(theta)

end function rotate_y
!==========================================================
subroutine read_spectrum(k_pt,eig_f,eig_v)

implicit none

integer :: i,j,k,m,n
complex(8) :: eig_v(60,60,12,12)
real(8) :: eig_f(60,60,12),k_pt(60,60,2)
real(8) :: eig,k_x,k_y,eig_r,eig_i



open(6900,file='eig_freq.dat',status='old',action='read')
do i=1,(60*60*12)
    read(6900,*) j,k,m,eig
    eig_f(j,k,m) = eig
enddo
close(6900)

open(6901,file='k_pt.dat',status='old',action='read')
do i=1,(60*60)
    read(6901,*) j,k,k_x,k_y
    k_pt(j,k,1) = k_x
    k_pt(j,k,2) = k_y
enddo
close(6901)

open(6902,file='eig_vec.dat',status='old',action='read')
do i=1,(60*60*12*12)
    read(6902,*) j,k,m,n,eig_r,eig_i
    eig_v(j,k,m,n) = eig_r + dcmplx(0,1)*eig_i
enddo
close(6902)

end subroutine read_spectrum
!==========================================================
subroutine t_matrix(t_mat,k_pt,eig_f,eig_v,k_x,k_y,branch,chain_length,phi_odd,D_chain)

implicit none

integer :: i,j,atom_i,atom_j
integer :: k_x,k_y,branch,chain_length
complex(8) :: dummy(3*(chain_length+1),3*(chain_length+1))
complex(8) :: res_chain(3*chain_length,3*chain_length),res_tot(3*(chain_length+1),3*(chain_length+1))
complex(8) :: eig_v(:,:,:,:),g_surf(3,3),g_chain(3*chain_length,3*chain_length),g_tot(3*(chain_length+1),3*(chain_length+1))
complex(8) :: V_surf(3,3),V_m(3,3),V_pert(3*(chain_length+1),3*(chain_length+1)),t_mat(3*(chain_length+1),3*(chain_length+1))
real(8) :: k_pt(:,:,:),eig_f(:,:,:),R_i(2),R_j(2),epsilon,mass,D_chain(3*chain_length,3*chain_length),phi_odd(3,3)

R_i = 0d0
R_j = 0d0
atom_i = 2
atom_j = 2
V_surf = 0d0
V_m = 0d0
V_pert = 0d0
t_mat = 0d0
res_chain = 0d0
res_tot = 0d0
g_surf = 0d0
g_chain = 0d0
g_tot = 0d0

mass = 30.973762d0

do i=1,3
    do j=1,3
        g_surf(i,j) = green_calc_sum(eig_f(k_x,k_y,branch),k_pt,eig_f,eig_v,R_i,R_j,atom_i,atom_j,i,j)
    enddo
enddo

epsilon = 1.0d-5

do i=1,3*chain_length
    do j=1,3*chain_length
        if (i.eq.j) then
            res_chain(i,j) = eig_f(k_x,k_y,branch) + dcmplx(0,1)*epsilon - D_chain(i,j)
        elseif (i.ne.j) then
            res_chain(i,j) = -D_chain(i,j)
        endif
    enddo
enddo

call z_inv(res_chain,g_chain,3*chain_length)

do i=1,3
    do j=1,3
        g_tot(i,j) = g_surf(i,j)
    enddo
enddo

do i=1,3*chain_length
    do j=1,3*chain_length
        g_tot(i+3,j+3) = g_chain(i,j)
    enddo
enddo

do i=1,3
    V_surf(i,i) = -sum(phi_odd(i,:))
    V_m(i,i) = -sum(phi_odd(:,i)) !transpose of the bond
enddo

do i=1,3
    do j=1,3
        V_pert(i,j+3) = phi_odd(i,j)/mass
        V_pert(i+3,j) = phi_odd(j,i)/mass
        V_pert(i,j) = V_surf(i,j)/mass
        V_pert(i+3,j+3) = V_m(i,j)/mass
    enddo
enddo

dummy = matmul_g(V_pert,g_tot)

do i=1,3*(chain_length+1)
    do j=1,3*(chain_length+1)
        if (i.eq.j) then
            res_tot(i,j) = 1d0 - dummy(i,j)
        elseif (i.ne.j) then
            res_tot(i,j) = -dummy(i,j)
        endif
    enddo
enddo

call z_inv(res_tot,dummy,3*(chain_length+1))

t_mat = matmul_g(dummy,V_pert)


end subroutine t_matrix
!==========================================================
function scat_rate(k_pt,eig_f,eig_v,k_x,k_y,branch,chain_length) result(scat)

implicit none

integer :: i
integer :: k_x,k_y,branch,chain_length
real(8) :: eig_f(:,:,:),scat,k_pt(:,:,:)
real(8) :: conv_amu,conv_bohr,conv_ryd,conv,omeg
complex(8) :: t_mat(3*(chain_length+1),3*(chain_length+1)),k_vec(3*(chain_length+1)),k_dum(3*(chain_length+1))
complex(8) :: eig_v(:,:,:,:),forward_amp
real(8) :: phi_odd(3,3),phi_even(3,3),D_chain(3*chain_length,3*chain_length)
real(8) :: mass


mass = 30.973762d0 !mass of chain (phosphorene for now)
conv_ryd = 2.1798741d-18
conv_bohr = 5.2917721092d-11
conv_amu = 1.66054d-27

conv = conv_ryd/(conv_bohr)**2/conv_amu

omeg = dsqrt(conv*eig_f(k_x,k_y,branch)) !rad/sec

do i=1,3*(chain_length+1)
    if (i.le.3) then
        k_vec(i) = eig_v(k_x,k_y,branch,i+3) !i+3 specifies 4:6 which is atom2 in the unit cell
    elseif (i.gt.3) then
        k_vec(i) = 0d0
    endif
enddo

call generate_bond(phi_odd,phi_even,D_chain,chain_length,mass)

call t_matrix(t_mat,k_pt,eig_f,eig_v,k_x,k_y,branch,chain_length,phi_odd,D_chain)

k_dum = matvec_z(t_mat,k_vec)

!do i=1,3
!    print*,t_mat(i,1),t_mat(i,2),t_mat(i,3)
!enddo
!
!do i=1,3*(chain_length+1)
!    print*,k_vec(i),k_dum(i)
!enddo
!
forward_amp = 0d0

do i=1,3*(chain_length+1)
    forward_amp = forward_amp + dconjg(k_vec(i))*k_dum(i)
enddo

scat = -conv/omeg*dimag(forward_amp) !scattering rate per (functionalized sites per unit cell)


end function scat_rate
!==========================================================
function lattice_vec(mm1,mm2) result(l_vec)

implicit none

integer :: mm1,mm1
real(8) :: a,b,a1(2),a2(2),l_vec(2)

a = 3.297d0
b = 1.4006161d0*a

a1(1) = a
a1(2) = 0d0

a2(1) = 0d0
a2(2) = b

if (mm1.gt.5) then
    mm1 = mm1-8
endif

if (mm2.gt.5) then
    mm2 = mm2-8
endif

l_vec = real((mm1-1))*a1 + real((mm2-1))*a2

end function lattice_vec
!==========================================================
function dyn_mat(k_pt,mass) result(dyn)

implicit none

integer :: m1,m2,ai,aj,xi,xj,i,j
complex(8) :: dyn(12,12)
real(8) :: ifc_r(8,8,1,4,4,3,3),vec(2),mass,k_pt(2)

dyn = 0d0

call ifc_phosphorene(ifc_r)

do m1=1,8
    do m2=1,8
        vec = lattice_vec(m1,m2)
        do ai=1,4
            do aj=1,4
                do xi=1,3
                    do xj=1,3
dyn(3*(ai-1)+xi,3*(aj-1)+xj) = dyn(3*(ai-1)+xi,3*(aj-1)+xj)&
+ 1d0/mass*ifc_r(m1,m2,1,ai,aj,xi,xj)*zexp(dcmplx(0,1)*dot_product(k_pt,vec))
                    enddo
                enddo
            enddo
        enddo
    enddo
enddo

!symmetrization step
do i=1,12
    do j=1,i
        if (i.eq.j) then
            dyn(i,j) = dreal(dyn(i,j))
        elseif (i.ne.j) then
            dyn(i,j) = dconjg(dyn(j,i))
        endif
    enddo
enddo

end function dyn_mat
!==========================================================
end module my_functions












