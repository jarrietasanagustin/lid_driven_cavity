subroutine stream_function(omega_in,psi_in,psi_out,n_x_i,n_z_i)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!Solver linea por linea para determinar la
!la funcion de corriente y el campo de velocidades
!en la primera etapa del RK-3
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
use definiciones
implicit none

real(8)::result_1,result_2
real(8),dimension((n_x_i-1)*(n_z_i-1))::p,residuo,kk,residuo_prima,z_cg,vector_psi
real(8),dimension(n_x_i+1,n_z_i+1),intent(in)::omega_in,psi_in
real(8),dimension(n_x_i+1,n_z_i+1),intent(out)::psi_out
integer(4)::i,j,k,k_t,iter,n_x_i,n_z_i
real(8)::a_gauss,b_gauss,c_gauss,d_gauss,f_gauss,e_gauss,inv_d
real(8)::alpha,beta
real(8)::max_res
real(8)::producto,producto_1,producto_residuo_prima
integer(4)::num_fila
integer(4)::num_threads,save_nt


!include 'mkl.fi'
!num_threads=4
!save_nt = mkl_set_num_threads_local(num_threads)
num_fila=(nx_i-1)*(nz_i-1)
!---------------------------------------
!Condiciones de contorno
!Valor de la funcion de corriente =0 en las
!paredes
!---------------------------------------
!write(*,*)omega_in

psi_out=psi_in
do j=1,n_z_i+1
i=1
!psi_out(i,j)=0.0d0
psi_out(i,j)=0.0d0
i=n_x_i+1
!psi_out(i,j)=0.0d0
psi_out(i,j)=0.0d0
end do
do i=1,n_x_i+1
j=1
!psi_out(i,j)=0.0d0
psi_out(i,j)=0.0d0
j=n_z_i+1
!psi_out(i,j)=0.0d0
psi_out(i,j)=0.0d0
end do

!%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!Coeficientes de la matriz
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b_gauss=1.0d0/delta_z**2

a_gauss=1.0d0/delta_z**2

d_gauss=-2.0d0/delta_z**2-2.0d0/delta_x**2

e_gauss=1.0d0/delta_x**2

f_gauss=1.0d0/delta_x**2

inv_d=1.0d0/d_gauss
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!Esquina inferior izquierda
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=1
j=1
c_gauss=-omega_in(i+1,j+1)
c_gauss=c_gauss-f_gauss*psi_out(i,j+1)-b_gauss*psi_out(i+1,j)
residuo((i-1)*(n_z_i-1)+j)=(c_gauss-e_gauss*psi_out(i+2,j+1)-a_gauss*psi_out(i+1,j+2)-d_gauss*psi_out(i+1,j+1))


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!Pared lateral izquierda
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
do  j=2,n_z_i-2
c_gauss=-omega_in(i+1,j+1)
c_gauss=c_gauss-f_gauss*psi_out(i,j+1)
residuo((i-1)*(n_z_i-1)+j)=(c_gauss-b_gauss*psi_out(i+1,j)-a_gauss*psi_out(i+1,j+2)-&
        e_gauss*psi_out(i+2,j+1)-d_gauss*psi_out(i+1,j+1))
end do

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!Esquina superior izquierda
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
j=n_z_i-1
c_gauss=-omega_in(i+1,j+1)
c_gauss=c_gauss-f_gauss*psi_out(i,j+1)-a_gauss*psi_out(i+1,j+2)
residuo((i-1)*(n_z_i-1)+j)=(c_gauss-b_gauss*psi_out(i+1,j)-e_gauss*psi_out(i+2,j+1)-d_gauss*psi_out(i+1,j+1))

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%Interior del dominio
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
do i=2,n_x_i-2

j=1
c_gauss=-omega_in(i+1,j+1)
c_gauss=c_gauss-b_gauss*psi_out(i+1,j)
residuo((i-1)*(n_z_i-1)+j)=(c_gauss-a_gauss*psi_out(i+1,j+2)-e_gauss*psi_out(i+2,j+1)-f_gauss*psi_out(i,j+1)-&
        d_gauss*psi_out(i+1,j+1))

    do j=2,n_z_i-2
    c_gauss=-omega_in(i+1,j+1)
    residuo((i-1)*(n_z_i-1)+j)=(c_gauss-b_gauss*psi_out(i+1,j)-a_gauss*psi_out(i+1,j+2)&
    -e_gauss*psi_out(i+2,j+1)-f_gauss*psi_out(i,j+1)-d_gauss*psi_out(i+1,j+1))
    end do
j=n_z_i-1
c_gauss=-omega_in(i+1,j+1)
c_gauss=c_gauss-a_gauss*psi_out(i+1,j+2)
residuo((i-1)*(n_z_i-1)+j)=(c_gauss-b_gauss*psi_out(i+1,j)-e_gauss*psi_out(i+2,j+1)&
-f_gauss*psi_out(i,j+1)-d_gauss*psi_out(i+1,j+1))
end do

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%x=L
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=n_x_i-1
j=1
c_gauss=-omega_in(i+1,j+1)
c_gauss=c_gauss-e_gauss*psi_out(i+2,j+1)-b_gauss*psi_out(i+1,j)
residuo((i-1)*(n_z_i-1)+j)=(c_gauss-a_gauss*psi_out(i+1,j+2)-f_gauss*psi_out(i,j+1)-d_gauss*psi_out(i+1,j+1))

do j=2,n_z_i-2
c_gauss=-omega_in(i+1,j+1)
c_gauss=c_gauss-e_gauss*psi_out(i+2,j+1)
residuo((i-1)*(n_z_i-1)+j)=(c_gauss-b_gauss*psi_out(i+1,j)-a_gauss*psi_out(i+1,j+2)-&
                f_gauss*psi_out(i,j+1)-d_gauss*psi_out(i+1,j+1))
end do

j=n_z_i-1
c_gauss=-omega_in(i+1,j+1)
c_gauss=c_gauss-e_gauss*psi_out(i+2,j+1)-a_gauss*psi_out(i+1,j+2)
residuo((i-1)*(n_z_i-1)+j)=(c_gauss-b_gauss*psi_out(i+1,j)-f_gauss*psi_out(i,j+1)-d_gauss*psi_out(i+1,j+1))


iter=1
max_res=maxval(dabs(residuo))
z_cg=inv_d*residuo
p(1:(n_x_i-1)*(n_z_i-1))=z_cg(1:(n_x_i-1)*(n_z_i-1))
producto=dot_product(residuo,z_cg)
!producto=ddot(num_fila,residuo,1,z_cg,1)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
do j=1,n_z_i-1
    do i=1,n_x_i-1
    vector_psi((i-1)*(n_z_i-1)+j)=psi_out(i+1,j+1)
    end do
end do
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


do while(max_res>eps)
!$OMP PARALLEL PRIVATE(i,j)  NUM_THREADS(4)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i=1
j=1
kk((i-1)*(n_z_i-1)+j)=-(-e_gauss*p((i)*(n_z_i-1)+j)-a_gauss*p((i-1)*(n_z_i-1)+j+1)-d_gauss*p((i-1)*(n_z_i-1)+j))

do  j=2,n_z_i-2
kk((i-1)*(n_z_i-1)+j)=-(-b_gauss*p((i-1)*(n_z_i-1)+j-1)-a_gauss*p((i-1)*(n_z_i-1)+j+1)-e_gauss*p((i)*(n_z_i-1)+j)-&
                    d_gauss*p((i-1)*(n_z_i-1)+j))
end do

j=n_z_i-1
kk((i-1)*(n_z_i-1)+j)=-(-b_gauss*p((i-1)*(n_z_i-1)+j-1)-e_gauss*p((i)*(n_z_i-1)+j)-d_gauss*p((i-1)*(n_z_i-1)+j))

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%Interior del dominio
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!$OMP DO  
do i=2,n_x_i-2

j=1
kk((i-1)*(n_z_i-1)+j)=-(-a_gauss*p((i-1)*(n_z_i-1)+j+1)-e_gauss*p((i)*(n_z_i-1)+j)-&
f_gauss*p((i-2)*(n_z_i-1)+j)-d_gauss*p((i-1)*(n_z_i-1)+j))

   do j=2,n_z_i-2

kk((i-1)*(n_z_i-1)+j)=-(-b_gauss*p((i-1)*(n_z_i-1)+j-1)-a_gauss*p((i-1)*(n_z_i-1)+j+1)&
-e_gauss*p((i)*(n_z_i-1)+j)-f_gauss*p((i-2)*(n_z_i-1)+j)-d_gauss*p((i-1)*(n_z_i-1)+j))

    end do

j=n_z_i-1
kk((i-1)*(n_z_i-1)+j)=-(-b_gauss*p((i-1)*(n_z_i-1)+j-1)-e_gauss*p((i)*(n_z_i-1)+j)-f_gauss*p((i-2)*(n_z_i-1)+j)-&
        d_gauss*p((i-1)*(n_z_i-1)+j))
end do
!$OMP END DO
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%x=L
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=n_x_i-1
j=1
kk((i-1)*(n_z_i-1)+j)=-(-a_gauss*p((i-1)*(n_z_i-1)+j+1)-f_gauss*p((i-2)*(n_z_i-1)+j)&
-d_gauss*p((i-1)*(n_z_i-1)+j))

do j=2,n_z_i-2
kk((i-1)*(n_z_i-1)+j)=-(-b_gauss*p((i-1)*(n_z_i-1)+j-1)-a_gauss*p((i-1)*(n_z_i-1)+j+1)-f_gauss*p((i-2)*(n_z_i-1)+j)&
-d_gauss*p((i-1)*(n_z_i-1)+j))
end do

j=n_z_i-1
kk((i-1)*(n_z_i-1)+j)=-(-b_gauss*p((i-1)*(n_z_i-1)+j-1)-f_gauss*p((i-2)*(n_z_i-1)+j)-d_gauss*p((i-1)*(n_z_i-1)+j))
!$OMP END PARALLEL
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

producto_1=dot_product(p,kk)

alpha=producto/producto_1


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%
!%Calculo del nuevo residuo
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!call daxpy(num_fila, alpha, kk, 1, vector_psi, 1)  
vector_psi=vector_psi+alpha*p
residuo_prima=residuo-alpha*kk
z_cg=inv_d*residuo_prima
producto_residuo_prima=dot_product(residuo_prima,z_cg)
beta=producto_residuo_prima/producto
p=z_cg+beta*p
residuo=residuo_prima
producto=producto_residuo_prima
iter=iter+1;
max_res=maxval(dabs(residuo_prima))

end do
do j=1,n_z_i-1
    do i=1,n_x_i-1
    psi_out(i+1,j+1)=vector_psi((i-1)*(n_z_i-1)+j)
    end do
end do


!write(*,'(e17.9)')max_res
!write(*,*)"total iteration",iter


end subroutine stream_function
