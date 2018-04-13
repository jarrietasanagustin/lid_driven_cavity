program main
!------------------------------------
!Este programa es una version mejorada
!para resolver el problema de bioconveccion
!con un foco luminoso.
!
!El metodo numerico consta de un esquema
!Runge-Kutta de tres etapas
!para la evolucion temporal y diferencias finitas de
!segundo orden para la discretizacion espacial
!
!           IMEDEA, Septiembre  2015
!
!   Autores: Jorge Arrieta
!
!------------------------------------
use definiciones
implicit none
integer(4)::computing_time,computing_time_final,clock_rate,clock_max

integer(4)::ke,i_time
real(8)::comp_time,suma_n
real(8)::v_s,l,delta_rho,g,n_0,r_cell,pi,diff_cell
real(8)::nu
i_write=600
!i_print=1000
n_0=1.0d6!in cells/ml
g=980.0d0!in cm/s^2
delta_rho=0.05d0!delta_rho/rho
v_s=0.5d0*78d-4!0.14d0*78d-4!58.0d-4;
l=1.0d0
r_cell=5.0d-4!in cm
pi=4.0d0*atan(1.0d0)
nu=1.0d-2
diff_cell=3.9d-4!5575.0d-8

gr=1.0d6!n_0*g*delta_rho*(4.0d0/3.0d0)*pi*(r_cell)**3*(l**3.0d0)/nu**2
sc=1.0d0
v_fall=sqrt(delta_rho*g*((4.0d0/3.0d0)*pi*(r_cell)**3)*n_0*l)
alpha_phot=v_s/v_fall


call read_data

call initialize

call mesh

call initial_conditions

ke=1
i_time=1

!ra=1.5d0*ra
!re=re!/10.0d0
!pe=pe!/10.0d0
write(*,'(f25.10)')gr,sc,v_fall
write(*,*)"altura del foco",z_foco
!stop
!call mkl_set_num_threads( 2)
call system_clock(computing_time,clock_rate,clock_max)
bucle_temporal:do

call proceso_iterativo(ke,i_time)

!----------------------------------------
!Escritura en pantalla del pasto temporal
!----------------------------------------




if(i_write==600)then
suma_n=maxval(dabs(omega_new-omega_old))!sum(n_new)/((nx_i+1)*(nz_i+1))
write(*,*)i_time,i_write
write(*,'(e17.9)')dabs(suma_n-1.0d0)
call system_clock(computing_time_final,clock_rate,clock_max)
write(*,*)"elapsed iteration time=",real(computing_time_final-computing_time)/(real(clock_rate))

i_write=0
computing_time=computing_time_final
end if

i_time=i_time+1
i_write=i_write+1
!call output_data(ke,i_time)
if(maxval(dabs(omega_new-omega_old))<=1.0d-11)then

exit
end if


n_old(1:nx_i+1,1:nz_i+1)=n_new(1:nx_i+1,1:nz_i+1)
omega_old(1:nx_i+1,1:nz_i+1)=omega_new(1:nx_i+1,1:nz_i+1)
psi_old(1:nx_i+1,1:nz_i+1)=psi_new(1:nx_i+1,1:nz_i+1)
u_old(1:nx_i+1,1:nz_i+1)=u_new(1:nx_i+1,1:nz_i+1)
w_old(1:nx_i+1,1:nz_i+1)=w_new(1:nx_i+1,1:nz_i+1)


ke=1








end do bucle_temporal
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!A partir de este instante se apaga la luz
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!w_photo=0.0d0
!u_photo=0.0d0
!bucle_temporal_off:do

!call proceso_iterativo(ke,i_time)

!----------------------------------------
!Escritura en pantalla del pasto temporal
!----------------------------------------

!write(*,*)i_time,i_write


!if(i_write==200)then
!call output_data(ke,i_time)
!i_write=0
!end if

!i_time=i_time+1
!i_write=i_write+1




!n_old(1:nx_i+1,1:nz_i+1)=n_new(1:nx_i+1,1:nz_i+1)
!omega_old(1:nx_i+1,1:nz_i+1)=omega_new(1:nx_i+1,1:nz_i+1)
!psi_old(1:nx_i+1,1:nz_i+1)=psi_new(1:nx_i+1,1:nz_i+1)
!u_old(1:nx_i+1,1:nz_i+1)=u_new(1:nx_i+1,1:nz_i+1)
!w_old(1:nx_i+1,1:nz_i+1)=w_new(1:nx_i+1,1:nz_i+1)


!ke=1


!if(i_time>=nt_i+1)then
!exit
!end if

!end do bucle_temporal_off
call cpu_time(comp_time)

write(*,*)comp_time
contains

!--------------------------------------
!Declaracion de las variables del problema
!lectura del fichero de datos para generar vectores
!y definir los parametros del problema
!--------------------------------------

subroutine read_data

use definiciones
implicit none

namelist /mesh/beta_z,h_z,nz,nz_i,tau,h_x,nx,nx_i
namelist /tiempo/t_max,nt,nt_i,ne_i
namelist /tolerance/eps
!------------------------------------
!sch=schidmt number
!gr=grashof number
!------------------------------------
namelist /bioconvec/sc,re,pe,ra,pr
namelist /phototaxis/i_max,x_foco,z_foco

open(unit=1,file='data.dat',status='old')

read(unit=1,nml=mesh)
read(unit=1,nml=tiempo)
read(unit=1,nml=tolerance)
!read(unit=1,nml=bioconvec)
read(unit=1,nml=phototaxis)

close(1)


end subroutine read_data


!------------------------------------
!Declaracion del tamanho de los vectores y
!matrices que se utilizaran en el problema
!Inicializacion de todas las variables a 0
!------------------------------------

subroutine initialize

use definiciones
implicit none
integer(4)::i,j,k

num_files=(nt_i/i_write)
write(*,*)num_files
!stop


allocate(time(num_files+1),z(nz_i+1),x(nx_i+1))
allocate(u_old(nx_i+1,nz_i+1),w_old(nx_i+1,nz_i+1),psi_old(nx_i+1,nz_i+1),omega_old(nx_i+1,nz_i+1),&
        n_old(nx_i+1,nz_i+1))
allocate(u_new(nx_i+1,nz_i+1),w_new(nx_i+1,nz_i+1),psi_new(nx_i+1,nz_i+1),omega_new(nx_i+1,nz_i+1),&
        n_new(nx_i+1,nz_i+1))

!allocate(u_star(nx_i+1,nz_i+1),w_star(nx_i+1,nz_i+1),omega_star(nx_i+1,nz_i+1),psi_star(nx_i+1,nz_i+1),&
!        n_star(nx_i+1,nz_i+1))
allocate(u_aprox(nx_i+1,nz_i+1),w_aprox(nx_i+1,nz_i+1),omega_aprox(nx_i+1,nz_i+1),&
        n_aprox(nx_i+1,nz_i+1))
allocate(u_photo(nx_i+1,nz_i+1),w_photo(nx_i+1,nz_i+1),intensity(nx_i+1,nz_i+1))
allocate(u_aprox_1(nx_i+1,nz_i+1),w_aprox_1(nx_i+1,nz_i+1),omega_aprox_1(nx_i+1,nz_i+1),&
n_aprox_1(nx_i+1,nz_i+1))
allocate(u_aprox_2(nx_i+1,nz_i+1),w_aprox_2(nx_i+1,nz_i+1),omega_aprox_2(nx_i+1,nz_i+1),&
n_aprox_2(nx_i+1,nz_i+1))


do i=1,nx_i+1
    do j=1,nz_i+1
intensity(i,j)=0.0d0
u_photo(i,j)=0.0d0
w_photo(i,j)=0.0d0
    end do
end do


do i=1,num_files+1
time(i)=0.0d0
end do

do i=1,nz_i+1
z(i)=0.0d0
end do

do i=1,nx_i+1
x(i)=0.0d0
end do




end subroutine initialize

!-----------------------------------
!Definicion de la malla del problema
!Se define una malla uniforme en el dominio
![0,L]x[0,L]
!-----------------------------------
subroutine mesh

use definiciones
implicit none
integer(4)::i,j,k
real(8)::sigma_l
integer(4)::n_elements_i
real(8)::n_elements,h_x_prima,h_z_prima
real(8)::alpha,beta!Parameters of the mesh
real(8)::h_x_segunda,h_z_segunda
real(8)::n_focos,delta_foco
integer(4)::n_focos_i
real(8),dimension(nz_i+1,nx_i+1)::d_uphoto
!beta_z=1.02d0
n_focos_i=8
n_focos=8.0d0
delta_x=h_x/nx

do i=1,nx_i+1
x(i)=(i-1)*delta_x
end do

delta_z=(h_z)/nz

do j=1,nz_i+1
z(j)=(j-1)*delta_z
end do

delta_t=t_max/nt


do k=1,num_files+1
    time(k)=real(i_write)*(k-1)*delta_t
end do

open(unit=10,file='time.dat',status='unknown')
open(unit=11,file='x.dat',status='unknown')
open(unit=13,file='z.dat',status='unknown')


write(10,'(e17.9)')time
write(11,'(e17.9)')x
write(13,'(e17.9)')z


close(10)
close(11)
close(12)
close(13)
close(14)


!-----------------------------------------
!Definicion del campo de intensidades luminicas
!y del campo de velocidades debida a la fototaxis
!-----------------------------------------

!z_foco=0.8d0/13.0d0
sigma_l=0.667d0*1.0d+3/(1.0d+4*1.0d0)!2*0.125d0/15.0d0!*2.0d0/4.0d0!sqrt(0.5d0)
do i=1,nx_i+1
    do j=1,nz_i+1

        intensity(i,j)=exp(-((x(i)-x_foco)**2+(z(j)-z_foco)**2)/(2.0d0*(sigma_l)**2))
        u_photo(i,j)=alpha_phot*sqrt((x(i)-x_foco)**2+(z(j)-z_foco)**2)/(519d-4)*&
                        (-exp(0.5d0))/sigma_l*(x(i)-x_foco)*exp(-((x(i)-x_foco)**2&
                    +(z(j)-z_foco)**2)/(2.0d0*sigma_l**2))
        w_photo(i,j)=-alpha_phot*sqrt((x(i)-x_foco)**2+(z(j)-z_foco)**2)/(519d-4)*&
                        exp(0.5d0)/sigma_l*(z(j)-z_foco)*exp(-((x(i)-x_foco)**2&
                    +(z(j)-z_foco)**2)/(2.0*sigma_l**2))
!        intensity(i,j)=exp(-((x(i)-(x_foco-0.0725d0))**2+(z(j)-z_foco)**2)/(2.0d0*(sigma_l)**2))+&
!                        exp(-((x(i)-x_foco+0.0725d0)**2+(z(j)-z_foco)**2)/(2.0d0*(sigma_l)**2))
!        u_photo(i,j)=-exp(0.5d0)/sigma_l*(x(i)-(x_foco-0.0725d0))*exp(-((x(i)-(x_foco-0.0725d0))**2&
!                    +(z(j)-z_foco)**2)/(2.0d0*sigma_l**2))&
!                      -exp(0.5d0)/sigma_l*(x(i)-(x_foco+0.0725d0))*exp(-((x(i)-(x_foco+0.0725d0))**2&
!                    +(z(j)-z_foco)**2)/(2.0d0*sigma_l**2))
!        w_photo(i,j)=-exp(0.5d0)/sigma_l*(z(j)-z_foco)*exp(-((x(i)-(x_foco-0.0725d0))**2&
!                    +(z(j)-z_foco)**2)/(2.0*sigma_l**2))&
!                        -exp(0.5d0)/sigma_l*(z(j)-z_foco)*exp(-((x(i)-(x_foco+0.0725d0))**2&
!                    +(z(j)-z_foco)**2)/(2.0*sigma_l**2))
        end do
end do

do i=2,nx_i
    do j=2,nz_i

!  d_uphoto(j,i)=(w_photo(i,j+1)-w_photo(i,j-1))/(2*delta_z)+(u_photo(i-1,j)-u_photo(i+1,j))/(2*delta_x)
   d_uphoto(j,i)=(w_photo(i,j+1)-w_photo(i,j-1))/(2*delta_z)
!        intensity(i,j)=exp(-((x(i)-(x_foco-0.0725d0))**2+(z(j)-z_foco)**2)/(2.0d0*(sigma_l)**2))+&
!                        exp(-((x(i)-x_foco+0.0725d0)**2+(z(j)-z_foco)**2)/(2.0d0*(sigma_l)**2))
!        u_photo(i,j)=-exp(0.5d0)/sigma_l*(x(i)-(x_foco-0.0725d0))*exp(-((x(i)-(x_foco-0.0725d0))**2&
!                    +(z(j)-z_foco)**2)/(2.0d0*sigma_l**2))&
!                      -exp(0.5d0)/sigma_l*(x(i)-(x_foco+0.0725d0))*exp(-((x(i)-(x_foco+0.0725d0))**2&
!                    +(z(j)-z_foco)**2)/(2.0d0*sigma_l**2))
!        w_photo(i,j)=-exp(0.5d0)/sigma_l*(z(j)-z_foco)*exp(-((x(i)-(x_foco-0.0725d0))**2&
!                    +(z(j)-z_foco)**2)/(2.0*sigma_l**2))&
!                        -exp(0.5d0)/sigma_l*(z(j)-z_foco)*exp(-((x(i)-(x_foco+0.0725d0))**2&
!                    +(z(j)-z_foco)**2)/(2.0*sigma_l**2))
        end do
end do
open(unit=15,file='intensity_n_0_1d_6.dat',status='unknown')
open(unit=16,file='u_photo_n_0_1d_6.dat',status='unknown')
open(unit=17,file='w_photo_n_0_1d_6.dat',status='unknown')
open(unit=44,file='d_uphoto.dat',status='unknown')

write(15,'(129e22.12E4,2x)')(intensity(i,1:nz_i+1),i=1,nx_i+1)
write(16,'(129e22.12E4,2x)')(u_photo(i,1:nz_i+1),i=1,nx_i+1)
write(17,'(129e22.12E4,2x)')(w_photo(i,1:nz_i+1),i=1,nx_i+1)
write(44,'(129e22.12E4,2x)')(d_uphoto(i,1:nx_i+1),i=1,nz_i+1)
close(15)
close(16)
close(17)
close(44)
!stop
deallocate(intensity,time)
end subroutine mesh


!-----------------------------------------
!Condiciones inciales, corresponden con
!fluido en reposo cuyo movimiento es iniciado
!por la conveccion natural en el interior de la
!cavidad. Por tanto no hay vorticidad y la funcion
!de corriente tiene valor constante
!-----------------------------------------

subroutine initial_conditions

use definiciones
implicit none
integer(4)::i,j,k

k=1

do i=1,nx_i+1
    do j=1,nz_i+1
        u_old(i,j)=0.0d0
        w_old(i,j)=0.0d0
        omega_old(i,j)=0.0d0
        psi_old(i,j)=0.0d0
!        n_old(i,j)=1.0d0
    end do
end do
call random_number(n_old)
n_old=1.0d0!+1.0d-3*n_old


end subroutine initial_conditions



subroutine first_stage(k,k_t)

use definiciones
implicit none
integer(4)::i,j,k,k_t



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!Primera etapa del esquema RK3
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!$OMP PARALLEL PRIVATE(i,j)  NUM_THREADS(6)
!!$OMP DO
do i=2,nx_i

    do j=2,nz_i

    k_omega(1,1)=(-u_old(i,j)*(omega_old(i+1,j)-omega_old(i-1,j))/(2.0d0*delta_x)-&
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w_old(i,j)*(omega_old(i,j+1)-omega_old(i,j-1))/(2.0d0*delta_z)+&
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
(omega_old(i+1,j)-2.0d0*omega_old(i,j)+omega_old(i-1,j))/(sqrt(gr)*delta_x**2)+&
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        (omega_old(i,j+1)-2.0d0*omega_old(i,j)+omega_old(i,j-1))/(sqrt(gr)*delta_z**2))

!    k_n(1,1)=-(u_old(i,j)+u_photo(i,j))*(n_old(i+1,j)-n_old(i-1,j))/(2.0d0*delta_x)-&
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!(w_old(i,j)+w_photo(i,j))*(n_old(i,j+1)-n_old(i,j-1))/(2.0d0*delta_z)+&
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!(n_old(i+1,j)-2.0d0*n_old(i,j)+n_old(i-1,j))/(sc*sqrt(gr)*(delta_x**2))+&
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!        (n_old(i,j+1)-2.0d0*n_old(i,j)+n_old(i,j-1))/(sc*sqrt(gr)*(delta_z**2))-&
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!        (n_old(i,j))*(u_photo(i+1,j)-u_photo(i-1,j))/(2.0d0*delta_x)-&
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                (n_old(i,j))*(w_photo(i,j+1)-w_photo(i,j-1))/(2.0d0*delta_z)



omega_aprox(i,j)=omega_old(i,j)+delta_t/2.0d0*k_omega(1,1)



    end do

end do
!!$OMP END DO
!!$OMP END PARALLEL
end subroutine first_stage

subroutine second_stage(k,k_t)
use definiciones
implicit none
integer(4)::i,j,k,k_t


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!Segunda etapa del esquema RK3
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

do i=2,nx_i

    do j=2,nz_i

k_omega(1,1)=-u_old(i,j)*(omega_old(i+1,j)-omega_old(i-1,j))/(2.0d0*delta_x)-&
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w_old(i,j)*(omega_old(i,j+1)-omega_old(i,j-1))/(2.0d0*delta_z)+&
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
(omega_old(i+1,j)-2.0d0*omega_old(i,j)+omega_old(i-1,j))/(sqrt(gr)*delta_x**2)+&
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
(omega_old(i,j+1)-2.0d0*omega_old(i,j)+omega_old(i,j-1))/(sqrt(gr)*delta_z**2)

k_omega(2,1)=-u_aprox(i,j)*(omega_aprox(i+1,j)-omega_aprox(i-1,j))/(2.0d0*delta_x)-&
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w_aprox(i,j)*(omega_aprox(i,j+1)-omega_aprox(i,j-1))/(2.0d0*delta_z)+&
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
(omega_aprox(i+1,j)-2.0d0*omega_aprox(i,j)+omega_aprox(i-1,j))/(sqrt(gr)*delta_x**2)+&
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
(omega_aprox(i,j+1)-2.0d0*omega_aprox(i,j)+omega_aprox(i,j-1))/(sqrt(gr)*delta_z**2)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






!omega_aprox_1(i,j)=omega_old(i,j)+delta_t/2.0d0*(k_omega(2,1))
!n_aprox_1(i,j)=n_old(i,j)+delta_t/2.0d0*(k_n(2,1))

omega_aprox_1(i,j)=omega_old(i,j)+delta_t*(-k_omega(1,1)+2.0d0*k_omega(2,1))


    end do

end do

end subroutine second_stage


!----------------------------------------
!Subrrutina para iterar en (x,z) para resolver
! en cada paso temporal el campo de velocidades
! funcion de corriente y vorticidad
!----------------------------------------

subroutine proceso_iterativo(k,k_t)

use definiciones
implicit none

integer(4)::iter,i,j,k,k_t
integer(4)::computing_time,computing_time_final,clock_rate,clock_max
integer(4)::n_itera
real(8)::a_top,b_top,c_top,a_bot,b_bot,c_bot
real(8)::a_left,b_left,c_left,a_right,b_right,c_right
real(8)::suma_n
!-----------------------------------
!Para la primera iteracion
!usamos el paso de la iteracion anterior
!-----------------------------------
!call system_clock(computing_time,clock_rate,clock_max)


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!Primera etapa del esquema Runge-Kutta
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
call first_stage(k,k_t)


    call stream_function(omega_aprox,psi_old,psi_new,nx_i,nz_i)


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!Aplicamos condiciones de contorno antes de la segunda
!etapa
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

call bc_first_stage(k,k_t)



    call second_stage(k,k_t)

   call stream_function(omega_aprox_1,psi_old,psi_new,nx_i,nz_i)

    call bc_second_stage(k,k_t)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!Tercera etapa Runge-Kutta
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

call third_stage(k,k_t)

call stream_function(omega_new,psi_old,psi_new,nx_i,nz_i)



call bc_third_stage(k,k_t)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!Cuarta etapa Runge-Kutta
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!call fourth_stage(k,k_t)

!call stream_function(omega_new,psi_old,psi_new,nx_i,nz_i)

!call bc_fourth_stage(k,k_t)




!call system_clock(computing_time_final,clock_rate,clock_max)
!write(*,*)"elapsed iteration time=",real(computing_time_final-computing_time)/(real(clock_rate))


end subroutine proceso_iterativo

subroutine output_data(k,k_t)

use definiciones
implicit none

integer(4)::i,j,k,k_t
character*8 counter

 write(counter,'(I0)')k_t



open(unit=21,file='./data_psi/psi_0'//trim(counter)//'.dat',status='unknown')
!open(unit=22,file='./data_n_n_0_1d_6/n_0'//trim(counter)//'.dat',status='unknown')
open(unit=23,file='./data_u/u_0'//trim(counter)//'.dat',status='unknown')
open(unit=24,file='./data_w/w_0'//trim(counter)//'.dat',status='unknown')
open(unit=25,file='./data_vort/omega_0'//trim(counter)//'.dat',status='unknown')


write(21,'(81e17.9,1x)')(psi_old(i,1:nz_i+1),i=1,nx_i+1)
!write(22,'(129e17.9,1x)')(n_old(i,1:nz_i+1),i=1,nx_i+1)
write(23,'(81e17.9,1x)')(u_old(i,1:nz_i+1),i=1,nx_i+1)
write(24,'(81e17.9,1x)')(w_old(i,1:nz_i+1),i=1,nx_i+1)
write(25,'(81e17.9,1x)')(omega_old(i,1:nz_i+1),i=1,nx_i+1)


close(21)
!close(22)
close(23)
close(24)
close(25)

end subroutine output_data
end program main
