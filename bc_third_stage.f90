subroutine bc_third_stage(k,k_t)

use definiciones
implicit none

integer(4)::iter,i,j,k,k_t
real(8)::a_top,b_top,c_top,a_bot,b_bot,c_bot
real(8)::a_left,b_left,c_left,a_right,b_right,c_right


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!Calculo de velocidades y aplicacion de condiciones de
!contorno
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!----------------------------------------
!Condiciones de contorno para la vorticidad
!---------------------------------------

do j=1,nz_i+1
i=1
!omega_new(i,j)=(psi_new(i+2,j)*(x(i+1)-x(i))**3-psi_new(i+1,j)*(x(i+2)-x(i))**3)/&
!(0.5d0*((x(i+2)-x(i))**3*(x(i+1)-x(i))**2-(x(i+1)-x(i))**3*(x(i+2)-x(i))**2))
!omega_aprox_2(i,j)=-2.0d0*psi_new(i+1,j)*(h_x**2)/delta_x**2
!omega_aprox_2(i,j)=(psi_new(i+2,j)-8.0d0*psi_new(i+1,j))/(2.0d0*delta_x**2)
omega_new(i,j)=(psi_new(i+2,j)-8.0d0*psi_new(i+1,j))/(2.0d0*delta_x**2)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=nx_i+1
!omega_new(i,j)=(psi_new(i-2,j)*(x(i)-x(i-1))**3-psi_new(i-1,j)*(x(i)-x(i-2))**3)/&
!(0.5d0*((x(i)-x(i-2))**3*(x(i)-x(i-1))**2-(x(i)-x(i-1))**3*(x(i)-x(i-2))**2))
!omega_aprox_2(i,j)=-2.0d0*psi_new(i-1,j)*(h_x**2)/delta_x**2
!omega_aprox_2(i,j)=(psi_new(i-2,j)-8.0d0*psi_new(i-1,j))/(2.0d0*delta_x**2)
omega_new(i,j)=(psi_new(i-2,j)-8.0d0*psi_new(i-1,j))/(2.0d0*delta_x**2)
end do

do i=1,nx_i+1
j=1
!omega_new(i,j)=(psi_new(i,j+2)*(z(j+1)-z(j))**3-psi_new(i,j+1)*(z(j+2)-z(j))**3)/&
!(0.5d0*((z(j+2)-z(j))**3*(z(j+1)-z(j))**2-(z(j+1)-z(j))**3*(z(j+2)-z(j))**2))
!omega_aprox_2(i,j)=-2.0d0*psi_new(i,j+1)*(h_z**2)/delta_z**2
!omega_aprox_2(i,j)=(psi_new(i,j+2)-8.0d0*psi_new(i,j+1))/(2.0d0*delta_z**2)
omega_new(i,j)=(psi_new(i,j+2)-8.0d0*psi_new(i,j+1))/(2.0d0*delta_z**2)
j=nz_i+1
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!omega_new(i,j)=(psi_new(i,j-2)*(z(j)-z(j-1))**3-psi_new(i,j-1)*(z(j)-z(j-2))**3)/&
!(0.5d0*((z(j)-z(j-2))**3*(z(j)-z(j-1))**2-(z(j)-z(j-1))**3*(z(j)-z(j-2))**2))
!omega_aprox_2(i,j)=-2.0d0*psi_new(i,j-1)*(h_z**2)/delta_z**2
!omega_aprox_2(i,j)=(psi_new(i,j-2)-8.0d0*psi_new(i,j-1))/(2.0d0*delta_z**2)
omega_new(i,j)=(psi_new(i,j-2)-8.0d0*psi_new(i,j-1)-6.0d0*delta_z)/(2.0d0*delta_z**2)
end do
do j=1,nz_i
i=1
!u_aprox_2(i,j)=0.0d0
!w_aprox_2(i,j)=0.0d0
u_new(i,j)=0.0d0
w_new(i,j)=0.0d0
i=nx_i+1
!u_aprox_2(i,j)=0.0d0
!w_aprox_2(i,j)=0.0d0
u_new(i,j)=0.0d0
w_new(i,j)=0.0d0
end do
do i=1,nx_i+1
j=1
!u_aprox_2(i,j)=0.0d0
!w_aprox_2(i,j)=0.0d0
u_new(i,j)=0.0d0
w_new(i,j)=0.0d0
j=nz_i+1
!u_aprox_2(i,j)=0.0d0
!w_aprox_2(i,j)=0.0d0
u_new(i,j)=1.0d0
w_new(i,j)=0.0d0
end do

!--------------------------------------
!Velocidad en el interior del dominio
!--------------------------------------

do i=2,nx_i
do j=2,nz_i
!u_aprox_2(i,j)=(psi_new(i,j+1)-psi_new(i,j-1))/(2.0d0*delta_z)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!w_aprox_2(i,j)=-(psi_new(i+1,j)-psi_new(i-1,j))/(2.0d0*delta_x)
u_new(i,j)=(psi_new(i,j+1)-psi_new(i,j-1))/(2.0d0*delta_z)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w_new(i,j)=-(psi_new(i+1,j)-psi_new(i-1,j))/(2.0d0*delta_x)
end do
end do


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!Fin condiciones de contorno
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end subroutine bc_third_stage
