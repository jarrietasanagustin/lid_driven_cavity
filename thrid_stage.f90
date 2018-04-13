subroutine third_stage(k,k_t)

use definiciones
implicit none
integer(4)::i,j,k,k_t


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!Tercera etapa del esquema RK3
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

k_omega(3,1)=-u_aprox_1(i,j)*(omega_aprox_1(i+1,j)-omega_aprox_1(i-1,j))/(2.0d0*delta_x)-&
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w_aprox_1(i,j)*(omega_aprox_1(i,j+1)-omega_aprox_1(i,j-1))/(2.0d0*delta_z)+&
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
(omega_aprox_1(i+1,j)-2.0d0*omega_aprox_1(i,j)+omega_aprox_1(i-1,j))/(sqrt(gr)*delta_x**2)+&
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
(omega_aprox_1(i,j+1)-2.0d0*omega_aprox_1(i,j)+omega_aprox_1(i,j-1))/(sqrt(gr)*delta_z**2)



omega_new(i,j)=omega_old(i,j)+delta_t/6.0d0*(k_omega(1,1)+4.0d0*k_omega(2,1)+k_omega(3,1))


end do

end do



end subroutine third_stage
