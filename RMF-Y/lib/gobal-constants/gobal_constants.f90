!*************************************************************************
!      > File Name: gobal_constants.f90
!      > Author: marigold
!      > Mail: mflovelky418@gmail.com 
!      > Created Time: 2022年01月07日 星期五 17时45分29秒
! ************************************************************************
module gobal_constants 
    implicit none
!meson field var and effect mass--------------------------------------------
         real(kind=8)::sigma
         real(kind=8)::omega
         real(kind=8)::rho_3
         real(kind=8)::effect_nucl
         real(kind=8)::effect_lam
         real(kind=8)::effect_sig
         real(kind=8)::effect_ksi
!particle density-----------------------------------------------------
         real(kind=8)::rho_b
         real(kind=8)::rho_n
         real(kind=8)::rho_p
         real(kind=8)::rho_lam
         real(kind=8)::rho_sig1
         real(kind=8)::rho_sig0
         real(kind=8)::rho_sig_1
         real(kind=8)::rho_ksi0
         real(kind=8)::rho_ksi_1
         real(kind=8)::rho_e
         real(kind=8)::rho_u
!fermi energy------------------------------------------------------------
         real(kind=8),save::kf_n
         real(kind=8),save::kf_p
         real(kind=8),save::kf_lam
         real(kind=8),save::kf_sig1
         real(kind=8),save::kf_sig0
         real(kind=8),save::kf_sig_1
         real(kind=8),save::kf_ksi0
         real(kind=8),save::kf_ksi_1
         real(kind=8),save::kf_e
         real(kind=8),save::kf_u
!chemical potentials and chemical min------------------------------------
         real(kind=8),save::che_n
         real(kind=8),save::che_p
         real(kind=8),save::che_lam
         real(kind=8),save::che_sig1
         real(kind=8),save::che_sig0
         real(kind=8),save::che_sig_1
         real(kind=8),save::che_ksi0
         real(kind=8),save::che_ksi_1
         real(kind=8),save::che_e
         real(kind=8),save::che_u		 
         real(kind=8),save::che_lam_min
         real(kind=8),save::che_sig1_min
         real(kind=8),save::che_sig0_min
         real(kind=8),save::che_sig_1_min
         real(kind=8),save::che_ksi0_min
         real(kind=8),save::che_ksi_1_min
!isopin 3rd--------------------------------------------------- 
         real(kind=8),parameter::iso_n= - 0.5
         real(kind=8),parameter::iso_p= 0.5
         real(kind=8),parameter::iso_lam=0.0
         real(kind=8),parameter::iso_sig1=1
         real(kind=8),parameter::iso_sig0=0.0
         real(kind=8),parameter::iso_sig_1= - 1
         real(kind=8),parameter::iso_ksi0= 0.5
         real(kind=8),parameter::iso_ksi_1= - 0.5
!population rho_n,rho_p
         real(kind=8),save::alpha_n
         real(kind=8),save::alpha_p
!neutron star ----------------------------------------------------
         real(kind=8),public,save::M_neutron(400)
         real(kind=8),public,save::R_neutron(400)
!other constants--------------------------------------------------		 
         real(kind=8),parameter::Pi=3.1415926535
         real(kind=8),parameter::planck=197.33
         real(kind=8),parameter::Mb=939
         real(kind=8),parameter::M_sig=1193
         real(kind=8),parameter::M_lam=1116
         real(kind=8),parameter::M_ksi=1318
         real(kind=8),parameter::M_e=0.511
         real(kind=8),parameter::M_u=105.7
         real(kind=8),parameter::gama=2
         real(kind=8),public,save::aa(4),xx(4)
         real(kind=8),public,save::Energy(3000),Pressure(3000)
         data xx /- 0.8611363,- 0.3398810,0.3398810,0.8611363/
         data  aa /0.3478543,0.6521452,0.6521452,0.3478543/
contains
end module 

