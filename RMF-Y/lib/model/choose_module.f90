!*************************************************************************
!      > File Name: choose_module.f90
!      > Author: marigold
!      > Mail: mflovelky418@gmail.com 
!      > Created Time: 2022年01月07日 星期五 17时31分13秒
! ************************************************************************
module choose_module
	implicit none
! FSUGold module constants---------------------------    
    real(kind=8)::g_sigma
    real(kind=8)::g_rho
    real(kind=8)::g_omega
    real(kind=8)::M_meson_sigma
    real(kind=8)::M_meson_omega
    real(kind=8)::M_meson_rho
    real(kind=8)::xi
    real(kind=8)::kappa
    real(kind=8)::lambda
    real(kind=8)::lambda_nu
    real(kind=8)::rho0
!phi meson
    real(kind=8)::M_meson_phi
    real(kind=8)::g_phi_lam
    real(kind=8)::g_phi_sig
    real(kind=8)::g_phi_ksi
    real(kind=8)::g_phi
    real(kind=8)::X_phi_lam
    real(kind=8)::X_phi_sig
    real(kind=8)::X_phi_ksi
    real(kind=8)::phi

!hyperon_parameter 
         real(kind=8)::X_omega
         real(kind=8)::X_rho
         real(kind=8)::U_lam
         real(kind=8)::U_sig
         real(kind=8)::U_ksi
         real(kind=8)::sigma_zero
         real(kind=8)::omega_zero
!hyperon meson coupli:s----------
         real(kind=8)::g_ome_lam
         real(kind=8)::g_ome_sig
         real(kind=8)::g_ome_ksi
         real(kind=8)::g_rho_lam
         real(kind=8)::g_rho_sig
         real(kind=8)::g_rho_ksi
         real(kind=8)::g_sig_lam
         real(kind=8)::g_sig_sig  	 
         real(kind=8)::g_sig_ksi

!sigma cut parameter
	 	 real(kind=8)::alpha
         real(kind=8)::beta
         real(kind=8)::f_zero
         real(kind=8)::c_sigma
         real(kind=8)::f_cut


contains

subroutine module_parameter()
	implicit none 
!	character(len=32)::arg
	integer::par,ierr
	
	interactive_loop: do
		write(*,*) " -- choose a module parameter   1[FSUGOLD]  2[IUFSU]  3[NL3]   4[TM1]   5[Bigapple]"
		read(*, *,iostat=ierr) par
		if(ierr /= 0 ) then 
			write(*,*) " -- Error, invalid input"
			cycle interactive_loop 
		end if 
		if(ierr .eq. 0) exit interactive_loop
	end do interactive_loop

	if(par.eq.1) then ! FSUGOLD
	write(*,*) " -- set module parameter FSUGOLD"
    rho0=0.149
	g_sigma=10.5924
	g_rho=11.7673
	g_omega=14.302
	M_meson_sigma=491.5
	M_meson_omega=782.5
	M_meson_rho=763
	xi=0.06
	kappa=1.42
	lambda=0.0238
	lambda_nu=0.03
    sigma_zero=34.507135
    omega_zero=20.594168
    end if
 
    if(par.eq.2) then	! IUFSU
	write(*,*) " -- set module parameter IUFSU"
    rho0=0.1545
	g_sigma=sqrt(99.4266)
	g_rho=sqrt(184.6877)
	g_omega=sqrt(169.8349)
	M_meson_sigma=491.5
	M_meson_omega=783
	M_meson_rho=763
	xi=0.03
	kappa=3.3808
	lambda=0.000296
	lambda_nu=0.046
    sigma_zero=36.7791
    omega_zero=22.5406
end if
    if(par.eq.3) then	! NL3
	write(*,*) " -- set module parameter NL3"
	g_sigma=sqrt(104.3871)
	g_rho=sqrt(79.6)
	g_omega=sqrt(165.5854)
	M_meson_sigma=508.194
	M_meson_omega=782.501
	M_meson_rho=763
	xi=0
	kappa=3.8599
	lambda= - 0.015905
	lambda_nu=0
    sigma_zero = 37.168022
    omega_zero = 23.89912605
    end if
    if(par.eq.4) then	! TM1
	write(*,*) " -- set module parameter TM1"
	g_sigma=10.0289
	g_rho=9.2644
	g_omega=12.6139
	M_meson_sigma=511.198
	M_meson_omega=783
	M_meson_rho=770
	xi=6*71.3075/g_omega**4
	kappa=197.33*7.2325*2/g_sigma**3
	lambda=6 * 0.6183/g_sigma**4
	lambda_nu=0
    sigma_zero = 34.15911495
    omega_zero = 21.729707
    end if

    if(par.eq.5) then !Bigapple
        write(*,*) " -- set module parameter Bigapple"
        g_sigma = 9.6699
        g_omega = 12.316
        g_rho = 14.1618
        M_meson_sigma = 429.73
        M_meson_omega = 782.5
        M_meson_rho   = 763.
        kappa = 5.20326
        lambda = - 0.021739
        xi = 0.00070
        lambda_nu = 0.047471
        sigma_zero =37.9535 
        omega_zero =23.895
    end if

end subroutine 

subroutine hyperon_parameter()
	implicit none
!hyperon potential------------------------------------------
         X_omega=2./3.
         X_rho=2.
         U_lam= - 30
         U_sig=  30
         U_ksi= - 14
        M_meson_phi = 1020
        X_phi_lam = - sqrt(2.)/3.
        X_phi_sig = - sqrt(2.)/3.
        X_phi_ksi = - 2*sqrt(2.)/3
        g_phi = 0

!hyperon -----------------------
         g_ome_lam=X_omega*g_omega
         g_ome_sig=X_omega*g_omega
         g_ome_ksi=g_omega/3.
         g_rho_lam=0.0
         g_rho_sig=X_rho*g_rho
         g_rho_ksi=g_rho
         g_sig_lam=(g_ome_lam*omega_zero - U_lam)/sigma_zero
         g_sig_sig=(g_ome_sig*omega_zero - U_sig)/sigma_zero
         g_sig_ksi=(g_ome_ksi*omega_zero - U_ksi)/sigma_zero
end subroutine


!sigma cut paramter ---------------
subroutine sigmacut()
	implicit none
	real(kind=8)::scut
	write(*,*) "--please set sigma cut parameter"
	read (*,*) scut
    if(scut.eq.0) then
        write(*,*) "--without sigma-cut"
        alpha=0
        beta=0
        f_zero=0
        c_sigma=0
		f_cut=f_zero + c_sigma * (1-f_zero)
    end if 
	if(abs(scut- 0.1).lt.0.000001) then
		write(*,*) "--set sigma_cut parameter = 0.1"
		alpha=3.7487767271173e8
		beta=120
		f_zero=0.392595
		c_sigma=0.1
		f_cut=f_zero + c_sigma * (1-f_zero)
	end if 
	if (abs(scut - 0.15).lt.0.000001) then 
		write(*,*) "--set sigma cut parameter = 0.15"
		alpha=3.7487767271173e8
		beta=120
		f_zero=0.392595
		c_sigma=0.15
		f_cut=f_zero + c_sigma * (1-f_zero)
	end if 
	if (abs(scut - 0.2).lt.0.000001) then 
		write(*,*) "--set sigma cut parameter = 0.2"
		alpha=3.7487767271173e8
		beta=120
		f_zero=0.392595
		c_sigma=0.2
		f_cut=f_zero + c_sigma * (1-f_zero)
	end if 
	end subroutine

end module 
