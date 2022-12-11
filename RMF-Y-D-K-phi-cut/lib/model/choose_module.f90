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
    real(kind=8)::g_phi_kaon
    real(kind=8)::g_phi
    real(kind=8)::X_phi_lam
    real(kind=8)::X_phi_sig
    real(kind=8)::X_phi_ksi
    real(kind=8)::X_phi_delta
    real(kind=8)::phi
    real(kind=8)::g_phi_d

!hyperon_parameter
         real(kind=8)::X_sig_lam
         real(kind=8)::X_sig_sig
         real(kind=8)::X_sig_ksi
         real(kind=8)::X_sig_delta
         real(kind=8)::X_sig_kaon
         real(kind=8)::X_ome_lam
         real(kind=8)::X_ome_sig
         real(kind=8)::X_ome_ksi
         real(kind=8)::X_ome_delta
         real(kind=8)::X_ome_kaon
         real(kind=8)::X_rho_lam
         real(kind=8)::X_rho_sig
         real(kind=8)::X_rho_ksi
         real(kind=8)::X_rho_delta
         real(kind=8)::X_rho_kaon
         real(kind=8)::U_lam
         real(kind=8)::U_sig
         real(kind=8)::U_ksi
         real(kind=8)::U_kaon
         real(kind=8)::U_delta
         real(kind=8)::M_delta
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
		 real(kind=8)::g_ome_delta
		 real(kind=8)::g_rho_delta
		 real(kind=8)::g_sig_delta

         real(kind=8)::g_sig_kaon
         real(kind=8)::g_ome_kaon
         real(kind=8)::g_rho_kaon

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
!	call sym_energy_IUFSU()
    rho0=0.1545
	g_sigma=9.9713
	g_rho=13.5899
	g_omega=13.0321
	M_meson_sigma=491.5
	M_meson_omega=783
	M_meson_rho=763
	xi=0.03
	kappa=3.37685
	lambda=0.000268
    sigma_zero=36.779127
    omega_zero=22.540617
	lambda_nu=0.046
    end if
    if(par.eq.3) then	! NL3
	write(*,*) " -- set module parameter NL3"
    rho0 = 0.148
	g_sigma=10.217
	g_rho=8.948
	g_omega=12.868
	M_meson_sigma=508.194
	M_meson_omega=782.501
	M_meson_rho=763
	xi=0
	kappa=3.8599
	lambda= - 0.015905
	lambda_nu=0
    sigma_zero=37.168023
    omega_zero=23.899126
    end if
    if(par.eq.4) then	! TM1
	write(*,*) " -- set module parameter TM1"
!	call sym_energy_TM1()
    rho0 = 0.145
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
    sigma_zero=34.159115
    omega_zero=21.719708
    end if
    if(par.eq.5) then !Bigapple
        write(*,*) " -- set module parameter Bigapple"
        rho0 = 0.1545
        g_sigma = 9.6699
        g_omega = 12.316
        g_rho = 14.1618
        M_meson_sigma = 429.73
        M_meson_omega = 782.5
        M_meson_rho   = 763
        kappa = 5.20326
        lambda = - 0.021739
        xi = 0.00070
        lambda_nu = 0.047471
        sigma_zero =37.9535
        omega_zero =23.895
    end if


end subroutine 

!kaon meson
subroutine kaon_meson()
    implicit none
    real(kind=8)::U
    write(*,*)" --  choose the optical potential of kaon meson"
    read(*,*) U
    if(U.eq.-120) then
        write(*,*) "-- set U_kaon = -120"
        U_kaon = -120
    end if
    if(U.eq.-140) then
        write(*,*) "-- set U_kaon = -140"
        U_kaon = -140
    end if
    if(U.eq.-160) then
        write(*,*) "-- set U_kaon = -160"
        U_kaon = -160
    end if
         X_ome_kaon=1./3.
         X_rho_kaon=1
         g_phi_kaon = 4.27
         g_ome_kaon=X_ome_kaon * g_omega
         g_rho_kaon=X_rho_kaon * g_rho
         X_sig_kaon=( -U_kaon - g_ome_kaon*omega_zero)/sigma_zero/g_sigma
         g_sig_kaon=X_sig_kaon*g_sigma

   end subroutine

subroutine choose_delta()
    implicit none
    integer::a,b
    write(*,*) "--choose X_sig_delta---"
    read(*,*) a
    if(a.eq. 1) then
        X_sig_delta= 1.05
    end if
    if(a.eq. 2) then
        X_sig_delta= 1.1
    end if
    if(a.eq.3) then
        U_delta= 1.15
    end if
    write(*,*) "--choose M_delta----"
    read(*,*) b
    if(b.eq.1112) then
        M_delta=1112
    end if
    if(b.eq.1232) then
        M_delta=1232
    end if
    if(b.eq.1352) then
        M_delta=1352
    end if
    end subroutine


subroutine hyperon_parameter()
	implicit none
!	character(len=32)::yn
!
!	write(*,*) " -- Is there a baryon octet state?  y[n]"
!	read(*,'(a1)') yn
!   if (yn=='y' .or. yn=='Y') then

!hyperon potential------------------------------------------
         X_ome_lam=2./3.
         X_ome_sig=2./3.
         X_ome_ksi=1./3.
         X_ome_delta=1.1
         X_rho_lam= 0
         X_rho_sig=2
         X_rho_ksi=1
         X_rho_delta=1
         U_lam= - 30
         U_sig=  30
         U_ksi= - 14


    M_meson_phi = 1020
    X_phi_lam = - sqrt(2.)/3.
    X_phi_sig = - sqrt(2.)/3.
    X_phi_ksi = - 2*sqrt(2.)/3
    X_phi_delta = 0

!hyperon -----------------------
         g_ome_lam=X_ome_lam*g_omega
         g_ome_sig=X_ome_sig*g_omega
         g_ome_ksi=X_ome_ksi*g_omega

         g_rho_lam=X_rho_lam * g_rho
         g_rho_sig=X_rho_sig * g_rho
         g_rho_ksi=X_rho_ksi * g_rho

         g_sig_lam=(g_ome_lam*omega_zero - U_lam)/sigma_zero
         g_sig_sig=(g_ome_sig*omega_zero - U_sig)/sigma_zero
         g_sig_ksi=(g_ome_ksi*omega_zero - U_ksi)/sigma_zero

         g_phi_lam = X_phi_lam * g_omega
         g_phi_sig = X_phi_sig * g_omega
         g_phi_ksi = X_phi_ksi * g_omega
         g_phi_d=0
         g_phi=0

         g_ome_delta=X_ome_delta*g_omega
         g_rho_delta=X_rho_delta*g_rho
         g_sig_delta=X_sig_delta*g_sigma
end subroutine

subroutine sym_energy_IUFSU()
	implicit none 
	real(kind=8)::L
	write(*,*) "--please set symmetry energy parameter L  [47.2] [50] [60] [70] [80] [90] [100] [110]"
	read(*,*) L
	if(abs(L - 47.2).lt.0.00001) then 
		write(*,*) "--set L = 47.2"
		g_rho = 13.59
		lambda_nu = 0.046
	end if 
	if(L.eq.50) then 
		write(*,*) "--set L = 50"
		g_rho = 12.8202
		lambda_nu = 0.042
	end if 
	if(L.eq.60) then 
		write(*,*) "--set L = 60"
		g_rho = 11.1893
		lambda_nu = 0.0305
	end if 
	if(L.eq.70) then 
		write(*,*) "--set L = 70"
		g_rho = 10.3150
		lambda_nu = 0.0220
	end if 
	if(L.eq.80) then 
		write(*,*) "--set L = 80"
		g_rho = 9.7537
		lambda_nu = 0.0153 
	end if 
	if(L.eq.90) then 
		write(*,*) "--set L = 90"
		g_rho = 9.3559
		lambda_nu = 0.0098
	end if 
	if(L.eq.100) then 
		write(*,*) "--set L = 100"
		g_rho = 9.0558
		lambda_nu = 0.0051
	end if 
	if(L.eq.110) then 
		write(*,*) "--set L = 110"
		g_rho = 8.8192
		lambda_nu = 0.0011
	end if 
	end subroutine

subroutine sym_energy_TM1()
	implicit none
	real(kind=8)::L
	write(*,*) "--please set symmetry parameter L"
	read(*,*) L
	if (L.eq.40) then 
		write(*,*) "--set L = 40"
		g_rho = 13.9714
		lambda_nu = 0.0429
	end if
		if(L.eq.50) then 
		write(*,*) "--set L = 50"
		g_rho = 12.2413
		lambda_nu = 0.0327
	end if 
	if(L.eq.60) then 
		write(*,*) "--set L = 60"
		g_rho = 11.2610
		lambda_nu = 0.0248
	end if 
	if(L.eq.70) then 
		write(*,*) "--set L = 70"
		g_rho = 10.6142
		lambda_nu = 0.0182
	end if 
	if(L.eq.80) then 
		write(*,*) "--set L = 80"
		g_rho = 10.1484
		lambda_nu = 0.0128 
	end if 
	if(L.eq.90) then 
		write(*,*) "--set L = 90"
		g_rho = 9.7933
		lambda_nu = 0.0080
	end if 
	if(L.eq.100) then 
		write(*,*) "--set L = 100"
		g_rho = 9.5114
		lambda_nu = 0.0039
	end if 
	if(abs(L - 110.8).lt.0.000001) then 
		write(*,*) "--set L = 110"
		g_rho = 9.2644
		lambda_nu = 0
	end if 
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
