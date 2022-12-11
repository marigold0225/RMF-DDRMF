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
!phi
    real(kind=8)::M_meson_phi
    real(kind=8)::phi
    real(kind=8)::g_phi
    real(kind=8)::g_phi_lam
    real(kind=8)::g_phi_sig
    real(kind=8)::g_phi_ksi
    real(kind=8)::X_phi_lam
    real(kind=8)::X_phi_sig
    real(kind=8)::X_phi_ksi
!hyperon_parameter 
         real(kind=8)::X_sig_lam
         real(kind=8)::X_ome_lam
         real(kind=8)::X_rho_lam
         real(kind=8)::X_sig_sig
         real(kind=8)::X_ome_sig
         real(kind=8)::X_rho_sig
         real(kind=8)::X_sig_ksi
         real(kind=8)::X_ome_ksi
         real(kind=8)::X_rho_ksi
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

!RRDMF
         real(kind=8)::g_sigma_0
         real(kind=8)::g_omega_0
         real(kind=8)::g_rho_0
         real(kind=8)::a_sig
         real(kind=8)::a_ome
         real(kind=8)::a_rho
         real(kind=8)::b_sig
         real(kind=8)::b_ome
         real(kind=8)::b_rho
         real(kind=8)::c_sig
         real(kind=8)::c_ome
         real(kind=8)::c_rho
         real(kind=8)::d_sig
         real(kind=8)::d_ome
         real(kind=8)::d_rho
         real(kind=8)::sigma_r
         real(kind=8)::rho0
         real(kind=8)::sigma_r0

         real(kind=8)::sigma_r_n
         real(kind=8)::sigma_r_p
         real(kind=8)::sigma_r_lam
         real(kind=8)::sigma_r_sig1
         real(kind=8)::sigma_r_sig0
         real(kind=8)::sigma_r_sig_1
         real(kind=8)::sigma_r_ksi0
         real(kind=8)::sigma_r_ksi_1
contains

subroutine DDRMF()
    implicit none
    integer::par,ierr
    interactive_loop: do
        write(*,*) " -- choose a module   1[PKDD]  2[DDME1]  3[DDME2]  4[TW99]  5[DD-LZ1]  6[DDMEX]  7[DD2]"
        read(*, *,iostat=ierr) par
        if(ierr /= 0 ) then
            write(*,*) " -- Error, invalid input"
            cycle interactive_loop
        end if
        if(ierr .eq. 0) exit interactive_loop
    end do interactive_loop
    if(par.eq.1) then
        write(*,*) "---------choose PKDD----------"
        rho0 = 0.1496
    a_sig = 1.327
    a_ome = 1.3421
    a_rho = 0.1833
    b_sig = 0.4351
    b_ome = 0.3711
    c_sig = 0.6916
    c_ome = 0.6113
    d_sig = 0.6942
    d_ome = 0.7383
    g_sigma_0 = 10.73
    g_omega_0 = 13.14
    g_rho_0 = 8.5996
    M_meson_sigma = 555
    M_meson_omega = 783
    M_meson_rho = 763
!    X_sig_lam = 0.610412
!    X_sig_sig = 0.461907
!    X_sig_ksi = 0.302709
    kappa = 0
    lambda = 0
    xi = 0
    lambda_nu = 0
    sigma_zero = 37.566379
    omega_zero = 24.633646
    sigma_r0=9.620949
    end if
    if (par.eq.2) then
        write(*,*) "------choose  DDME1--------"
        rho0 = 0.152
    a_sig = 1.3854
    a_ome = 1.3879
    a_rho = 0.5008
    b_sig = 0.9781
    b_ome = 0.8525
    c_sig = 1.5342
    c_ome = 1.3566
    d_sig = 0.4661
    d_ome = 0.4957
    g_sigma_0 = 10.4434
    g_omega_0 = 12.8939
    g_rho_0 = 7.6106
    M_meson_sigma = 549.5
    M_meson_omega = 783
    M_meson_rho = 763
!    X_sig_lam = 0.608602
!    X_sig_sig = 0.457163
!    X_sig_ksi = 0.301777
    kappa = 0
    lambda = 0
    xi = 0
    lambda_nu = 0
    sigma_zero = 37.983772
    omega_zero = 24.595118
    sigma_r0=9.463733
    end if
     if (par.eq.3) then
        write(*,*) "------choose DDME2 --------"
        rho0 = 0.152
    a_sig = 1.3881
    a_ome = 1.3892
    a_rho = 0.5647
    b_sig = 1.0943
    b_ome = 0.9240
    c_sig = 1.7057
    c_ome = 1.4620
    d_sig = 0.4421
    d_ome = 0.4775
    g_sigma_0 = 10.53
    g_omega_0 = 13.01
    g_rho_0 = 7.36
    M_meson_sigma = 550
    M_meson_omega = 783
    M_meson_rho = 763
!    X_sig_lam = 0.620035
!   X_sig_sig = 0.470799
!    X_sig_ksi = 0.315064
    kappa = 0
    lambda = 0
    xi = 0
    lambda_nu = 0
    sigma_zero = 38.165262
    omega_zero = 24.800062
    sigma_r0=7.85165

    end if
     if (par.eq.4) then
        write(*,*) "------choose TW99 --------"
        rho0 = 0.153
    a_sig = 1.365469
    a_ome = 1.402488
    a_rho = 0.515
    b_sig = 0.226061
    b_ome = 0.172577
    c_sig = 0.409704
    c_ome = 0.344293
    d_sig = 0.901995
    d_ome = 0.983955
    g_sigma_0 = 10.7285
    g_omega_0 = 13.2902
    g_rho_0 = 7.32196
    M_meson_sigma = 550
    M_meson_omega = 783
    M_meson_rho = 763
!    X_sig_lam = 0.612049
!    X_sig_sig = 0.468796
!    X_sig_ksi = 0.303632
    kappa = 0
    lambda = 0
    xi = 0
    lambda_nu = 0
    sigma_zero = 38.955149
    omega_zero = 25.484705
    sigma_r0=8.266248
    end if
      if (par.eq.5) then
        write(*,*) "------choose DD-LZ1 --------"
        rho0 = 0.1581
    a_sig = 1.062748
    b_sig = 1.763627
    c_sig = 2.308928
    d_sig = 0.379957
    a_ome = 1.059181
    b_ome = 0.418273
    c_ome = 0.538663
    d_ome = 0.786649
    a_rho = 0.776095
    g_sigma_0 = 12.001429
    g_omega_0 = 14.292525
    g_rho_0 = 15.150934
    M_meson_sigma = 538.619216
    M_meson_omega = 783
    M_meson_rho = 769
    kappa = 0
    lambda = 0
    xi = 0
    lambda_nu = 0
    sigma_zero = 40.2726322
    omega_zero = 25.75736
    sigma_r0=8.115498
    end if
      if (par.eq.6) then
        write(*,*) "------choose DDMEX --------"
        rho0 = 0.152
    a_sig = 1.3970
    b_sig = 1.3350
    c_sig = 2.0671
    d_sig = 0.4016
    a_ome = 1.3936
    b_ome = 1.0191
    c_ome = 1.6060
    d_ome = 0.4556
    a_rho = 0.6202
    g_sigma_0 = 10.7076
    g_omega_0 = 13.3388
    g_rho_0 = 7.2380
    M_meson_sigma = 547.3327
    M_meson_omega = 783
    M_meson_rho = 763
    kappa = 0
    lambda = 0
    xi = 0
    lambda_nu = 0
    sigma_zero = 39.014539
    omega_zero = 25.410509
    sigma_r0=3.2519
    end if
      if (par.eq.7) then
        write(*,*) "------choose DD2 --------"
        rho0 = 0.149
    a_sig = 1.357630
    b_sig = 0.634442
    c_sig = 1.005358
    d_sig = 0.575810
    a_ome = 1.369718
    b_ome = 0.496475
    c_ome = 0.817753
    d_ome = 0.638452
    a_rho = 0.518903
    g_sigma_0 = 10.686681
    g_omega_0 = 13.342362
    g_rho_0 = 7.25388
    M_meson_sigma = 546.212459
    M_meson_omega = 783
    M_meson_rho = 763
    kappa = 0
    lambda = 0
    xi = 0
    lambda_nu = 0
    sigma_zero = 38.42188
    omega_zero = 24.915838
    sigma_r0=7.79411
    end if


  end


subroutine module_parameter()
	implicit none 
!	character(len=32)::arg
	integer::par,ierr
	
	interactive_loop: do
		write(*,*) " -- choose a module parameter   1[FSUGOLD]  2[IUFSU]  3[NL3]   4[TM1]"
		read(*, *,iostat=ierr) par
		if(ierr /= 0 ) then 
			write(*,*) " -- Error, invalid input"
			cycle interactive_loop 
		end if 
		if(ierr .eq. 0) exit interactive_loop
	end do interactive_loop
    
	if(par.eq.1) then ! FSUGOLD
	write(*,*) " -- set module parameter FSUGOLD"
	g_sigma=sqrt(112.1996)
	g_rho=sqrt(138.4701)
	g_omega=sqrt(204.5469)
	M_meson_sigma=491.5
	M_meson_omega=782.5
	M_meson_rho=763
	xi=0.06
	kappa=1.42
	lambda=0.0238
	lambda_nu=0.03
    end if
    if(par.eq.2) then	! IUFSU
	write(*,*) " -- set module parameter IUFSU"
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
    end if
end subroutine 

subroutine hyperon_parameter()
	implicit none

!hyperon potential------------------------------------------		 
         X_ome_lam= 2./3.
         X_rho_lam=0.

         X_ome_sig=2./3.
         X_rho_sig=2.

         X_ome_ksi=1./3.
         X_rho_ksi=1.

         M_meson_phi = 1020
         X_phi_lam = -sqrt(2.)/3.
         X_phi_sig = -sqrt(2.)/3.
         X_phi_ksi = -2*sqrt(2.)/3.

         U_lam= - 30
         U_sig=  30
         U_ksi= - 14
!hyperon -----------------------
         g_ome_lam=X_ome_lam*g_omega
         g_ome_sig=X_ome_sig*g_omega
         g_ome_ksi=X_ome_ksi*g_omega

         g_rho_lam=X_rho_lam * g_rho
         g_rho_sig=X_rho_sig * g_rho
         g_rho_ksi=X_rho_ksi * g_rho
         X_sig_lam=(X_ome_lam*g_omega_0*omega_zero - U_lam+sigma_r0)/sigma_zero/g_sigma_0
         X_sig_sig=(X_ome_sig*g_omega_0*omega_zero - U_sig+sigma_r0)/sigma_zero/g_sigma_0
         X_sig_ksi=(X_ome_ksi*g_omega_0*omega_zero - U_ksi+sigma_r0)/sigma_zero/g_sigma_0
         g_sig_lam= X_sig_lam*g_sigma
         g_sig_sig= X_sig_sig*g_sigma
         g_sig_ksi= X_sig_ksi*g_sigma
!         X_sig_lam = g_sig_lam/g_sigma
!         X_sig_sig = g_sig_sig/g_sigma
!         X_sig_ksi = g_sig_ksi/g_sigma
         g_phi_lam = X_phi_lam * g_omega
         g_phi_sig = X_phi_lam * g_omega
         g_phi_ksi = X_phi_ksi * g_omega

end subroutine
end module 
