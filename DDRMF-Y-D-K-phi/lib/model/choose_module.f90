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
    real(kind=8)::M_delta
    real(kind=8)::xi
    real(kind=8)::kappa
    real(kind=8)::lambda
    real(kind=8)::lambda_nu
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
         real(kind=8)::X_omega
         real(kind=8)::X_rho
         real(kind=8)::U_lam
         real(kind=8)::U_sig
         real(kind=8)::U_ksi
!         real(kind=8)::U_delta
!         real(kind=8)::U_n
         real(kind=8)::X_d
         real(kind=8)::sigma_zero
         real(kind=8)::omega_zero
!kaon meson
         real(kind=8)::X_sig_kaon
         real(kind=8)::X_ome_kaon
         real(kind=8)::X_rho_kaon
         real(kind=8)::U_kaon
         real(kind=8)::g_sig_kaon
         real(kind=8)::g_ome_kaon
         real(kind=8)::g_rho_kaon

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
         real(kind=8)::X_sig_delta
         real(kind=8)::X_ome_delta
         real(kind=8)::X_rho_delta
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

!sigma cut parameter
	 	 real(kind=8)::alpha
         real(kind=8)::beta
         real(kind=8)::f_zero
         real(kind=8)::c_sigma
         real(kind=8)::f_cut

!DDRMF parameter
         real(kind=8)::rho0
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
         real(kind=8)::sigma_r0
         real(kind=8)::sigma_r
         real(kind=8)::sigma_r_n
         real(kind=8)::sigma_r_p
         real(kind=8)::sigma_r_lam
         real(kind=8)::sigma_r_sig1
         real(kind=8)::sigma_r_sig0
         real(kind=8)::sigma_r_sig_1
         real(kind=8)::sigma_r_ksi0
         real(kind=8)::sigma_r_ksi_1
         real(kind=8)::sigma_r_delta2
         real(kind=8)::sigma_r_delta1
         real(kind=8)::sigma_r_delta0
         real(kind=8)::sigma_r_delta_1
         real(kind=8)::M_n
         real(kind=8)::M_p


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
        M_n = 939.5731
        M_p = 938.2796
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
        M_n = 939
        M_p = 939
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
        M_n = 939
        M_p = 939
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
        M_n = 939
        M_p = 939
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
        M_n = 938.9
        M_p = 938.9
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
        M_n = 939
        M_p = 939
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
        M_n = 939.56536
        M_p = 938.27203
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
!kaon meson
subroutine kaon_meson()
    implicit none
    real(kind=8)::U
    write(*,*)" --  choose the optical potential of kaon meson"
    read(*,*) U
    if(U.eq. -120) then
        write(*,*) "-- set U_kaon = -120"
        U_kaon = -120
    end if
    if(U.eq. -140) then
        write(*,*) "-- set U_kaon = -140"
        U_kaon = -140
    end if
    if(U.eq. -160) then
        write(*,*) "-- set U_kaon = -160"
        U_kaon = -160
    end if
 
    X_ome_kaon = 1./3.
    X_rho_kaon = 1
    g_phi_kaon = 4.27
    g_ome_kaon=X_ome_kaon * g_omega_0
    g_rho_kaon=X_rho_kaon * g_rho_0
    X_sig_kaon=( - U_kaon - g_ome_kaon*omega_zero+sigma_r0)/sigma_zero/g_sigma_0
    g_sig_kaon = g_sigma_0*X_sig_kaon

end subroutine


subroutine choose_delta()
    implicit none
    integer::a,b
    write(*,*) "--choose X_sig_delta---"
    read(*,*) a
    if(a.eq.1) then
        X_sig_delta= 1.05
    end if
    if(a.eq.2) then
        X_sig_delta=1.1
    end if
    if(a.eq.3) then
        X_sig_delta=1.15
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

    U_lam = - 30
    U_sig = 30
    U_ksi = - 14
    X_rho_lam=0
    X_rho_sig=2
    X_rho_ksi=1

    X_rho_delta =1
    X_ome_delta = 1.1

    X_ome_lam=2./3.
    X_ome_sig=2./3.
    X_ome_ksi=1./3.


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

         X_sig_lam = (g_omega_0*X_ome_lam*omega_zero-U_lam+sigma_r0)/sigma_zero/g_sigma_0
         X_sig_sig = (g_omega_0*X_ome_sig*omega_zero-U_sig+sigma_r0)/sigma_zero/g_sigma_0
         X_sig_ksi = (g_omega_0*X_ome_ksi*omega_zero-U_ksi+sigma_r0)/sigma_zero/g_sigma_0

         g_sig_lam= X_sig_lam * g_sigma
         g_sig_sig= X_sig_sig * g_sigma
         g_sig_ksi= X_sig_ksi * g_sigma

         g_phi_lam = X_phi_lam * g_omega
         g_phi_sig = X_phi_sig * g_omega
         g_phi_ksi = X_phi_ksi * g_omega
         g_phi_d=0
         g_phi=0

         g_ome_delta =X_ome_delta*g_omega
         g_rho_delta =X_rho_delta*g_rho
!         X_sig_delta =(X_ome_delta*g_omega_0*omega_zero-U_delta+sigma_r0)/sigma_zero/g_sigma_0
         g_sig_delta = X_sig_delta*g_sigma
!      end if
end subroutine

subroutine sym_energy_IUFSU()
	implicit none 
	real(kind=8)::L
	write(*,*) "--please set symmetry energy parameter L"
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
