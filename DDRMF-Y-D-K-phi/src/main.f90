program main
	use TOV
	use MFT
	use CMD_Progress
    implicit none
    real(kind=8)::up_k,down_k
    character(len=32)::par
	type (CLS_CMD_Progress)::Progress
    integer::num
    call choose_delta()
    call DDRMF()
!    call sigmacut()
	call Progress % Set(N=228,L=50)
	Progress % Prefix = "calculate :"
	Progress % M = "#"
	Progress % O = "."
	    rho_n=0
	    rho_p=0
	    rho_e=0
	    rho_u=0
	    rho_lam=0
	    rho_sig1=0
	    rho_sig0=0
	    rho_sig_1=0
	    rho_ksi0=0
	    rho_ksi_1=0
        rho_delta2=0
        rho_delta1=0
        rho_delta0=0
        rho_delta_1=0
	    kf_n=0
	    kf_p=0
	    kf_e=0
	    kf_u=0
	    kf_lam=0
	    kf_sig1=0
	    kf_sig0=0
	    kf_sig_1=0
	    kf_ksi0=0
	    kf_ksi_1=0
        kf_delta2=0
        kf_delta1=0
        kf_delta0=0
        kf_delta_1=0
	
	    open(1,file='density.txt')
	    open(2,file='EOS.dat')
        write(2,*) "  dentity  pressure  energy"
	    open(3,file='parameter.txt')
	    open(4,file='M_R.txt')
	    100 format(1X,16f15.6)
	    write(*,*) "include kaon mesons?   ---y[n]"
        read(*,'(a1)') par
        if(par=='y'.or.par=='Y') then
        call kaon_meson()
        do num=1,228
			call Progress % Put(num,CMD_PROGRESS_ABSOLUTE)
	        rho_b=num * 0.005
            up_k=1
            down_k=0.0
            do while(.true.)
                rho_kaon = 0.5*(up_k+down_k)
	            alpha_n=cal_alpha(rho_b)
                call cal_che_kaon()
                if(che_e.lt.che_kaon) then
                    up_k = rho_kaon
                else
                    down_k=rho_kaon
                end if
                if((abs(che_e-che_kaon).lt.0.00001).or.(abs(up_k-down_k).lt.0.00001)) exit
                end do
                if (abs(che_e-che_kaon).gt.0.6) then
                    rho_kaon=0
                end if
                
            alpha_n=cal_alpha(rho_b)
	        rho_n=alpha_n*rho_b
	        rho_p=kf_p**3/kf_n**3 * rho_n
	        rho_e=kf_e**3/kf_n**3 * rho_n
	        rho_u=kf_u**3/kf_n**3 * rho_n
	        rho_lam=kf_lam**3/kf_n**3 * rho_n
	        rho_sig1=kf_sig1**3/kf_n**3 * rho_n
	        rho_sig0=kf_sig0**3/kf_n**3 * rho_n
	        rho_sig_1=kf_sig_1**3/kf_n**3 * rho_n
	        rho_ksi0=kf_ksi0**3/kf_n**3 * rho_n
	        rho_ksi_1=kf_ksi_1**3/kf_n**3 * rho_n
            rho_delta2=2*kf_delta2**3/(kf_n**3) * rho_n
            rho_delta1=2*kf_delta1**3/(kf_n**3) * rho_n
            rho_delta0=2*kf_delta0**3/(kf_n**3) * rho_n
            rho_delta_1=2*kf_delta_1**3/(kf_n**3) * rho_n
	
	        Energy(num)=cal_energy()/planck**3
	        Pressure(num)=cal_pressure()/planck**3
	        call cal_TOV(num,M_neutron,R_neutron)
	        write(1,100)rho_b/rho0,rho_n/rho_b,rho_p/rho_b,rho_e/rho_b,rho_u/rho_b, &
	        rho_lam/rho_b,rho_sig1/rho_b,rho_sig0/rho_b,rho_sig_1/rho_b,rho_ksi0/rho_b,rho_ksi_1/rho_b, &
            rho_delta2/rho_b,rho_delta1/rho_b,rho_delta0/rho_b,rho_delta_1/rho_b,rho_kaon/rho_b
	        write(2,'(4e15.6)')rho_b,Pressure(num),Energy(num)
	        write(3,100)rho_b/rho0,g_sig_delta,g_sig_kaon,rho_kaon,che_e,che_kaon,&
                (M_meson_kaon-g_sig_kaon*sigma)/M_meson_kaon,effect_n/M_n
	        write(4,*)rho_b,Energy(num)*1.7827e12,M_neutron,R_neutron
	        end do
        end if
	    if(par=='n'.or.par=='N') then
                rho_kaon=0
                che_kaon=0
            do num=1,228
			call Progress % Put(num,CMD_PROGRESS_ABSOLUTE)
	        rho_b=num * 0.005
            alpha_n=cal_alpha(rho_b)
	        rho_n=alpha_n*rho_b
	        rho_p=kf_p**3/kf_n**3 * rho_n
	        rho_e=kf_e**3/kf_n**3 * rho_n
	        rho_u=kf_u**3/kf_n**3 * rho_n
	        rho_lam=kf_lam**3/kf_n**3 * rho_n
	        rho_sig1=kf_sig1**3/kf_n**3 * rho_n
	        rho_sig0=kf_sig0**3/kf_n**3 * rho_n
	        rho_sig_1=kf_sig_1**3/kf_n**3 * rho_n
	        rho_ksi0=kf_ksi0**3/kf_n**3 * rho_n
	        rho_ksi_1=kf_ksi_1**3/kf_n**3 * rho_n
            rho_delta2=2*kf_delta2**3/(kf_n**3) * rho_n
            rho_delta1=2*kf_delta1**3/(kf_n**3) * rho_n
            rho_delta0=2*kf_delta0**3/(kf_n**3) * rho_n
            rho_delta_1=2*kf_delta_1**3/(kf_n**3) * rho_n
	
	        Energy(num)=cal_energy()/planck**3
	        Pressure(num)=cal_pressure()/planck**3
	        call cal_TOV(num,M_neutron,R_neutron)
            write(1,100)rho_b/rho0,rho_n/rho_b,rho_p/rho_b,rho_e/rho_b,rho_u/rho_b, &
	        rho_lam/rho_b,rho_sig1/rho_b,rho_sig0/rho_b,rho_sig_1/rho_b,rho_ksi0/rho_b,rho_ksi_1/rho_b, &
           rho_delta2/rho_b,rho_delta1/rho_b,rho_delta0/rho_b,rho_delta_1/rho_b,rho_kaon/rho_b
            write(2,'(4e15.6)')rho_b,Pressure(num),Energy(num)
	        write(3,100)rho_b/rho0,g_sig_delta/g_sigma,effect_n/M_n,effect_p/M_p
	        write(4,*)rho_b,Energy(num)*1.7827e12,M_neutron,R_neutron
	        end do
        end if
 
            write(*,*) "done!"
    end
