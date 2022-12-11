program main
	use TOV
	use MFT
	use CMD_Progress
    implicit none
	type (CLS_CMD_Progress)::Progress
    integer::num
!	call module_parameter()
!	call hyperon_parameter()
    call DDRMF()
    call Progress % Set(N=140,L=50)
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
        sigma_r = 0
        sigma=0
        omega=0
        rho_3=0
	    open(1,file='density.txt')
	    open(2,file='EOS.dat')
        write(2,*) "dentity  pressure  energy"
	    open(3,file='parameter.txt')
	    open(4,file='M_R.txt')
	    do num=1,140
			call Progress % Put(num,CMD_PROGRESS_ABSOLUTE)
	        rho_b=num * 0.01
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
	        Energy(num)=cal_energy()/planck**3
	        Pressure(num)=cal_pressure()/planck**3
!            Pressure1(num)=cal_pressure1()/planck**3
            call cal_TOV(num,M_neutron,R_neutron)
           100 format(1X,13f15.6)
	        write(1,100)rho_b/rho0,rho_n/rho_b,rho_p/rho_b,rho_e/rho_b,rho_u/rho_b, &
	        rho_lam/rho_b,rho_sig1/rho_b,rho_sig0/rho_b,rho_sig_1/rho_b,rho_ksi0/rho_b,rho_ksi_1/rho_b 
	        write(2,'(3e15.6)')rho_b,Pressure(num),Energy(num)
	        write(3,100)rho_b/rho0,effect_nucl/Mb,X_sig_lam,X_sig_sig,X_sig_ksi
	        write(4,'(4e15.6)')rho_b,Energy(num)*1.7827e12,M_neutron,R_neutron
 !           write(3,100)rho_b,sigma,omega,rho_3,sigma_r,Energy(num)/rho_b-Mb
    end do
	write(*,*) "done!"
    end   
