program main
	use TOV
	use MFT
	use CMD_Progress
    implicit none
	type (CLS_CMD_Progress)::Progress
    real(kind=8)::kf_b,t,cs
    integer::num
	call module_parameter()
	call hyperon_parameter()
    call sigmacut()
	call Progress % Set(N=240,L=50)
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
	    open(1,file='density.txt')
	    open(2,file='EOS.dat')
        write(2,*) "  density  pressure  energy"
	    open(3,file='parameter.txt')
	    open(4,file='M_R.txt')
	    do num=1,240
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
	        Energy(num)=cal_energy()/planck**3
	        Pressure(num)=cal_pressure()/planck**3
	        call cal_TOV(num,M_neutron(num),R_neutron(num))
	        100 format(1X,13f15.6)
	        write(1,100)rho_b/rho0,rho_n/rho_b,rho_p/rho_b,rho_e/rho_b,rho_u/rho_b, &
	        rho_lam/rho_b,rho_sig1/rho_b,rho_sig0/rho_b,rho_sig_1/rho_b,rho_ksi0/rho_b,rho_ksi_1/rho_b 
	        write(2,'(3e15.6)')rho_b,Pressure(num),Energy(num)!,(Pressure(num+1)-Pressure(num))/(Energy(num+1)-Energy(num))
	        write(3,100)rho_b/rho0,sigma,effect_nucl/Mb,g_sigma*sigma/Mb,g_sig_lam/g_sigma,g_sig_sig/g_sigma,g_sig_ksi/g_sigma
	        write(4,'(4e15.6)')rho_b/rho0,Energy(num)*1.7827e12,M_neutron(num),R_neutron(num)
	    end do
	write(*,*) "done!"
    end   
