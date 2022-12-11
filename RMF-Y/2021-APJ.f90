module lepton
	implicit none
!meson coupling constants-------------------------------------------------
	real(kind=8),parameter::g_sigma=10.217
	real(kind=8),parameter::g_omega=12.868
	real(kind=8),parameter::g_rho_0=8.948
	real(kind=8),parameter::g2 = - 10.431 * 197.33
	real(kind=8),parameter::g3 = - 28.885
	real(kind=8),parameter::M_meson_sigma=508.194
	real(kind=8),parameter::M_meson_omega=782.501
	real(kind=8),parameter::M_meson_rho=763
	real(kind=8),parameter::rho_0=0.148
!meson field--------------------------------------------------------------
	real(kind=8)::sigma
	real(kind=8)::omega
	real(kind=8)::rho
!density denpend g_rho---------------------------------------------------
	real(kind=8)::g_rho_rho
	real(kind=8),parameter::alpha_rho=0.1578
	real(kind=8)::alpha_n
	real(kind=8)::S
!effect mass--------------------------------------------------------------
	real(kind=8)::effect_nucl
!particle density---------------------------------------------------------
	real(kind=8)::rho_b
	real(kind=8)::rho_n
	real(kind=8)::rho_p
	real(kind=8)::rho_e
	real(kind=8)::rho_u
!fermi energy--------------------------------------------------------------
	real(kind=8)::kf_n
	real(kind=8)::kf_p
	real(kind=8)::kf_e
	real(kind=8)::kf_u
!chemical potentials------------------------------------------------------
	real(kind=8)::che_n
	real(kind=8)::che_p
	real(kind=8)::che_e
	real(kind=8)::che_u
!sigma term---------------------------------------------------------------
	real(kind=8)::sigma_r
!neutron star-------------------------------------------------------------
	real(kind=8)::M_neutron
	real(kind=8)::R_neutron

!other constants----------------------------------------------------------
	real(kind=8),parameter::Pi=3.1415926535
	real(kind=8),parameter::planck=197.33
	real(kind=8),parameter::Mb=939.
	real(kind=8),parameter::M_e=0.511
	real(kind=8),parameter::M_u=105.7
	real(kind=8),parameter::gama=2
	real(kind=8),public,save::aa(4),xx(4)
	real(kind=8),public,save::Energy(3000),Pressure(3000)
	data xx /- 0.8611363, - 0.3398810, 0.3398810, 0.8611363/
	data aa /0.3478543, 0.6521452, 0.6521452, 0.3478543/

contains

!rho-->fermi energy-------------------------------------------------------
	real function fermi_energy(rho_baryon) result(kf)
		implicit none
		real(kind=8)::rho_baryon
		if (rho_baryon.lt.0.0) then
			kf=0.0
		else
			kf=rho_baryon * 3 * Pi**2
			kf=planck * kf**(1./3.)
		end if
		return
	end function
!sigma-->effect mass----------------------------------------------------
	subroutine cal_effect_mass()
		implicit none
		effect_nucl=Mb + g_sigma * sigma
	end subroutine
!g_rho_rho-->density-------------------------------------------------------------
	real function cal_g_rho_rho(rho_baryon)
		implicit none 
		real(kind=8)::rho_baryon
		cal_g_rho_rho=g_rho_0 * exp( - alpha_rho * (rho_baryon/rho_0 - 1))
		return
	end function
!rho_baryon--->sigma_r----------------------------------------------------------
	real function cal_sigma_r(rho_baryon,rho_neutron,rho_proton)
		implicit none
		real(kind=8)::rho_baryon,rho_neutron,rho_proton
		cal_sigma_r= - 0.5 * alpha_rho * cal_g_rho_rho(rho_baryon) * &
			(rho_proton - rho_neutron) * rho/rho_0
		return
	end function
!rho,effect_mass-->sclar density-----------------------------------------------	
	real function cal_sclar_density(rho_baryon,effect_mass)
		implicit none
		real(kind=8)::rho_baryon,effect_mass,kf_fermi
		if ((rho_baryon.gt.0.0).and.(effect_mass.gt.0.0)) then
			kf_fermi=fermi_energy(rho_baryon)
			cal_sclar_density=effect_mass/(2*Pi**2)*(kf_fermi*sqrt(kf_fermi**2 + &
				effect_mass**2)- effect_mass**2 * &
				log((kf_fermi + sqrt(kf_fermi**2 + effect_mass**2))/effect_mass))
		else
			cal_sclar_density=0.0
		end if
		return
	end function 
!sigma meson field function -->nucl
	real function cal_sigma_nucl(rho_neutron,rho_proton)
		implicit none
		real(kind=8)::rho_neutron,rho_proton,left,right,up1,up,down,effect_mass,err0,err1,err2
		up=0.0
		down= - 100.
		sigma=up
		effect_mass=Mb + g_sigma * sigma
		left=M_meson_sigma**2 * sigma + g2 * sigma**2 + g3 * sigma**3
		right= - g_sigma * (cal_sclar_density(rho_neutron,effect_mass) + &
			cal_sclar_density(rho_proton,effect_mass))
		err0=left - right
		do while(.true.)
			sigma=down
			effect_mass=Mb + sigma * g_sigma
			if (effect_mass.gt.0) then
				left=M_meson_sigma**2 * sigma + g2 * sigma**2 + g3 * sigma**3
				right= - g_sigma * (cal_sclar_density(rho_neutron,effect_mass) + &
					cal_sclar_density(rho_proton,effect_mass))
				err1=left - right
			else
				up1=down
				down=(up+down)/2
				cycle
			end if
			if (err0 * err1.lt.0) then
				exit
			else
				up=down
				down=up1
			end if
		end do
		do while(.true.)
			sigma=(up+down)/2
			effect_mass=Mb + g_sigma * sigma
			left=M_meson_sigma**2 * sigma + g2 * sigma**2 + g3 * sigma**3
			right= - g_sigma * (cal_sclar_density(rho_neutron,effect_mass) + &
				cal_sclar_density(rho_proton,effect_mass))
			err2=left-right
			if (err2 * err0.lt.0.0) then
				down=sigma
			else
				up=sigma
			end if
			if (abs(down-up).lt.0.00001) exit
		end do
		cal_sigma_nucl=(up+down)/2.
	end function
!omega meson field equation------------------------------------------------------
	subroutine cal_omega()
		implicit none
		real(kind=8)::left,right
		right=g_omega * (rho_n + rho_p)
		left=M_meson_omega**2
		omega=planck**3*right/left
	end subroutine
!rho meson field equation--------------------------------------------------------
	subroutine cal_rho()
		implicit none
		real(kind=8)::right,left
		right=g_rho_rho * (rho_p - rho_n)/2.
		left=M_meson_rho**2
		rho=planck**3*right/left
	end subroutine
!rho,effect_mass,omega,rho-->chemical potentials of baryon-----------------------
	real function chemical_neutron(rho_neutron,rho_proton,rho_baryon,effect_mass,omega_m,rho_m)
		implicit none 
		real(kind=8)::rho_neutron,rho_proton,rho_baryon,effect_mass,omega_m,rho_m
		if(rho_neutron.gt.0.0) then
			chemical_neutron=sqrt(fermi_energy(rho_neutron)**2 + effect_mass**2) + &
				g_omega * omega_m + cal_g_rho_rho(rho_baryon) * rho_m * (-0.5)  + &
				cal_sigma_r(rho_baryon,rho_neutron,rho_proton)
		else 
			chemical_neutron=0.0
		end if
		return
	end function
!proton chemical potentials ------------------------------------------------------
	real function chemical_proton(rho_proton,rho_neutron,rho_baryon,effect_mass,omega_m,rho_m)
		implicit none
		real(kind=8)::rho_proton,rho_neutron,rho_baryon,effect_mass,omega_m,rho_m
		if (rho_proton.gt.0.0) then
			chemical_proton=sqrt(fermi_energy(rho_proton)**2 + effect_mass**2) + g_omega * &
				omega_m + cal_sigma_r(rho_baryon,rho_neutron,rho_proton) + &
				cal_g_rho_rho(rho_baryon) * rho_m * 0.5
		else
			chemical_proton=0.0
		end if
		return
	end function
!neutron partical density------------------------------------------------------------
	real function cal_alpha_n(rho_baryon)
		implicit none
		real(kind=8)::rho_baryon,alpha1,up,down,up1,down1,diff
		up=1.
		down=0.0
		alpha_n=(up+down)/2.
		do while(.true.)
			rho_n=alpha_n * rho_baryon
			kf_n=fermi_energy(rho_n)
			rho_p=rho_baryon - rho_n
			sigma=cal_sigma_nucl(rho_n,rho_p)
			g_rho_rho=cal_g_rho_rho(rho_baryon)
			call cal_omega()
			call cal_rho()
			sigma_r=cal_sigma_r(rho_baryon,rho_n,rho_p)
			rho_e=rho_p
			call cal_effect_mass()
			che_n=chemical_neutron(rho_n,rho_p,rho_baryon,effect_nucl,omega,rho)
			che_p=chemical_proton(rho_p,rho_n,rho_baryon,effect_nucl,omega,rho)
			if (rho_e.gt.0.0) then
				che_e=sqrt(fermi_energy(rho_e) * fermi_energy(rho_e) + M_e**2)
			else
				che_e=0.0
				kf_e=0.0
			end if
			if (che_e.lt.M_u) then
				che_u=0.0
				kf_u=0.0
				rho_u=0.0
				diff=che_n - che_e - che_p
				if (abs(up-down).lt.0.0001.or.abs(diff).lt.0.0001) then
					exit
				else
					if(diff.gt.0.0) then 
						up=alpha_n
						alpha_n=(up+down)/2.
					else
						down=alpha_n
						alpha_n=(up+down)/2.
					end if
				end if
			else
				up1=1.
				down1=0.
				alpha1=(up1+down1)/2
				do while(.true.)
					rho_e=rho_p * alpha1
					rho_u=rho_p - rho_e
					che_e=sqrt(fermi_energy(rho_e) * fermi_energy(rho_e) + M_e**2)
					if (rho_u.gt.0.0) then
						che_u=sqrt(fermi_energy(rho_u) * fermi_energy(rho_u) + M_u**2)
					else
						che_u=0.0
					end if
					if(abs(up1-down1).lt.0.001.or.abs(che_e-che_u).lt.0.001) then
						exit
					else
						if(che_e.lt.che_u) then
							down1=alpha1
							alpha1=(up1+down1)/2.
						else
							up1=alpha1
							alpha1=(up1+down1)/2.
						end if
					end if
				end do
				diff=che_n - che_p - che_e
				if(abs(up-down).lt.0.0001.or.abs(diff).lt.0.0001) then
					exit
				else
					if(diff.gt.0.) then
						up=alpha_n
						alpha_n=(down + up)/2.
					else
						down=alpha_n
						alpha_n=(up+down)/2.
					end if
				end if
			end if
		end do		
		cal_alpha_n=alpha_n
		kf_n=fermi_energy(rho_n)
		kf_p=fermi_energy(rho_p)
		kf_e=fermi_energy(rho_e)
		kf_u=fermi_energy(rho_u)
	end function
!calculate energy---------------------------------------------------------------
	real function cal_energy()
		implicit none
		real(kind=8)::E
		E=gama/(2*Pi)**3 * (int_baryon(kf_n,effect_nucl) + & 
			int_baryon(kf_p,effect_nucl))
		E=E + (int_lepton(kf_e,M_e) + int_lepton(kf_u,M_u))/Pi**2
		E=E + 0.5 * M_meson_sigma**2 * sigma**2 + g2*sigma**3/3. + g3*sigma**4/4. + &
			M_meson_omega**2 * omega**2/2. + M_meson_rho**2 * rho**2/2.
		cal_energy=E
		return 
	end function
!calculate Pressure----------------------------------------------------------------
	real function cal_pressure()
		implicit none
		real(kind=8)::P
		P=gama/3/(2*Pi)**3*(int_P_baryon(kf_n,effect_nucl) + &
			int_P_baryon(kf_p,effect_nucl))
		P=P+(int_P_lepton(kf_e,M_e)+int_P_lepton(kf_u,M_u))/3/Pi**2 + sigma_r * rho_b * planck**3
		P=P - M_meson_sigma**2 * sigma**2/2. - g2*sigma**3/3. - g3*sigma**4/4.
		P=P + M_meson_omega**2 * omega**2/2. + M_meson_rho**2 * rho**2/2.
		cal_pressure=P
	end function		
! gauss integrate baryon energy------------------------------------------
                 real function int_baryon(dk,effect_mass) 
					  implicit none
					  real(kind=8):: dk,effect_mass
					  real(kind=8):: up,down,t,k,re,f
					  integer:: i
					  f(k)=4*Pi*k**2*sqrt(k**2+effect_mass**2)
					  if (dk.gt.0.) then
						  re=0.
						  up=dk
						  down=0.
						  do i=1,4,1
							  t=(up-down)*xx(i)/2+(up+down)/2
							  re=re+aa(i)*f(t)
						  end do
						  re=re*(up-down)/2.
						  int_baryon=re
					  else
						  int_baryon=0.
					  end if
					  return
			     end function
!gauss integrate lepton energy
                real function int_lepton(dk,effect_mass)
						implicit none
						real(kind=8):: dk,effect_mass
						real(kind=8):: up,down,t,k,re,f
						integer:: i
						f(k)=k**2*sqrt(effect_mass**2+k**2)
						if (dk.gt.0.) then
							re=0.
							up=dk
							down=0.
							do i=1,4,1
								t=(up-down)*xx(i)/2+(up+down)/2
								re=re+aa(i)*f(t)
							end do
							re=re*(up-down)/2.
							int_lepton=re
						else
							int_lepton=0.
						end if
						return
				end function				
!gauss integrate baryon press
                 real function int_P_baryon(dk,effect_mass) 
						 implicit none
						 real(kind=8):: dk,effect_mass
						 real(kind=8):: up,down,re,t,k,f
						 integer:: i
						 f(k)=4*Pi*k**4/sqrt(effect_mass**2+k**2)
						 if (dk.gt.0.) then
							 re=0.
							 up=dk
							 down=0.
							 do i=1,4,1
								 t=(up-down)*xx(i)/2+(up+down)/2.
								 re=re+aa(i)*f(t)
							 end do
							 re=re*(up-down)/2
							 int_P_baryon=re
						 else
							 int_P_baryon=0.
						 end if
						 return
				 end function				 
!gauss integrate lepton press
                 real function int_P_lepton(dk,effect_mass)
						 implicit none
						 real(kind=8):: dk,effect_mass
						 real(kind=8):: up,down,re,t,k,f
						 integer ::i
						 f(k)=k**4/sqrt(effect_mass**2+k**2)
						 if (dk.gt.0) then
							 re=0.
							 up=dk
							 down=0.
							 do i=1,4,1
								 t=(up-down)*xx(i)/2+(up+down)/2
								 re=re+aa(i)*f(t)
							 end do
							 re=re*(up-down)/2
						     int_P_lepton=re
						 else
							 int_P_lepton=0.
						 end if
						 return
				 end function
!caclulate TOV equation----------------------------------------------------------------------------
                 subroutine TOV(n,M,R)
						 implicit none
						 real(kind=8):: temp_P,temp_E,dp,dr,M,R
                         real(kind=8),parameter::M_sun=1.9891e33
						 real(kind=8),parameter::G=6.673e-8
						 real(kind=8),parameter::unit_E=1.7827e12
						 real(kind=8),parameter::unit_P=1.59162e33
						 real(kind=8),parameter::v_light=2.998e10
						 integer ::n,loop
						 loop=n
						 temp_E=Energy(loop)
						 temp_P=Pressure(loop)
						 R=1.
						 dR=1.
						 M=4.*Pi*R**3*temp_E*unit_E/3.
						 do while(.true.)
							 dp= - (G/R**2 * M * temp_E * unit_E) * (1 + temp_P/temp_E)
							 dp=dp/(1 - 2 * G * M / (R * v_light**2)) * (1 + 4 * Pi * R**3 * temp_P*unit_E / M)
							 temp_P=temp_P + dp * dR/unit_P
							 if(temp_P.lt.1e-4) exit
							 do while(.true.)
								 if(loop.eq.0) exit
								 if(temp_P.gt.Pressure(loop)) then
									 exit
								 else
									 loop=loop-1
								 end if
							 end do
							 if (loop.eq.0) then
								 temp_E=(Energy(loop+1)/Pressure(loop+1))*temp_P
							 else
								 temp_E=Energy(loop) + (Energy(loop+1) - Energy(loop))/(Pressure(loop+1) &
										 - Pressure(loop)) * (temp_P-Pressure(loop))
							 end if
							 R=R+dR
							 M=M + 4 * Pi * R**2 * temp_E * unit_E * dR
						 end do
						 M=M/M_sun
						 R=R/100000
				 end subroutine

 end
program main
	use lepton
	implicit none
	character (len=32)::arg
	real(kind=8)::rho_1,h,kf_b
	integer::num,n1
	call get_command_argument (1,arg)
	read (arg,*) rho_1
	call get_command_argument (2,arg)
	read (arg,*) n1
	call get_command_argument (3,arg)
	read (arg,*) h
	rho_n=0
	rho_p=0
	rho_e=0
	rho_u=0
	kf_n=0
	kf_p=0
	kf_e=0
	kf_u=0
	open(1,file='bar_energy.csv')
	open(2,file='EOS.csv')
	open(3,file='g_rho_rho.csv')
	open(4,file='M_R.csv')
	open(5,file='S.csv')
	write(1,*) "density(fm^-3),BE(MeV/fm-3)"
	write(2,*) "Energy(MeV/fm-3),Pressure(MeV/fm-3)"
	write(3,*) "g_rho_rho,density"
	write(4,*) "density,Energy,M,R"
	write(5,*) "density,symmetry-Energy"
	do num=1,n1
		rho_b=rho_1 + num * h
		kf_b=rho_b * 3 * Pi**2
		kf_b=(kf_b/2.)**(1./3.)
		kf_b=planck * kf_b
		g_rho_rho=cal_g_rho_rho(rho_b)
		alpha_n=cal_alpha_n(rho_b)
		rho_n=rho_b * alpha_n
		rho_p=kf_p**3/kf_n**3 * rho_n
		rho_e=kf_e**3/kf_n**3 * rho_n
		rho_u=kf_u**3/kf_n**3 * rho_u
		S=kf_b**2/(6 * sqrt(kf_b**2 + effect_nucl**2)) + g_rho_rho**2 * rho_b * planck**3 /(8 * M_meson_rho**2)
		Energy(num)=cal_energy()/planck**3
		Pressure(num)=cal_pressure()/planck**3
		call TOV(num,M_neutron,R_neutron)
		100 format(1X,13(g15.6,','))
		write(1,100)rho_b,Energy(num)/rho_b - Mb
		write(2,100)Energy(num),Pressure(num)
		write(3,100)rho_b,sigma_r
		write(4,100)rho_b,Energy(num)*1.7827e12,M_neutron,R_neutron
		write(5,100)rho_b,S
	end do
	end
	
