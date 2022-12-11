!*************************************************************************
!      > File Name: MFT.f90
!      > Author: marigold
!      > Mail: mflovelky418@gmail.com 
!      > Created Time: 2022年01月07日 星期五 17时39分10秒
! ************************************************************************
module MFT 
	use choose_module
	use gobal_constants
    implicit none
contains
!rho------------------------------------>fermi energy---------------------------------------------------------------------------------
                real function fermi_energy(rho_baryon) result(kf)
                        implicit none
                        real(kind=8)::rho_baryon
                        if (rho_baryon.lt.0.0) then
                            kf=0.0
                        else
                            kf=rho_baryon * 3 * Pi**2
                            kf=planck*kf**(1./3.)
                        end if
                        return
                end function
!sigma---------------------------------->effect mass(baryon)---------------------------------------------------------------------------
                subroutine cal_effect_mass() 
                        implicit none
                        effect_nucl=Mb - g_sigma*sigma
                        effect_lam=M_lam - g_sig_lam*sigma
                        effect_sig=M_sig - g_sig_sig*sigma
                        effect_ksi=M_ksi - g_sig_ksi*sigma
                end subroutine
!rho,effect_mass------------------------>sigma*g_sigma_hyperon------------------------------------------------------------------------
                real function sigma_middle(rho_baryon,effect_mass) 
                        implicit none
                        real(kind=8)::k_fermi,effect_mass,rho_baryon
                        if ((rho_baryon.gt.0.0).and.(effect_mass.gt.0.0)) then
                            k_fermi=fermi_energy(rho_baryon)
                            sigma_middle=effect_mass/(2*Pi**2)*(k_fermi * sqrt(k_fermi**2 + effect_mass**2) - effect_mass**2 * &
                                log((k_fermi + sqrt(k_fermi**2 + effect_mass**2))/effect_mass))
                        else
                            sigma_middle=0.0
                        end if
                        return
                end function
!sigma field equation self_consistent------------------------------------------------------------------
                real function cal_sigma()
                        implicit none
                        real(kind=8)::right,left
                        right=g_sigma*(sigma_middle(rho_n,effect_nucl)+sigma_middle(rho_p,effect_nucl))
                        right=right+g_sig_lam*(sigma_middle(rho_lam,effect_lam))
                        right=right+g_sig_sig*(sigma_middle(rho_sig1,effect_sig)+sigma_middle(rho_sig0,effect_sig)+ &
                                sigma_middle(rho_sig_1,effect_sig))
                        right=right+g_sig_ksi*(sigma_middle(rho_ksi0,effect_ksi)+sigma_middle(rho_ksi_1,effect_ksi))
                        right=right-alpha*g_sigma*beta*exp(beta*(g_sigma*sigma/Mb-f_cut))/Mb/(1+exp(beta*(g_sigma*sigma/Mb-f_cut)))
                        left=M_meson_sigma**2+kappa*g_sigma**3*sigma/2.+lambda*g_sigma**4*sigma**2/6.
                        cal_sigma=right/left
                        return
                end function
!omega field equation self_consistent----------------------------------------------------------
                real function cal_omega()
                        implicit none
                        real(kind=8)::right,left
                        right=g_omega*(rho_n+rho_p)+g_ome_lam*(rho_lam)+g_ome_sig*(rho_sig1+rho_sig0+rho_sig_1)
                        right=right+g_ome_ksi*(rho_ksi0+rho_ksi_1)
                        left=M_meson_omega**2 + xi *g_omega**4 *omega**2/6. + 2 * lambda_nu*g_rho**2 * g_omega**2*rho_3**2
                        cal_omega=planck**3*right/left
                        return
                end function
!rho field equation self_consistent------------------------------------------------------------
                real function cal_rho_3()
                        implicit none
                        real(kind=8)::right,left
                        right=g_rho*iso_n*rho_n + g_rho*iso_p*rho_p + g_rho_sig*iso_sig1*rho_sig1+g_rho_sig*iso_sig_1*rho_sig_1
                        right=right+g_rho_ksi*iso_ksi0*rho_ksi0+g_rho_ksi*iso_ksi_1*rho_ksi_1
                        left=M_meson_rho**2 + 2 * lambda_nu * g_rho**2 * g_omega**2 * omega**2
                        cal_rho_3=planck**3*right/left
                        return
                end function

!phi meson field
                real function cal_phi()
                    implicit none
                    real(kind=8)::right,left
                    right = g_phi_lam*rho_lam+g_phi_sig*(rho_sig1+rho_sig0+rho_sig_1)
                    right = right+ g_phi_ksi*(rho_ksi0+rho_ksi_1)
                    left = M_meson_phi**2
                    cal_phi = planck**3*right/left
                end function


!rho,effect_mass,omega,rho_3----------->chemical potentials of baryon---------------------------------
                 real function chemical_baryon(rho_baryon,effect_mass,g_ome_coup,omega_m,g_rho_coup,iso_3,rho_3_m) result(u)
                        implicit none
                        real(kind=8)::g_ome_coup,g_rho_coup,effect_mass,iso_3,rho_baryon
                        real(kind=8)::omega_m,rho_3_m
                        if (rho_baryon.gt.0.0) then
                            u=sqrt(fermi_energy(rho_baryon)**2 + effect_mass**2) + g_ome_coup*omega_m + g_rho_coup*iso_3*rho_3_m
                        else
                            u=0.0
                        end if
                        return
                 end function
!chemical potentials------------------->fermi energy
                 real function chemical_fermi(che,effect_mass,g_ome_coup,omega_m,g_rho_coup,iso_3,rho_3_m,g_phi_c,phi_c,che_min)
                         implicit none
                         real(kind=8)::che,effect_mass,g_ome_coup,g_rho_coup,iso_3,che_min,g_phi_c,phi_c
                         real(kind=8)::middle,rho_3_m,omega_m
                         if(che.gt.che_min) then
                             middle=che - g_ome_coup*omega_m - g_rho_coup*iso_3*rho_3_m-g_phi_c*phi_c
                             middle=middle**2
                             middle=middle - effect_mass**2
                             middle=sqrt(middle)
                             chemical_fermi=middle
                         else
                             chemical_fermi=0.0
                         end if
                         return
                 end function
!chemical potentials------------------->fermi energy
 subroutine chemical_kf_rho()
         implicit none
 kf_lam=chemical_fermi(che_lam,effect_lam,g_ome_lam,omega,g_rho_lam,iso_lam,rho_3,g_phi_lam,phi,che_lam_min)
 kf_sig1=chemical_fermi(che_sig1,effect_sig,g_ome_sig,omega,g_rho_sig,iso_sig1,rho_3,g_phi_sig,phi,che_sig1_min)
 kf_sig0=chemical_fermi(che_sig0,effect_sig,g_ome_sig,omega,g_rho_sig,iso_sig0,rho_3,g_phi_sig,phi,che_sig0_min)
 kf_sig_1=chemical_fermi(che_sig_1,effect_sig,g_ome_sig,omega,g_rho_sig,iso_sig_1,rho_3,g_phi_sig,phi,che_sig_1_min)
 kf_ksi0=chemical_fermi(che_ksi0,effect_ksi,g_ome_ksi,omega,g_rho_ksi,iso_ksi0,rho_3,g_phi_ksi,phi,che_ksi0_min)
 kf_ksi_1=chemical_fermi(che_ksi_1,effect_ksi,g_ome_ksi,omega,g_rho_ksi,iso_ksi_1,rho_3,g_phi_ksi,phi,che_ksi_1_min)
            if (che_e.gt.M_e) then
                kf_e=sqrt(che_e**2-M_e**2)
            else
                kf_e=0.0
            end if
            if (che_u.gt.M_u) then
                kf_u=sqrt(che_u**2-M_u**2)
            else
                kf_u=0.0
            end if
            rho_lam=kf_lam**3/kf_n**3 * rho_n
            rho_sig1=kf_sig1**3/kf_n**3 * rho_n
            rho_sig0=kf_sig0**3/kf_n**3 * rho_n
            rho_sig_1=kf_sig_1**3/kf_n**3 * rho_n
            rho_ksi0=kf_ksi0**3/kf_n**3 * rho_n
            rho_ksi_1=kf_ksi_1**3/kf_n**3 * rho_n
            rho_e=kf_e**3/kf_n**3 * rho_n
            rho_u=kf_u**3/kf_n**3 * rho_n
    end subroutine
!calculate che_min------------------------------------------------------
                 subroutine cal_che_min()
                         implicit none
                         che_lam=che_n
                         che_sig0=che_n
                         che_ksi0=che_n
                         che_sig1=che_p
                         che_sig_1=2*che_n - che_p
                         che_ksi_1=2*che_n - che_p
                         che_e=che_n-che_p
                         che_u=che_e
                         che_lam_min=effect_lam + omega*g_ome_lam + rho_3*g_rho_lam*iso_lam+g_phi_lam*phi
                         che_sig1_min=effect_sig + omega*g_ome_sig + rho_3*g_rho_sig*iso_sig1+g_phi_sig*phi
                         che_sig0_min=effect_sig + omega*g_ome_sig + rho_3*g_rho_sig*iso_sig0+g_phi_sig*phi
                         che_sig_1_min=effect_sig + omega*g_ome_sig + rho_3*g_rho_sig*iso_sig_1+g_phi_sig*phi
                         che_ksi0_min=effect_ksi + omega*g_ome_ksi + rho_3*g_rho_ksi*iso_ksi0+g_phi_ksi*phi
                         che_ksi_1_min=effect_ksi + omega*g_ome_ksi + rho_3*g_rho_ksi*iso_ksi_1+g_phi_ksi*phi
                         if(che_lam.lt.che_lam_min) che_lam=0.0
                         if(che_sig1.lt.che_sig1_min) che_sig1=0.0
                         if(che_sig0.lt.che_sig0_min) che_sig0=0.0
                         if(che_sig_1.lt.che_sig_1_min) che_sig_1=0.0
                         if(che_ksi0.lt.che_ksi0_min) che_ksi0=0.0
                 end subroutine	
!calculate beta equilibruim --------------------------------------------
                 subroutine betaequi(rho_bar)
                         implicit none
                         real(kind=8)::up_n,down_n,up_p,down_p,rho_bar
                         up_n=1.
                         down_n=0.
                         call cal_effect_mass()
                         do while(.true.)
                             alpha_n=0.5*(up_n+down_n)
                             rho_n=rho_bar*alpha_n
                             kf_n=fermi_energy(rho_n)
                             che_n=chemical_baryon(rho_n,effect_nucl,g_omega,omega,g_rho,iso_n,rho_3)
                             up_p=1.75
                             down_p=0.0
                             do while(.true.)
                                 alpha_p=0.5 * (up_p+down_p)
                                 rho_p=(1-alpha_n)*alpha_p*rho_bar
                                 kf_p=fermi_energy(rho_p)
                                 che_p=chemical_baryon(rho_p,effect_nucl,g_omega,omega,g_rho,iso_p,rho_3)
                                 call cal_che_min()
                                 call chemical_kf_rho()
                                 if((rho_p+rho_sig1-rho_e-rho_u-rho_sig_1-rho_ksi_1).lt.0.0) then
                                     down_p=alpha_p
                                 else
                                     up_p=alpha_p
                                 end if
                                 if(abs(down_p-up_p).lt.0.00001) exit
                             end do
                             if((rho_bar-rho_n-rho_p-rho_lam-rho_sig1-rho_sig0-rho_sig_1-rho_ksi0-rho_ksi_1).lt.0.0) then
                                 up_n=alpha_n
                             else
                                 down_n=alpha_n
                             end if
                             if(abs(down_n-up_n).lt.0.00001) exit
                         end do
                 end subroutine
!calculate alpha -------------------------------------------------------
                 real function cal_alpha(rho_bar)
                         implicit none
                         real(kind=8)::temp_sigma,temp_omega,temp_rho_3,temp_phi,rho_bar
                         real(kind=8)::up_sigma,down_sigma,up_omega,down_omega,up_rho_3,down_rho_3
                         real(kind=8)::up_phi,down_phi
                         integer::i
                         i=1
                         call betaequi(rho_bar)
                         temp_sigma=cal_sigma()
                         temp_omega=cal_omega()
                         temp_rho_3=cal_rho_3()
                         temp_phi=cal_phi()
                         call rank_initia(up_sigma,down_sigma,temp_sigma,sigma)
                         call rank_initia(up_omega,down_omega,temp_omega,omega)
                         call rank_initia(up_rho_3,down_rho_3,temp_rho_3,rho_3)
                         call rank_initia(up_phi,down_phi,temp_phi,phi)
                         do while(.true.)
                             sigma=0.5*(up_sigma+down_sigma)
                             omega=0.5*(up_omega+down_omega)
                             rho_3=0.5*(up_rho_3+down_rho_3)
                             phi=0.5*(up_phi+down_phi)
                             call betaequi(rho_bar)
                             temp_sigma=cal_sigma()
                             temp_omega=cal_omega()
                             temp_rho_3=cal_rho_3()
                             temp_phi=cal_phi()
                             call rank_centre(up_sigma,down_sigma,temp_sigma,sigma)
                             call rank_centre(up_omega,down_omega,temp_omega,omega)
                             call rank_centre(up_rho_3,down_rho_3,temp_rho_3,rho_3)
                             call rank_centre(up_phi,down_phi,temp_phi,phi)
                             if(abs(up_rho_3-down_rho_3).lt.0.00001) then
                                 if(abs(up_sigma-down_sigma).lt.0.00001) then
                                     if(abs(up_omega-down_omega).lt.0.00001) then
                                         if(abs(up_phi-down_phi).lt.0.00001) then
                                         if (i.le.6) then
                                         i=i+1
                                         call betaequi(rho_bar)
                                         temp_sigma=cal_sigma()
                                         temp_omega=cal_omega()
                                         temp_rho_3=cal_rho_3()
                                         temp_phi=cal_phi()
                                         call rank_initia(up_sigma,down_sigma,temp_sigma,sigma)
                                         call rank_initia(up_omega,down_omega,temp_omega,omega)
                                         call rank_initia(up_rho_3,down_rho_3,temp_rho_3,rho_3)
                                         call rank_initia(up_phi,down_phi,temp_phi,phi)
                                         cycle
                                         end if
                                     call betaequi(rho_bar)
                                     temp_sigma=cal_sigma()
                                     temp_omega=cal_omega()
                                     temp_rho_3=cal_rho_3()
                                     temp_phi=cal_phi()
                                     sigma=0.5*(sigma+temp_sigma)
                                     omega=0.5*(omega+temp_omega)
                                     rho_3=0.5*(rho_3+temp_rho_3)
                                     phi=0.5*(phi+temp_phi)
                                     call betaequi(rho_bar)
                                     cal_alpha=alpha_n
                                         exit
                                     end if
                                     end if
                                 end if
                             end if
                         end do
                         kf_n=fermi_energy(rho_n)
                         kf_p=fermi_energy(rho_p)
                         kf_lam=fermi_energy(rho_lam)
                         kf_sig1=fermi_energy(rho_sig1)
                         kf_sig0=fermi_energy(rho_sig0)
                         kf_sig_1=fermi_energy(rho_sig_1)
                         kf_ksi0=fermi_energy(rho_ksi0)
                         kf_ksi_1=fermi_energy(rho_ksi_1)
                         kf_e=fermi_energy(rho_e)
                         kf_u=fermi_energy(rho_u)
                 end function
!calculate energy------------------------------------------------------
                 real function cal_energy() 
                         implicit none
                         real(kind=8)::E
                         E=gama/(2*Pi)**3*(int_baryon(kf_n,effect_nucl)+int_baryon(kf_p,effect_nucl)+int_baryon(kf_lam,effect_lam))
                         E=E+gama/(2*Pi)**3*(int_baryon(kf_sig1,effect_sig)+ &
							 int_baryon(kf_sig0,effect_sig)+int_baryon(kf_sig_1,effect_sig))
                         E=E+gama/(2*Pi)**3*(int_baryon(kf_ksi0,effect_ksi)+int_baryon(kf_ksi_1,effect_ksi))
                         E=E+(int_lepton(kf_e,M_e)+int_lepton(kf_u,M_u))/Pi**2+ &
							 kappa*g_sigma**3*sigma**3/6.+lambda*g_sigma**4*sigma**4/24.
                         E=E+M_meson_sigma**2*sigma**2/2.+M_meson_omega**2*omega**2/2.+M_meson_rho**2*rho_3**2/2.
                         E=E+ xi * g_omega**4*omega**4/8. + 3 * lambda_nu * g_rho**2 * g_omega**2 * omega**2 * rho_3**2
                         E=E+M_meson_phi**2*phi**2/2.
                         E=E+ alpha*log(1+exp(beta*(g_sigma*sigma/Mb-f_cut)))
                         cal_energy=E
                 end function
!calculate pressure-------------------------------------------------------
                 real function cal_pressure()
                         implicit none
                         real(kind=8)::P
                         P=gama/3/(2*Pi)**3*(int_P_baryon(kf_n,effect_nucl)+int_P_baryon(kf_p,effect_nucl))
                         P=P+gama/3/(2*Pi)**3*(int_P_baryon(kf_lam,effect_lam)+int_P_baryon(kf_sig1,effect_sig))
                         P=P+gama/3/(2*Pi)**3*(int_P_baryon(kf_sig0,effect_sig)+int_P_baryon(kf_sig_1,effect_sig))
                         P=P+gama/3/(2*Pi)**3*(int_P_baryon(kf_ksi0,effect_ksi)+int_P_baryon(kf_ksi_1,effect_ksi))
                         P=P+(int_P_lepton(kf_e,M_e)+int_P_lepton(kf_u,M_u))/3/Pi**2
                         P=P+M_meson_omega**2*omega**2/2 + M_meson_rho**2*rho_3**2/2 - M_meson_sigma**2*sigma**2/2 
                         P=P - kappa*g_sigma**3*sigma**3/6. - lambda*g_sigma**4*sigma**4/24.
                         P=P+xi*g_omega**4*omega**4/24. + lambda_nu * g_rho**2 * g_omega**2 * omega**2 * rho_3**2
                         P=P+M_meson_phi**2*phi**2/2.
                         P=P - alpha*log(1+exp(beta*(g_sigma*sigma/Mb-f_cut)))
                         cal_pressure=P
                 end function


!rank1-------------------------------------------------------------------
                 subroutine rank_centre(up,down,temp,val)
                         implicit none
                         real(kind=8)::up,down,temp,val
                         if(temp.gt.val) then
                             if(temp.gt.up) then
                                 down=val
                             else
                                 up=temp
                                 down=val
                             end if
                         else
                             if(temp.lt.down) then
                                 up=val
                             else
                                 up=val
                                 down=temp
                             end if
                         end if
                 end subroutine
!rank2--------------------------------------------------------------------
                 subroutine rank_initia(up,down,temp,val)
                         implicit none
                         real(kind=8)::up,down,temp,val
                         if(temp.gt.val) then
                             up=temp
                             down=val
                         else
                             up=val
                             down=temp
                         end if
                 end subroutine

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
end module 
