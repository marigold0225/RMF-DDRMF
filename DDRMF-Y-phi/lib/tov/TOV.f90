!*************************************************************************
!      > File Name: TOV.f90
!      > Author: marigold
!      > Mail: mflovelky418@gmail.com 
!      > Created Time: 2022年01月07日 星期五 17时41分53秒
! ************************************************************************
!caclulate TOV equation----------------------------------------------------------------------------
module TOV 
	use gobal_constants
	implicit none 
contains
subroutine cal_TOV(n,M,R)
        real(kind=8):: temp_P,temp_E,dp,dr,M,R
        real(kind=8),parameter::M_sun=1.9891e33
        real(kind=8),parameter::G=6.673e-8
        real(kind=8),parameter::unit_E=1.7827e12
        real(kind=8),parameter::unit_P=1.59162e33
        real(kind=8),parameter::v_light=3.0e10
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
end module 
