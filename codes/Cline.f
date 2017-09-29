      !~ 2012/05/16
      subroutine C_line
      
	use AP_comput
      implicit none
      
      include 'common_global.dat'
      
      real*8:: u_l,c_l,d_l,x_l
      real*8:: u_l1,u_l2,c_l1,c_l2
      real*8:: b_l,s0_l,nMa_l,a_l,pwet_l,rH_l,sfx_l,src_l
      real*8:: rtnewt_d,x1,x2,coef1,coef2,coef3,coef4
      real*8:: rtnewt_theta
      
      open(101,file= 'Cline.out')
      
      !~
      call bottom(3,3)
      
      if(slopm <= 1.e-6)then
        d(3)= a(3)/Bwid(3) 
      else
c        d(3)= 0.5*( -Bwid(3) + 
c     &         dsqrt(Bwid(3)*Bwid(3) + 4.*slopm*a(3)) )
c     &         /slopm
        theta(3)= rtnewt_theta(a(3), Dia, Amax) !for pipe
        d(3)= 0.5*(1-cos(0.5*theta(3)))*Dia
      end if
      
c     Twid(3)= Bwid(3) + 2.*slopm*d(3)
      Twid(3)= Dia*dsin(0.5*theta(3)) !for pipe
      u(3)= Q(3)/a(3)
      c(3)= dsqrt( g*a(3)/Twid(3) )
      !*
      
      !~ ul,cl & dl --> xl
c     Twid(1)= Bwid(1) + 2.*slopm*d(1)
      Twid(1)= Dia*dsin(0.5*theta(1)) !for pipe
      c(1)= dsqrt( g*a(1)/Twid(1) )
      
      u_l1= u(1) - dt*(  c(3)*u(1) - c(1)*u(3) )/dx0
      u_l2= 1. - dt*( u(1) - u(3) - c(1) + c(3) )/dx0
      u_l= u_l1/u_l2
      
      c_l1= c(1) + dt*u_l*( c(1) - c(3) )/dx0
      c_l2= 1. + dt*( c(1) - c(3) )/dx0
      c_l= c_l1/c_l2
      
      d_l= d(1) + dt*( u_l - c_l )*( d(1) - d(3) )/dx0
            
      x_l= xinflow - ( u_l - c_l )*dt
      x(3)= x_l
      call divide(1,3,3)
      !*
      
      
cc      write(101,'(4(e12.5,1x))')x_l
cc      write(101,*)
      
      
      !~
      call bottom(3,3)
      b_l= Bwid(3) 
      s0_l= s0(3)
      nMa_l= nMa(3)
      !*
      
      
      !~
c     a_l= AP(d_l,b_l,slopm,1)
      a_l= AP(d_l,Dia,1)         !for pipe
c     pwet_l= AP(d_l,b_l,slopm,2)
      pwet_l= AP(d_l,Dia,2)      !for pipe 
      rH_l= a_l/pwet_l
      sfx_l= u_l*u_l*nMa_l*nMa_l/rH_l**(4./3.)
      src_l= g*(s0_l - sfx_l)*dt
      
	!~for rectangular/trapezoidal channel
c     coef1= -slopm*g/c_l
c     coef2= slopm*(-u_l + g*d_l/c_l -src_l) - Bwid(1)*g/c_l
c     coef3= Bwid(1)*(-u_l + g*d_l/c_l -src_l)
c     coef4= Qinf
      !*
      
c     x1= 0.73
c     x2= 5.
c     dinf= rtnewt_d(x1,x2,coef1,coef2,coef3,coef4)
      !~for rectangular/trapezoidal channel

	!~for circular pipe
	if(a(1) > 0.) then
	dinf= d_l + c_l/g*(Qinf/a(1)-u_l-src_l)
      else
      pause 'a(1) <= 0.'
      endif
      
	if(dinf >= Dia) then
	  write(101,*) time, dinf  
	  if(dinf <= Dia * 1.01) dinf = Dia
      endif  

      end