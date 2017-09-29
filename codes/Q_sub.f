      !~ 2012/07/23
      subroutine Q_sub   ![paper, 2005 2010 2011]
      
	use AP_comput      
      implicit none
      
      include 'common_global.dat'
      
      real*8:: u_r,c_r,d_r,x_r
      real*8:: u_r1,u_r2,c_r1,c_r2
      real*8:: b_r,s0_r,nMa_r,a_r,pwet_r,rH_r,sfx_r,src_r
      real*8:: dx2, rtnewt_theta
      
      open(24,file= 'depth.out')
                 
      
      nn= 0
      dx2= 1.e-6
      
      
      !~ The method of specified time interval - Qoutf
      !*******************************************
      !~ ghost grid 
      !***********************
      Qi= 0.
      ai= 0.
      sumw_out= dx2
      
      x(4)= xoutflow - dx0
      call divide(1,4,4)
      
      ii= iipar(4) !self
      if(nc(1,ii) > 0)then
        do k= 1,nc(1,ii)
          j= ibox(1,ii,k)
          if(typep(j) == 0)then
            dx= x(4) - x(j)
            r= abs(dx)
            hi= h(4)
            
            if(r <= 2.*hi)then
              call kernel(r,dx,hi,w,dwdr,dwdx,4,j)
              
              mj= mass(j)
              
              Qj= Q(j)
              aj= a(j)
              
              Qi= Qi + mj*Qj*w/aj
              ai= ai + mj*w
              sumw_out=  sumw_out + mj*w/aj
            end if
          end if
        end do
        
        iil= ii - 1 !neighbor-left
        if(iil > 0 .and. nc(1,iil) > 0)then
          do k= 1,nc(1,iil)
            j= ibox(1,iil,k)
            if(typep(j) == 0)then
              dx= x(4) - x(j)
              r= abs(dx)
              hi= h(4)
              
              if(r <= 2.*hi)then
                call kernel(r,dx,hi,w,dwdr,dwdx,4,j)
                
                mj= mass(j)
                Qj= Q(j)
                aj= a(j)
                
                Qi= Qi + mj*Qj*w/aj
                ai= ai + mj*w
                sumw_out=  sumw_out + mj*w/aj
              end if
            end if
          end do
        end if
        
        iir= ii + 1 !neighbor-right
        if(iir > 0 .and. nc(1,iir) > 0)then
          do k= 1,nc(1,iir)
            j= ibox(1,iir,k)
            if(typep(j) == 0)then
              dx= x(4) - x(j)
              r= abs(dx)
              hi= h(4)
              
              if(r <= 2.*hi)then
                call kernel(r,dx,hi,w,dwdr,dwdx,4,j)
                
                mj= mass(j)
                Qj= Q(j)
                aj= a(j)
                
                Qi= Qi + mj*Qj*w/aj
                ai= ai + mj*w
                sumw_out=  sumw_out + mj*w/aj
              end if
            end if
          end do
        end if
      end if
      
      Q(4)= Qi/sumw_out
      a(4)= ai/sumw_out
      
      call bottom(4,4)
      
      if(slopm <= 1.e-6)then
        d(4)= a(4)/Bwid(4) 
      else
c       d(4)= 0.5*( -Bwid(4) + 
c    &         dsqrt(Bwid(4)*Bwid(4) + 4.*slopm*a(4)) )
c    &         /slopm
        theta(4)= rtnewt_theta(a(4), Dia, Amax) !for pipe
        d(4)= 0.5*(1-cos(0.5*theta(4)))*Dia
      end if
      
c     Twid(4)= Bwid(4) + 2.*slopm*d(4)
      Twid(4)= Dia*dsin(0.5*theta(4)) !for pipe
      
      u(4)= Q(4)/a(4)
      c(4)= dsqrt(g*a(4)/Twid(4))
      !***********************
                 
      
      
      !~ ur,cr & dr --> xr
c     Twid(2)= Bwid(2) + 2.*slopm*d(2)
      Twid(2)= Dia*dsin(0.5*theta(2)) !for pipe
      c(2)= dsqrt(g*a(2)/Twid(2))
      
      u_r1= u(2) + dt*( -u(2)*c(4) + c(2)*u(4) )/dx0
      u_r2= 1. + dt*( u(2) - u(4) + c(2) - c(4) )/dx0
      u_r= u_r1/u_r2
      
      c_r1= c(2) - dt*u_r*( c(2) - c(4) )/dx0
      c_r2= 1. + dt*( c(2) - c(4) )/dx0
      c_r= c_r1/c_r2
      
      d_r= d(2) - dt*( u_r + c_r )*( d(2) - d(4) )/dx0
      
      x_r= xoutflow - ( u_r + c_r )*dt
      x(4)= x_r
      call divide(1,4,4)
      !*
      
      !~
      call bottom(4,4)
      b_r= Bwid(4)
      s0_r= s0(4)
      nMa_r= nMa(4)
      !*
      
      !~
c     a_r= AP(d_r,b_r,slopm,1)
      a_r= AP(d_r,Dia,1)
c     pwet_r= AP(d_r,b_r,slopm,2)
      pwet_r= AP(d_r,Dia,2)
      rH_r= a_r/pwet_r
      sfx_r= u_r*u_r*nMa_r*nMa_r/rH_r**(4./3.)
      src_r= g*(s0_r - sfx_r)*dt
      
      u(2)= u_r - g*(doutf - d_r)/c_r + src_r

c     a(2)= AP(doutf,Bwid(2),slopm,1)
      a(2)= AP(doutf,Dia,1)
      Qoutf= a(2)*u(2)
      !*
      !*******************************************      
                  
            
      end