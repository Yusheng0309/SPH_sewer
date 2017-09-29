      !~ 2017/01/18
      function rtnewt_theta(area,Diameter,Areamax)
      !solve the wetted angle of the pipe by giving A  
      implicit none

      real, parameter:: pi= 3.14159 !4*atan(1.0))
      integer, parameter:: jmax= 200
      integer:: j
      real*8:: rtnewt_theta,x1,x2,xacc
      real*8:: area, Diameter, Areamax
      real*8:: df,dx,f

c      Dia = 
c      Amax = pi*Dia*Dia/4

      if(area >= Areamax) then
	 rtnewt_theta = 2*pi
	 return

      else  !Newton iteration
	!0 ~ 2pi
	x1=0
	x2=2*pi
         
      rtnewt_theta= 0.5*(x1+x2)
      xacc= 1.e-10
      do j= 1,jmax
	   call func_theta(rtnewt_theta,f,df,area,Areamax)
 
         if(df <= 0. ) then
	     dx = 0.  
	     rtnewt_theta = 2*pi
         
	    else
	     dx= f/df
           rtnewt_theta= rtnewt_theta - dx
        
           if((x1-rtnewt_theta)*(rtnewt_theta-x2) < 0.)then
		    write(*,*)rtnewt_theta
              pause 'rtnewt_theta jumped out of brackets'
	     end if
          endif 
	  
         if(abs(dx) < xacc) return
      end do
      
      pause 'rtnewt_theta exceeded maximum iterations'
      
	endif
      end
      
      
      
      subroutine func_theta(xw,fn,df,area,Areamax)
      
      implicit none
      include 'common_global.dat'
      real*8:: Areamax !Diameter
      real*8:: area,fn,df,xw,Ar
c      real, parameter:: pi= 3.14159 !4*atan(1.0))

      Ar = area / Areamax 
 
      !given A, solve x for A = D^2*(x-sinx)/8 
	!                  A/Amax =( D^2*(x-sinx)/8 )/Amax
	!                 2*pi*Ar = x-sinx

      fn= xw-dsin(xw)-2*pi*Ar
      
	df= 1-cos(xw)
      
      return
      
      end