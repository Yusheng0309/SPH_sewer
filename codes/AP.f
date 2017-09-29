      !~ 2012/05/14
	!~ 2017/01/18 for pipe 
      module AP_comput
        implicit none

	interface AP !computer area or hydraulic radius
	  module procedure AP_trape
	  module procedure AP_circular
	end interface

	contains
          
      function AP_trape(d,b,m,ind)
      !trapezoidal/rectangular
      implicit none
           
      real*8:: AP_trape,d,b,m
      integer:: ind
      
      if(ind == 1)then !area
       AP_trape= (m*d + b)*d 
      else if(ind == 2)then !hydraulic radius
        AP_trape= 2.*dsqrt(m*m + 1.)*d + b 
      end if
           
      end function AP_trape


      function AP_circular(d,diameter,ind)
      !circular
      implicit none
           
      real*8:: AP_circular,d,theta,diameter
      integer:: ind
      
	if( (1-2*d/diameter) > 1 .or. (1-2*d/diameter) <-1) then
c	write(*,*) 'd=',d,' dia=',diameter
	 if(d <= diameter*1.001) then
	  d = diameter
c	 write(*,*) 'd~diameter * 1.001'
	 else
	  pause 'acos: domain error'
	 endif
	endif
      theta=2*acos(1-2*d/diameter)
      
      if(ind == 1)then !area                        
	   AP_circular= diameter*diameter*(theta-dsin(theta))/8
      else if(ind == 2)then !hydraulic radius
	   AP_circular= theta*diameter/2
      end if
           
      end function AP_circular

	end module