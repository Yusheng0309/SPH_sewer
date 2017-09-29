      !~ 2011/07/05
      subroutine kernel(r,dx,sh,w,dwdr,dwdx,ia,ib)   
        
      implicit none
      
      include 'common_global.dat'
      
      real*8:: sh,s,factor
      integer:: ia,ib
      
      
      s= r/sh
      factor= 1./sh
                                              
      if(s >= 0 .and. s <= 1.)then          
        w= factor*(2./3.-s*s+s**3/2.)
        dwdr= factor*(-2*s+1.5*s*s)/sh
        dwdx= factor*dx*(-2.+3./2.*s)/sh**2  
      else if(s > 1. .and. s <= 2.)then        
        w= factor*(2.-s)**3/6.
        dwdr= factor*(-0.5*(2-s)**2)/sh
        dwdx= -factor*dx*(2.-s)**2/sh/r/2.
      else
        w= 0.
        dwdr= 0.
        dwdx= 0.
      end if
          
      
      end