      !~ 2012/05/14
      function rtnewt_d(x1,x2,coef1,coef2,coef3,coef4)
      
      implicit none
      
      integer, parameter:: jmax= 200
      integer:: j
      real*8:: rtnewt_d,x1,x2,xacc
      real*8:: coef1,coef2,coef3,coef4
      real*8:: df,dx,f
      
      
      rtnewt_d= 0.5*(x1+x2)
      xacc= 1.e-10
      do j= 1,jmax
        call funcd(rtnewt_d,f,df,coef1,coef2,coef3,coef4)
        
        dx= f/df
        rtnewt_d= rtnewt_d - dx
        
        if((x1-rtnewt_d)*(rtnewt_d-x2) < 0.)then
		  write(*,*)rtnewt_d
            pause 'rtnewt_d jumped out of brackets'
	  end if

        if(abs(dx) < xacc) return
      end do
      
      pause 'rtnewt_d exceeded maximum iterations'
      
      end
      
      
      
      subroutine funcd(x,fn,df,coef1,coef2,coef3,coef4)
      
      implicit none
      
      real*8:: coef1,coef2,coef3,coef4,df,fn,x
      
      
      fn= coef1*(x**3.) + 
     &    coef2*(x**2.) +
     &    coef3*(x**1.) +
     &    coef4
      df= 3.*coef1*(x**2.) + 
     &    2.*coef2*(x**1.) +
     &    1.*coef3
      
      return
      
      end