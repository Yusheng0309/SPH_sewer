      !~ 2012/01/30
      subroutine ini_divide(kindp)
      
      implicit none
      
      include 'common_global.dat'
      
      
      do i= 1,ncell
        nc(kindp,i)= 0
        ibox(kindp,i,1:nplink)= 0
      end do
      
      if(kindp == 1)then
        do i= 1,nmax
          iipar(i)= 0
        end do
      end if
  
      
      end