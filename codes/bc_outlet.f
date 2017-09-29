      !~ 2012/01/30
      subroutine bc_outlet(ip)
      
      implicit none
      
      include 'common_global.dat'
      
      integer:: ip
      

      typep(ip)= -5
      nstr= nstr + 1
      
      !~
      if(nstr > 500)then
        write(*,*)'- - - - - - - - - -'
        write(*,'("nstr")')
        write(*,'((i3,1x))')nstr
        write(*,*)'- - - - - - - - - -'
        pause 'nstr > 500 in sub bc_outlet'
      end if
      !*
      
      do i= 1,500
        if(astr(i) == 0)then
          astr(i)= 1
          istr(i)= ip
          x(ip)= xoutlet + 20.*dx0
          exit
        end if
      end do
      
      
      end