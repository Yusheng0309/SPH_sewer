      !~ 2012/02/07
      subroutine bc_outflow(ip)
      
      implicit none
      
      include 'common_global.dat'
      
      integer:: ip
      
      
      typep(ip)= 1
      Q(ip)= Qoutf
      call state(ip)
      
      end