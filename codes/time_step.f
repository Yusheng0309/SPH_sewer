      !~ 2012/07
      subroutine time_step
      
	use AP_comput
      implicit none
      
      include 'common_global.dat'
      
      real*8:: ci
      real*8:: dtp
      
      open(31,file= 'time_step.out')
      
      
      do i= nvirp1,np
        if(typep(i) >= -1)then
          Qi= Q(i)
          ai= a(i)
          di= d(i)
          
c         Twid(i)= Bwid(i) + 2.*slopm*d(i)
c         Twid(i)= Dia*sin(0.5*theta(i))
          
c         ci= dsqrt(g*ai/Twid(i)) !c is given by a when pressurised ![Sanders and Bradford, 2011]
          ci= dsqrt(g*di)         !c is tunable [Bourdarias and Gerbi, 2007, sec. 4]

          ui= Qi/ai
          dtp= courant*dx0/(ci + abs(ui))
          if(dtp < dt) dt= dtp
        end if
      end do 
          
      write
     & (31,'("time = ",e10.3," itime = ",i6,", time step = "
     & ,e10.3)')
     & time,itime,dt
      
      
      end