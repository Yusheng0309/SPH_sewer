      !~ 2012/09
      subroutine time_integration
      
	use AP_comput
      implicit none     
      
      include 'common_global.dat'
      
      real*8:: dQ,dx2
      
      open(27,file= 'intergration.out')
                     
      
      dx2= 1.e-6
      
      
      !~ dQdt at step n
      !*******************************************
      if(itime == 1)then
        call ini_divide(1)
        call divide(1,1,np)
        call ini_divide(2)
        call divide(2,1,nvirb)
        
        !~ Bwid, zb, & nMa
        call bottom(nvirp1,np)
        call bottom(1,2)
        do i= nvirp1,np
          if(typep(i) == -1)then
            Bwid(i)= Bwid(1)
            s0(i)= s0(1)
            nMa(i)= nMa(1)
          else if(typep(i) == 1)then
            Bwid(i)= Bwid(2)
            s0(i)= s0(2)
            nMa(i)= nMa(2)
          end if
        end do
        !*
        
        call ac_main
      else if(itime > 1)then
        call ini_divide(1)
        call divide(1,1,np) 
        
        !~ Bwid, zb, & nMa
        call bottom(nvirp1,np)
        call bottom(1,2)
        do i= nvirp1,np
          if(typep(i) == -1)then
            Bwid(i)= Bwid(1)
            s0(i)= s0(1)
            nMa(i)= nMa(1)
          else if(typep(i) == 1)then
            Bwid(i)= Bwid(2)
            s0(i)= s0(2)
            nMa(i)= nMa(2)
          end if
        end do
        !*
        
        !~ area
        if(udcon == 1)then
          call SWarea_bb
        else if(udcon == 2)then
          call SWarea_ss
        else if(udcon == 3)then
          call SWarea_sb
        end if        
        !*
        
        call ac_main
      end if
      !*******************************************
      
      
            
      !~ Q at step n + 1/2
      !*****************************************
      if(itime == 1)then
        do i= nvirp1,np
          if(typep(i) == 0)then   !inner particles
            Q(i)= Q(i) + 0.5*dt*dQdt(i)
            Q_min(i)= Q(i)
          end if                         
        end do                           
      else if(itime > 1)then
        do i= nvirp1,np
          if(typep(i) == 0)then   !inner particles
            Q(i)= Q_min(i) + dt*dQdt(i)
            Q_min(i)= Q(i)
          end if
        end do
      end if
      
      !~ Qoutf
      if(udcon == 1 .or. udcon == 3)then
        call Q_sub
      else if(udcon == 2)then
        call Q_super
      end if
      !*
      
      do i= nvirp1,np
        if(typep(i) == -1)then   !inflow particles
          Q(i)= Qinf
          Q_min(i)= Q(i)
        else if(typep(i) == 1)then   !outflow particles
          Q(i)= Qoutf
          Q_min(i)= Q(i)
        end if
      end do
      !*****************************************
      
      
      
      !~ x at step n+1
      !*******************************************
      do i= nvirp1,np
        if(typep(i) >= -1)then
          x(i)= x(i) + dt*Q(i)/a(i)
        end if
      end do
      !*******************************************
      
      
      
      !~ Q at step n+1
      !*******************************************
      do i= nvirp1,np
        if(typep(i) == 0)then   !inner particles
          Q(i)= Q(i) + 0.5*dt*dQdt(i)
        end if
      end do
      

      !~ Qoutf
      if(udcon == 1 .or. udcon == 3)then
        call Q_sub
      else if(udcon == 2)then
        call Q_super
      end if
      !*
      
      do i= nvirp1,np
        if(typep(i) == -1)then   !inflow particles
          Q(i)= Qinf
        else if(typep(i) == 1)then   !outflow particles
          Q(i)= Qoutf
        end if
      end do
      !*******************************************
      
      
      
      !~ variable time step
      dt= dt_ini
      call time_step 
      !*
          
      end