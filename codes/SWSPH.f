      !~ 2012/09
      program SWSPH
      
	use AP_comput      
      implicit none     
      
      include 'common_global.dat'
      
      real*8:: s1,s2
      character:: supp*8,name*100,inifilname*300
  
      
      write(*,*)'Keyin the filename of basic'
      read(*,*)inifilname
      open(11,file= inifilname)
      write(*,*)
      write(*,*)'Keyin the filename of set_info'
      read(*,*)inifilname
      open(12,file= inifilname)
      write(*,*)
      open(21,file= 'sph.out')
      
      
      call cpu_time(s1)
      
      
      !*******************************************
      read(11,*)
      read(11,*)dx0
      read(11,*)
      read(11,*)tout
      read(11,*)
      read(11,*)courant
      read(11,*)
      read(11,*)dt_ini
      read(11,*)
      read(11,*)slopm
      read(11,*)
      read(11,*)udcon
      dt= dt_ini
  
      hmax= 8.*dx0
      read(12,*)
      read(12,*)xinlet,xinflow,xoutflow,xoutlet
      xmax= xoutlet
      xmin= xinlet
      nct= int(0.5*(xmax - xmin)/hmax) + 1
      
      read(12,*)
      read(12,*)Qinf,dinf,Qoutf,doutf
      !*******************************************

      !geometry of pipe 
	!Diameter and Areamax are initialized in common_global.dat
	
	!~ input data
      call input_geo
      !*
            
      write(*,*)'Keyin the final simulated time'
      read(*,*)time_max
      write(*,*)
      
      
      time= 0.
      itime= 0
      ngrab= 0
      nstr= 0
                    	
      do while(time <= time_max)
        itime= itime + 1
        time= time + dt
        
        !~ output
        if(itime == 1)then
          write(*,'(i7,1x,a7,i7,a8,e12.5)')
     &      ngrab,"step = ",itime," time = ",time
          write(*,*)
          
          write(supp,'(i4.4)')ngrab+1
          name= 'output/'//supp
          
          open(23,file= name)
          call output(23)
          close(23)
        end if
        !*
        

        !~
        if((time-ngrab*tout) >= tout)then
          iout= 1
          ngrab= ngrab + 1
        else
          iout= 0
        end if
        !*
        
        !~
        call time_integration
        !*
        
        
        !~ in/out-flow algorithm
        nadd= 0
        do i= nvirp1,np
          if(typep(i) == 0)then
            if(x(i) >= xoutflow .and. x(i) < xoutlet)
     &        call bc_outflow(i)     !detect inner particles if over the outflow boundary
          else if(typep(i) == -1)then
            if((x(i)-xinflow) >= 1.e-10)then     !detect  inflow particles if over the inflow boundary
              nadd= nadd + 1
              typep(i)= 0
              xinf(nadd)= x(i)
            end if
          else if(typep(i) == 1)then
            if(x(i) >= xoutlet) call bc_outlet(i)   !detect outflow particles if over the outlet boundary
          end if
        end do
        
        if(nadd > 0) call bc_inflow
        !*
        
        
        !~ output
        if(iout == 1)then
          write(*,'(i7,1x,a7,i7,a8,e12.5)')
     &      ngrab,"step = ",itime," time = ",time
          write(*,*)
          
          write(supp,'(i4.4)')ngrab+1
          name= 'output/'//supp
          
          open(23,file= name)
          call output(23)
          close(23)
        end if
        !*
      end do
      !*
      
            
      write(*,*)
      call cpu_time(s2)      
      write(*,*)'Elapsed CPU time = ', s2-s1
      !write(21,'("Elapsed CPU time = ",e13.6)')s2-s1
      
      pause 'well done'
      
      
      end