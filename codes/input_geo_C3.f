      !~ 2012/09
      subroutine input_geo
      
	use AP_comput
      implicit none     
      
      include 'common_global.dat'
      
      real*8:: xx,zz,dzz,ak1,ak2,ak3,dxb,zb
      integer:: nf
      character*100:: filname_zb
          
      write(*,*)'Keyin the input filename of bed elevation'
      read(*,*)filname_zb
      write(*,*)
      
      open(13,file= filname_zb)
      open(22,file= 'virbed.out')
      
          
      write(22,*)'variables= "x"'
      write(22,*)'"d"'
      write(22,*)'"b"'
      write(22,*)'"zbso"'
              
      ak1= -0.111051
      ak2=  0.026876
      ak3= -0.217567
      
      dxb= dx0
      
      !~ virtual bed particles
      !*******************************************
      read(13,*)nvirb

      do i= 1,nvirb
        read(13,*)xvirb(i),zb,s0virb(i),zz
        
cc      nn= 0
cc      do i= -4,404
cc        nn= nn + 1
cc        xvirb(nn)= dxb*i
cc        xx= xvirb(nn)
cc        
cc        if(xx < 0.)then
cc        	xx= 0.
cc        	zz=  0.723449*( 1. - dtanh(0.001*xx-0.3) )
cc          dzz= -0.001*0.723449*( 1. - dtanh(0.001*xx-0.3)**2 )
cc        else if(xx >= 0. .and. xx <= 300.)then
cc          zz=  0.723449*( 1. - dtanh(0.001*xx-0.3) )
cc          dzz= -0.001*0.723449*( 1. - dtanh(0.001*xx-0.3)**2 )
cc        else if(xx > 300. .and. xx <= 600.)then
cc          zz=  0.723449*( 1. - dtanh(6.*(0.001*xx-0.3))/6. )
cc          dzz= -0.001*0.723449*( 1. - dtanh(6.*(0.001*xx-0.3))**2  )
cc        else if(xx > 600. .and. xx <= 1000.)then
cc          zz=  0.75 + 0.6*dexp(0.001*xx-1.) +
cc     &         ak1*dexp(-20.*(0.001*xx-0.6)) +
cc     &         ak2*dexp(-40.*(0.001*xx-0.6)) +
cc     &         ak3*dexp(-60.*(0.001*xx-0.6))
cc          dzz= 0.001*0.6*dexp(0.001*xx-1.)        -
cc     &         0.020*ak1*dexp(-20.*(0.001*xx-0.6)) -
cc     &         0.040*ak2*dexp(-40.*(0.001*xx-0.6)) -
cc     &         0.060*ak3*dexp(-60.*(0.001*xx-0.6))         	            
cc        else i f(xx > 1000.)then
cc          xx= 1000.
cc          zz=  0.75 + 0.6*dexp(0.001*xx-1.) +
cc     &         ak1*dexp(-20.*(0.001*xx-0.6)) +
cc     &         ak2*dexp(-40.*(0.001*xx-0.6)) +
cc     &         ak3*dexp(-60.*(0.001*xx-0.6))
cc          dzz= 0.001*0.6*dexp(0.001*xx-1.)        -
cc     &         0.020*ak1*dexp(-20.*(0.001*xx-0.6)) -
cc     &         0.040*ak2*dexp(-40.*(0.001*xx-0.6)) -
cc     &         0.060*ak3*dexp(-60.*(0.001*xx-0.6))  
cc        end if
cc                
cc        zbsovirb(nn)=
cc     &    ( 1. - 400.*(10.+2.*zz)/g/(10.+zz)**3/zz**3 )*dzz +
cc     &     0.16 * ( 10. + 2.*sqrt(2.)*zz )**(4./3.) /
cc     &    (10.+zz)**(10./3.) / zz**(10./3.)
          
        
        Bwidvirb(i)= 10. !-> Dia = 20 
        massvirb(i)= dxb
        nMavirb(i)= 0 !0.02
        hvirb(i)= 1.3*dxb
        
        write(22,'(4(e12.5,1x))')xvirb(i),zz,Bwidvirb(i),s0virb(i)
      end do 
cc      nvirb= nn
      !*******************************************
      
      
      
      !~ fluid particles
      !*******************************************
      !typep = -5, storage particles 
      !typep = -3, virtual particles
      !typep = -1, inflow particles                                
      !typep =  0, inner particles                                
      !typep =  1, outflow particles                                 
            
      !~ virtual particles
      nn= 0
      do i= 1,5
        nn= nn+1
        if(i == 1) x(nn)= xinflow
        if(i == 2) x(nn)= xoutflow
        if(i == 3) x(nn)= 0.
        if(i == 4) x(nn)= 0.
        if(i == 5) x(nn)= 0.
                
        typep(nn)= -3
        
        Q(nn)= Qinf
        d(nn)= doutf
	  h_pressure(nn)= 0.
c       a(nn)= AP(d(nn),dble(10.),slopm,1)
        theta(nn)= 2*acos(1-2*d(nn)/Dia)
        a(nn)= Dia*Dia*(theta(nn) - sin( theta(nn) ) )/8
        aini(nn)= a(nn)
	  if( (i==1) .or. (i==3) )then !inflow free surface
	  surcharge(nn)= 0
	  surcharge_pre(nn)= 0
	  vent(nn)= 0
	  else  !outflow pressurised
	  surcharge(nn)= 1
	  surcharge_pre(nn)= 1
	  vent(nn)= 1
	  endif
        u(nn)= Q(nn)/a(nn)
        mass(nn)= dx0*a(nn)
        h(nn)= 1.3*dx0
        hini(nn)= h(nn)
      end do
      nvirp= nn
      nvirp1= nn + 1
      !*
      
      !~ fluid particles 
      do i= 1,nvirb-1+8
        nn= nn + 1
        x(nn)= dx0*(i-5+0.5)
        xx= x(nn)
        
        if(xx <= xinflow)then
          typep(nn)= -1
          
          Q(nn)= Qinf
          d(nn)= doutf
	    h_pressure(nn)= 0. 
c         a(nn)= AP(d(nn),dble(10.),slopm,1)
          theta(nn)= 2*acos(1-2*d(nn)/Dia)
          a(nn)= Dia*Dia*(theta(nn) - sin( theta(nn) ) )/8
          aini(nn)= a(nn)
	    surcharge(nn)= 0 
	    surcharge_pre(nn)= 0
	    vent(nn)= 0   !remember to set the vent state in updating_vent.f
          mass(nn)= dx0*a(nn)
          h(nn)= 1.3*dx0
          hini(nn)= h(nn)
          cycle
        else if(xx > xinflow .and. xx < xoutflow)then
        	typep(nn)= 0
        	
          Q(nn)= Qinf
          d(nn)= doutf
	    h_pressure(nn)= 0.
c         a(nn)= AP(d(nn),dble(10.),slopm,1)
          theta(nn)= 2*acos(1-2*d(nn)/Dia)
          a(nn)= Dia*Dia*(theta(nn) - sin( theta(nn) ) )/8
          aini(nn)= a(nn)
	    surcharge(nn)= 0
	    surcharge_pre(nn)= 0
	    vent(nn)= 0
          mass(nn)= dx0*a(nn)
          h(nn)= 1.3*dx0
          hini(nn)= h(nn)
          cycle
        else if(xx >= xoutflow)then
        	typep(nn)= 1
          
          Q(nn)= Qinf
          d(nn)= doutf
	    h_pressure(nn)= 0.
c         a(nn)= AP(d(nn),dble(10.),slopm,1)
          theta(nn)= 2*acos(1-2*d(nn)/Dia)
          a(nn)= Dia*Dia*(theta(nn) - sin( theta(nn) ) )/8
          aini(nn)= a(nn)
	    surcharge(nn)= 1
	    surcharge_pre(nn)= 1
	    vent(nn)= 1   !remember to set the vent state in updating_vent.f
          mass(nn)= dx0*a(nn)
          h(nn)= 1.3*dx0
          hini(nn)= h(nn)
          cycle
        end if
      end do
      np= nn
      !*
      !*******************************************
      
      
      end