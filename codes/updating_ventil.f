      subroutine check_ventil
	!2017/02/14
      !for determining ventilation or not
	!if ventilation =0, if not ventilation/pressurised =1
	!bc the state of 0 can be spreaded, 1 cannot. (e.g.0*1=0)
	!
	!free surface (surcharge=0, vent=0)  H_pressure = 0 
	!    d=d(area),theta=theta(area) 
	!    output depth = d+0,  M. eqn. depth = d+0
	!-------------------------------------------------
	!depression   (surcharge=0, vent=1)  H_pressure < 0
	!    d=Dia, theta=theta(area)=theta(d+H)  <-???????????????
	!    output depth = Dia+H_press, M. eqn. depth = Dia+H_press
	!--------------------------------------------------
	!pressurised  (surcharge=1, vent=1) H_pressure > 0
	!    d=Dia, theta=2*pi, area=Amax 
	!    output depth = Dia, M. eqn. depth = Dia+h_press
	!-------------------------------------------------- 

      !algorithm and rules
	!if surcharge =1 => vent =1 cannot be changed. (pressurised)
	!if surcharge =0 => vent =1/0 changeable (free surface/depression)
      !in/outflow both s and v are imposed and cannot be changed
	!         
	!       free surface            depression
	!           s=0         s=1         s=0      s=1
      !           v=0        (v=1)        v=1     (v=1)
      !-----------------------------------------------------------
      !inflow|000000000000|11111111111|111111111|111111111|outflow
	!-------~~~~~~~~~~~~-------------+++++++++------------------
      !**if a (s=0) region is surrounded(up and downstream) by (s=1, v=1)
	!   => it becomes (s=0, v=1) (depression) (+++++++)
	!**if there is a point of v=0 in the (s=0) region, or the region is connected to v=0 region
	!   => it becomes (s=0, v=0) (free surface) (~~~~~~~)
	!** the searching range is ended 
	!       while the up and downstream of region both encounter s=1 or in/outflow boundary
      !                                              e.g.        +++++ /  ~~~~~~
	!
      !ps in/ouflow boundary: can only be (s=1, v=1), (s=0, v=0)
      implicit none
      
      include 'common_global.dat'

	integer:: kk, sid(np), vent_temp(np)
	real*8:: sx(np)
   
      !vent state at in/outflow boundary
      do i= nvirp1,np
        if(typep(i) == -1) vent(i) = 0
	  if(typep(i) ==  1) vent(i) = 1
      enddo  


      !give index bc particle id is not in sequence with position
      nn= 0
      do i= nvirp1,np
        if(typep(i) >= -1)then
          nn= nn + 1
          sx(nn)= x(i)
          sid(nn)= i
        end if
      end do
      call bubble_sort(sx(1:nn),sid(1:nn),nn)

      !give a temp (ordered) index for computing ventilation state
      do kk=1,nn
	 if( typep(sid(kk)) /= 0 ) then !in/outflow, vent is imposed.
	  vent_temp(kk) = vent( sid(kk) ) 
       else
	  vent_temp(kk) = 1 !others are initialized by 1
	 !**elseif 
       !special set up - ventilation shaft
	 !at particular location = x0 
	 !vent( at x0 ) = 0
       !to be included
	 endif
	enddo
      !scanning from up to downstream
	do kk=1,nn !only notsurcharge particle can change vent state
	 if( (typep(sid(kk))==0) .and. (surcharge(sid(kk))==0) ) then
       vent_temp(kk)=vent_temp(kk-1) * vent_temp(kk+1)
       endif
	enddo
      !scanning from down to upstream
	do kk=nn,1,-1
       if( (typep(sid(kk))==0) .and. (surcharge(sid(kk))==0) ) then
       vent_temp(kk)=vent_temp(kk-1) * vent_temp(kk+1)
       endif
	enddo

	do kk=1,nn
	  vent(sid(kk))=vent_temp(kk)
	enddo


      end