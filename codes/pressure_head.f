      !~2017/02/13
      real*8 function pressure_head(kk) 
	!compute pressure head
	implicit none
	
	include 'common_global.dat'

	integer:: kk
	real*8:: aws !celerity, acoustic wave-speed


      !*****free pressure condition-> formulation 
	! can computate pressure, despression, free surface
      aws=5 !7.5 !5 !1 !10 !imposed case by case

c	if( (surcharge(kk)==0) .and. (vent(kk)==0) )then
c         pressure_head = 0.0
c	else
c         pressure_head = (aws*aws/g) * (a(kk)-Amax) / Amax
c	endif
      
	!******forced pressure condition -> imposed constant pressure
      if( (surcharge(kk)==0) .and. (vent(kk)==0) )then
         pressure_head = 0.0
	else 
         pressure_head = 0.334
	endif


      return
	end