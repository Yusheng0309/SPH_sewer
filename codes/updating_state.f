      subroutine state(kk)
      !2017/02/12 
	!for determining the state of pressurised flow or free surface flow
	!surcharge(i) = 1 for a pressurised flow 
	!surcharge(i) = 0 for a free surface flow
	!compute for the kk th particle

      implicit none
      
      include 'common_global.dat'

	integer:: kk

      !algorithm of updatnig the sate E
	!Bourdarias and Gerbi, 2007

      if(typep(kk)==0) then !* !only state of inner particle can change

c	surcharge_pre(kk) = surcharge(kk)

c	 if(surcharge_pre(kk) == 0) then
	   if( a(kk) < Amax ) then 
	   surcharge(kk) = 0
	   else
	   surcharge(kk) = 1
	   endif
      
c	 else !(surcharge_pre ==1)
c	   if( a(kk) >= Amax ) then
c	   surcharge(kk) = 1
c	   else          !near inflow/outflow -> same?
c	   surcharge(kk)= surcharge_pre(kk-1) * surcharge_pre(kk+1)
c         endif
      
c	 endif

      elseif(typep(kk)==1) then  !outflow boundary condition
       surcharge(kk) = 1   !downstream: pressurised
	                     !sur = 1 and vent = 1 cannot be changed.

      elseif(typep(kk)==-1) then  !inflow boundary condition
       surcharge(kk) = 0   !upstream: free surface
	                     !sur = 0 and vent = 0 cannot be changed.
		
	endif !*  
 

	end