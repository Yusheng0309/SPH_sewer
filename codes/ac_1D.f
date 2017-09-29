      !~ 2012/09
      subroutine ac_main
      
	use AP_comput
      implicit none
      
      include 'common_global.dat'
      
      real*8:: pweti,rHi,sfxi
      
      open(32,file= 'ac_main.out')

      
      do i= nvirp1,np
        dQdt(i)= 0.
      end do
      
      nn= nct
      do i= 1,nn
        nn1= nc(1,i) 
        if(nn1 > 0)then
          iir= i+1
          if(iir <= nn .and. nn1 > 0)then
            kind_cell= 1
            call single_step(i,iir,kind_cell) !neighbor
          end if
          
          kind_cell= 2
          call single_step(i,i,kind_cell) !self
        end if    
      end do  
      
      !~ sf
      do i= nvirp1,np
        if(typep(i) == 0)then
          Qi= Q(i)
          ai= a(i)
c          pweti= AP(d(i),Bwid(i),slopm,2)
          if( vent(i) == 1 ) then !pipe
	    rHi= 0.25*Dia !pressurised/depression
	    else
          pweti= AP(d(i),Dia,2) !free surface
          rHi= ai/pweti
	    endif
          sfxi= nMa(i)**2*Qi*abs(Qi)/ai/ai/rHi**(4./3.)
          dQdt(i)= dQdt(i) + g*ai*( s0(i) - sfxi )
        end if
      end do
      !*
      
      
      end