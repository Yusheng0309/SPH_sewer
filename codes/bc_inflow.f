      !~ 2012/03/22
      subroutine bc_inflow
      
	use AP_comput
      implicit none
      
      include 'common_global.dat'
      real*8:: rtnewt_theta

      nn= np
      nn1= 0 
      loop1: do j= 1,nadd
        if(nstr > 0)then
          !~
          if(nstr > 500)then
            write(*,*)'- - - - - - - - - -'
            write(*,'("nstr")')
            write(*,'(1(i3,1x))')nstr
            write(*,*)'- - - - - - - - - -'
            pause 'nstr > 500 in sub bc_inflow'
          end if
          !*
          
          nstr= nstr - 1
          nn1= nn1 + 1
          loop2: do i= 1,500
            if(astr(i) == 1)then
              astr(i)= 0
              nn2= istr(i)
              typep(nn2)= -1
              x(nn2)= xinf(j) - xinflow + xinlet
              Q(nn2)= Qinf
c             a(nn2)= AP(dinf,Bwid(1),slopm,1)
              a(nn2)= AP(dinf,Dia,1)
              theta(nn2)= rtnewt_theta(a(nn2), Dia, Amax)
              aini(nn2)= a(nn2)
              mass(nn2)= dx0*a(nn2)
              h(nn2)= 1.3*dx0
              hini(nn2)= h(nn2)
	        call state(nn2) !updating surcharge state
              cycle loop1
            end if
          end do loop2
        else
          !~
          if(nstr /= 0)then
            write(*,*)'- - - - - - - - - -'
            write(*,'("nstr")')
            write(*,'((i3,1x))')nstr
            write(*,*)'- - - - - - - - - -'
            pause 'nstr /= 0  in sub bc_inflow'
          end if
          !*
          np= np + 1
          nn2= np
          typep(nn2)= -1
          x(nn2)= xinf(j) - xinflow + xinlet
          Q(nn2)= Qinf
c         a(nn2)= AP(dinf,Bwid(1),slopm,1)
          a(nn2)= AP(dinf,Dia,1)
          theta(nn2)= rtnewt_theta(a(nn2), Dia, Amax)
          aini(nn2)= a(nn2)      
          mass(nn2)= dx0*a(nn2)
          h(nn2)= 1.3*dx0
          hini(nn2)= h(nn2)
	    call state(nn2) !updating surcharge state
        end if
      end do loop1
      
      
      !~
      if(np /= nn + nadd - nn1)then
        write(*,*)'- - - - - - - - - -'
        write(*,'("nstr")')
        write(*,'(2(i3,1x))')np,nn + nadd - nn1 
        write(*,*)'- - - - - - - - - -'
        pause 'np error in sub bc_inflow'
      end if
      !*
      
      
      end