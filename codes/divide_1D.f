      !~ 2012/01/30
      subroutine divide(kindp,istart,iend)
      
      implicit none
      
      include 'common_global.dat'
      
      
      open(30,file= 'divide.out')
      
      
      if(kindp == 1)then !fluid particles
        do i= istart,iend
          if(x(i) >= xinlet .and. x(i) <= xoutlet)then
            dx= x(i) - xmin
            
            !~
            icell= int(0.5*dx/hmax)+1 
            !*
            
            if(icell < 1)then
              write(30,'(6(i6,1x))')itime,kindp,i,icell
              write(30,'(1(e12.5,1x))')x(i)
              pause 'icell < 1 in sub divide'
            else if(icell > ncell)then
              write(30,'(6(i6,1x))')itime,kindp,i,icell
              write(30,'(3(e12.5,1x))')x(i),dx,hmax
              pause 'icell > ncell in sub divide'
            end if
            
            !~
            nc(kindp,icell)= nc(kindp,icell) + 1
            nn= nc(kindp,icell)
            !*
            
            if(nn > nplink)then
              write(30,'(6(i6,1x))')itime,kindp,i,nn
              pause 'nc > nplink in sub divide'
            end if
            
            !~
            ibox(kindp,icell,nn)= i
            iipar(i)= icell
            !*
          else
            !write(30,'(2(i7,1x),e12.5)')i,typep(i),x(i)
          end if
        end do
      else if(kindp == 2)then !virtual bed particles
        do i= istart,iend
          dx= xvirb(i) - xmin
          
          !~
          icell= int(0.5*dx/hmax) + 1 
          !*
           
          if(icell < 1)then
            write(30,'(4(i6,1x))')itime,kindp,i,icell
            write(30,'(1(e12.5,1x))')xvirb(i)
            pause 'icell < 1 in sub divide'
          else if(icell > ncell)then
            write(30,'(4(i6,1x))')itime,kindp,i,icell
            write(30,'(1(e12.5,1x))')xvirb(i)
            pause 'icell > ncell in sub divide'
          end if
          
          !~
          nc(kindp,icell)= nc(kindp,icell) + 1
          nn= nc(kindp,icell)
          !*
          
          if(nn > nplink)then
            write(30,'(4(i6,1x))')itime,kindp,i,nn
            pause 'nc > nplink in sub divide'
          end if
          
          !~
          ibox(kindp,icell,nn)= i
          !*
        end do
      end if
      
      
      end