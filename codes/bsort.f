      !~ 2012/03/22
      subroutine bubble_sort(a,b,n)
      
      implicit none
      
      integer:: n,b(n),i,j,tempb
      real*8:: a(n),tempa
      
      
      do i= n-1,1,-1   
        do j= 1,i
          if(a(j) > a(j+1))then
            tempa= a(j)
            tempb= b(j)
            a(j)= a(j+1)
            b(j)= b(j+1)       
            a(j+1)= tempa
            b(j+1)= tempb
          end if
        end do
      end do
      
      return
      
      end