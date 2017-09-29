      !~ 2012/07
      subroutine output(nf)
      
	use AP_comput
      implicit none
       
      include 'common_global.dat'
      
      integer:: nf,sid(np)
      real*8:: fr,dx2,sx(np),zb
      
	open(51,file= 'par_evol.out')
            
      write(nf,*)'variables= "x"'
      write(nf,*)'"d"'
	write(nf,*)'"hs"'
      write(nf,*)'"Q"'
      write(nf,*)'"A"'
      write(nf,*)'"u"'
      write(nf,*)'"Fr"'
      write(nf,*)'"theta"' !B
      write(nf,*)'"nMa"'
      write(nf,*)'"zb"'
      write(nf,*)'"surcharge"'
	write(nf,*)'"ventilation"'
      write(nf,*)'"id"'
      write(nf,*)'"type"'
      
      nn= 0
      do i= nvirp1,np
        if(typep(i) >= -1)then
          nn= nn + 1
          sx(nn)= x(i)
          sid(nn)= i
        end if
      end do
      call bubble_sort(sx(1:nn),sid(1:nn),nn)
      
      
      write(51,'(e12.5,1x,i5)')time,nn
      
                 
      do j= 1,nn   
        i= sid(j)
        Qi= Q(i)
        ai= a(i)
        di= d(i)
        ui= Qi/ai
        
c       fr= ui/dsqrt(g*ai/Twid(i))
        if(di<=0.) then
	   fr= 0. 
	  else
         fr= ui/dsqrt(g*di)
	  endif
        if(x(i)<=2) then !for case T3
	   zb=di+0.694594
        elseif( (x(i)>2) .and. (x(i)<6) ) then
	   zb=di+(-0.173648*x(i)+1.04189)
	  else
	   zb=di+0.0
	  endif
        
        write(nf,1001)x(i),di,h_pressure(i),Qi,ai,ui,
c    &       fr,Bwid(i),nMa(i),s0(i),i,typep(i)
     &       fr,theta(i),nMa(i),zb,surcharge(i),vent(i),i,typep(i)
      end do
      
          
1001  format(10(e13.6,1x),4(i6,1x))
      
	
      end