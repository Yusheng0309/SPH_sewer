      !~ 2012/07/23
      subroutine SWarea_bb   ![paper, 2005 2010 2011]
            
      use AP_comput
	implicit none
      
      include 'common_global.dat'
      
      integer:: index_conv   
      real*8:: Aprev(nmax),sumA(nmax),apha(nmax),resA(nmax)
      real*8:: dx2, rtnewt_theta, pressure_head
      
      open(24,file= 'depth.out')
                 
      
      nn= 0
      dx2= 1.e-6
      
      !~ Characteristic method - dinf
      Qi= 0.
      ai= 0.
      sumw_in= dx2
      
      x(3)= xinflow + dx0
      call divide(1,3,3)
      
      ii= iipar(3) !self
      if(nc(1,ii) > 0)then
        do k= 1,nc(1,ii)
          j= ibox(1,ii,k)
          if(typep(j) == 0)then
            dx= x(3) - x(j)
            r= abs(dx)
            hi= h(3)
            
            if(r <= 2.*hi)then
              call kernel(r,dx,hi,w,dwdr,dwdx,3,j)
              
              mj= mass(j)
              
              Qj= Q(j)
              aj= a(j)
              
              Qi= Qi + mj*Qj*w/aj
              ai= ai + mj*w
              sumw_in=  sumw_in + mj*w/aj
            end if
          end if
        end do
        
        iil= ii - 1 !neighbor-left
        if(iil > 0 .and. nc(1,iil) > 0)then
          do k= 1,nc(1,iil)
            j= ibox(1,iil,k)
            if(typep(j) == 0)then
              dx= x(3) - x(j)
              r= abs(dx)
              hi= h(3)
              
              if(r <= 2.*hi)then
                call kernel(r,dx,hi,w,dwdr,dwdx,3,j)
                
                mj= mass(j)
                Qj= Q(j)
                aj= a(j)
                
                Qi= Qi + mj*Qj*w/aj
                ai= ai + mj*w
                sumw_in=  sumw_in + mj*w/aj
              end if
            end if
          end do
        end if
        
        iir= ii + 1 !neighbor-right
        if(iir > 0 .and. nc(1,iir) > 0)then
          do k= 1,nc(1,iir)
            j= ibox(1,iir,k)
            if(typep(j) == 0)then
              dx= x(3) - x(j)
              r= abs(dx)
              hi= h(3)
              
              if(r <= 2.*hi)then
                call kernel(r,dx,hi,w,dwdr,dwdx,3,j)
                
                mj= mass(j)
                Qj= Q(j)
                aj= a(j)
                
                Qi= Qi + mj*Qj*w/aj
                ai= ai + mj*w
                sumw_in=  sumw_in + mj*w/aj
              end if
            end if
          end do
        end if
      end if
      
      Q(3)= Qi/sumw_in
      a(3)= ai/sumw_in
      
      call C_line   ![PhD thesis, chap. 4, eq. (4.3)]
      
      d(1)= dinf
c     a(1)= AP(dinf,Bwid(1),slopm,1)
      a(1)= AP(dinf,Dia,1)
      u(1)= Qinf/a(1)
      !*
      
      
      
      !~ sumd(0),apha(0),resd(0)
      !*******************************************
      do i= nvirp1,np
        apha(i)= 0.
        sumA(i)= 0.
      end do
      
      do l= 1,nct
        if(nc(1,l) > 0)then !self cell
          do k1= 1,nc(1,l)-1
            i= ibox(1,l,k1)
            do k2= k1+1,nc(1,l)
              j= ibox(1,l,k2)
              if(typep(i) >= -1 .and. typep(j) >= -1)then
                dx= x(i) - x(j)
                r= abs(dx)
                hlar= dmax1(h(i),h(j))
                
                if(r <= 2.*hlar)then
                  call kernel(r,dx,h(i),w,dwdr,dwdx,i,j)
                  apha(i)= apha(i) - mass(j)*r*dwdr    
                  sumA(i)= sumA(i) + mass(j)*w
                  call kernel(r,-dx,h(j),w,dwdr,dwdx,j,i)
                  apha(j)= apha(j) - mass(i)*r*dwdr    
                  sumA(j)= sumA(j) + mass(i)*w
                end if
              end if
            end do
          end do
        end if
        
        if(nc(1,l) > 0 .and. (l+1) <= nct .and.
     &    nc(1,l+1) > 0)then !neighbor cell
          do k1= 1,nc(1,l)
            i= ibox(1,l,k1)
            do k2= 1,nc(1,l+1)
              j= ibox(1,l+1,k2)
              if(typep(i) >= -1 .and. typep(j) >= -1)then
                dx= x(i) - x(j)
                r= abs(dx)
                hlar= dmax1(h(i),h(j))
                
                if(r <= 2.*hlar)then
                  call kernel(r,dx,h(i),w,dwdr,dwdx,i,j)
                  apha(i)= apha(i) - mass(j)*r*dwdr    
                  sumA(i)= sumA(i) + mass(j)*w
                  call kernel(r,-dx,h(j),w,dwdr,dwdx,j,i)
                  apha(j)= apha(j) - mass(i)*r*dwdr    
                  sumA(j)= sumA(j) + mass(i)*w
                end if
              end if
            end do
          end do
        end if
      end do
      
      do i= nvirp1,np
        if(typep(i) == 0)then
          call kernel(0.,0.,h(i),w,dwdr,dwdx,i,i)
          sumA(i)= sumA(i) + mass(i)*w
          resA(i)= a(i) - sumA(i)
        end if
      end do
      !*******************************************
      
      
      
      !~ area(k+1),h(k+1)
      !*******************************************
      do while(nn <= 19)
        nn= nn+1
        index_conv= 1
        
        !~ area(k+1),h(k+1)
        do i= nvirp1,np
          if(typep(i) == 0)then
            Aprev(i)= a(i)
            ecrit= abs(resA(i))/Aprev(i)
            
            !~ check if convergence
            if(ecrit .ge. 1.e-10)then
              index_conv= index_conv*0
            else
              index_conv= index_conv*1
            end if
            !*
            
            a(i)= Aprev(i)*(1.-resA(i)/(resA(i)+apha(i)))
            h(i)= hini(i)*aini(i)/a(i)
            
            !~ check smoothing length
            if(h(i) > hmax)then
              write(24,'(3(i6,1x),5(e12.5,1x))')
     &          itime,nn,i,h(i),a(i),x(i),Aprev(i),apha(i)
              pause 'smoothing length over the maximum !!!'
            end if
            !*
          end if
        end do
        !*
        
        if(index_conv == 1) exit       
        
        !~ resA(k),sumA(k),apha(k)
        do i= nvirp1,np
        	apha(i)= 0.
          sumA(i)= 0.
        end do
        
        do l= 1,nct
          if(nc(1,l) > 0)then !self cell
            do k1= 1,nc(1,l)-1
              i= ibox(1,l,k1)
              do k2= k1+1,nc(1,l)
                j= ibox(1,l,k2)
                if(typep(i) >= -1 .and. typep(j) >= -1)then
                  dx= x(i) - x(j)
                  r= abs(dx)
                  hlar= dmax1(h(i),h(j))
                  
                  if(r <= 2.*hlar)then
                    call kernel(r,dx,h(i),w,dwdr,dwdx,i,j)
                    apha(i)= apha(i) - mass(j)*r*dwdr    
                    sumA(i)= sumA(i) + mass(j)*w
                    call kernel(r,-dx,h(j),w,dwdr,dwdx,j,i)
                    apha(j)= apha(j) - mass(i)*r*dwdr    
                    sumA(j)= sumA(j) + mass(i)*w
                  end if
                end if
              end do
            end do
          end if
          
          if(nc(1,l) > 0 .and. (l+1) <= nct .and.
     &      nc(1,l+1) > 0)then !neighbor cell
            do k1= 1,nc(1,l)
              i= ibox(1,l,k1)
              do k2= 1,nc(1,l+1)
                j= ibox(1,l+1,k2)
                if(typep(i) >= -1 .and. typep(j) >= -1)then
                  dx= x(i) - x(j)
                  r= abs(dx)
                  hlar= dmax1(h(i),h(j))
                  
                  if(r <= 2.*hlar)then
                    call kernel(r,dx,h(i),w,dwdr,dwdx,i,j)
                    apha(i)= apha(i) - mass(j)*r*dwdr    
                    sumA(i)= sumA(i) + mass(j)*w
                    call kernel(r,-dx,h(j),w,dwdr,dwdx,j,i)
                    apha(j)= apha(j) - mass(i)*r*dwdr    
                    sumA(j)= sumA(j) + mass(i)*w
                  end if
                end if
              end do
            end do
          end if
        end do
        
        do i= nvirp1,np
          if(typep(i) == 0)then
            call kernel(0.,0.,h(i),w,dwdr,dwdx,i,i)
            sumA(i)= sumA(i) + mass(i)*w
            resA(i)= a(i) - sumA(i)
          end if
        end do
        !*
      end do
      !*******************************************
      
      !Updating surcharge state
	do i= nvirp1,np
	call state(i)
	enddo
      
	!updating ventilation state
      call check_ventil

	!for pipe, give A, compute theta
      do i= nvirp1,np  
	 if(vent(i) ==1) then !(pressurised/depression)
	 theta(i)= 2*pi
	 else   !(free surface)
       theta(i)= rtnewt_theta(a(i), Dia, Amax)
	 endif
      end do

	!computing Hs pressure head
	do i= nvirp1,np  
	 h_pressure(i) = pressure_head(i) 
      enddo

      
      !~ inflow & outflow
      do i= nvirp1,np
        if(typep(i) == 0)then
          if(slopm <= 1.e-6)then
            d(i)= a(i)/Bwid(i) 
          else
c            d(i)= 0.5*( -Bwid(i) + 
c     &             dsqrt(Bwid(i)*Bwid(i) + 4.*slopm*a(i)) )
c     &             /slopm
c           theata(i)= rtnewt_theta(a(i), Dia, Amax) !for pipe
            d(i)= 0.5*(1-cos(0.5*theta(i)))*Dia
	        if(d(i)>Dia) then
	         write(*,*) 'd=',d(i),' dia=',Dia
	         if(d(i) <= Dia*1.001) then
	          d(i) = Dia
	          write(*,*) 'd~Dia * 1.001'
	         else
	          pause 'd(i) > Dia'
	         endif
	        endif	      
          end if
        else if(typep(i) == 1)then  !outflow particles 
c         a(i)= AP(doutf,Bwid(2),slopm,1)
          a(i)= AP(doutf,Dia,1) 
          theta(i)= rtnewt_theta(a(i), Dia, Amax)
          d(i)= doutf
          h(i)= hini(i)*aini(i)/a(i)
        else if(typep(i) == -1)then  !inflow particles
c         a(i)= AP(dinf,Bwid(1),slopm,1)
          a(i)= AP(dinf,Dia,1)   
          theta(i)= rtnewt_theta(a(i), Dia, Amax)    
          aini(i)= a(i)
          mass(i)= dx0*a(i)
          
          h(i)= 1.3*dx0
          hini(i)= h(i)
          d(i)= dinf
        end if
      end do
      !*
      
      
      write(24,'(2(i6,1x))')itime,nn
      
      if(nn > 20) write(24,'(2(i6,1x))')itime,nn
                  
      
      end