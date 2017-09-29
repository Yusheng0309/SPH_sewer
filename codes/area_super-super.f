      !~ 2012/03/22
      subroutine SWarea_ss   ![paper, 2005 2010 2011]
      
	use AP_comput      
      implicit none
      
      include 'common_global.dat'
      
      integer:: index_conv   
      real*8:: Aprev(nmax),sumA(nmax),apha(nmax),resA(nmax)
      real*8:: dx2
      
      open(24,file= 'depth.out')
                 
      
      nn= 0
      dx2= 1.e-6
      
      
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
              write(24,'(3(i6,1x),3(e12.5,1x))')
     &          itime,nn,i,h(i),a(i),x(i)
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
      
      
      !~ doutf
      Qi= 0.
      ai= 0.
      sumw_out= dx2
      
      ii= iipar(2) !self
      if(nc(1,ii) > 0)then
        do k= 1,nc(1,ii)
          j= ibox(1,ii,k)
          if(typep(j) == 0)then
            dx= x(2) - x(j)
            r= abs(dx)
            hi= h(2)
            
            if(r <= 2.*hi)then
              call kernel(r,dx,hi,w,dwdr,dwdx,2,j)
              
              mj= mass(j)
              Qj= Q(j)
              aj= a(j)
              
              Qi= Qi + mj*Qj*w/aj
              ai= ai + mj*w
              sumw_out=  sumw_out + mj*w/aj
            end if
          end if
        end do
        
        iil= ii - 1 !neighbor-left
        if(iil > 0 .and. nc(1,iil) > 0)then
          do k= 1,nc(1,iil)
            j= ibox(1,iil,k)
            if(typep(j) == 0)then
              dx= x(2) - x(j)
              r= abs(dx)
              hi= h(2)
              
              if(r <= 2.*hi)then
                call kernel(r,dx,hi,w,dwdr,dwdx,2,j)
                
                mj= mass(j)
                Qj= Q(j)
                aj= a(j)
                
                Qi= Qi + mj*Qj*w/aj
                ai= ai + mj*w
                sumw_out=  sumw_out + mj*w/aj
              end if
            end if
          end do
        end if
        
        iir= ii + 1 !neighbor-right
        if(iir > 0 .and. nc(1,iir) > 0)then
          do k= 1,nc(1,iir)
            j= ibox(1,iir,k)
            if(typep(j) == 0)then
              dx= x(2) - x(j)
              r= abs(dx)
              hi= h(2)
              
              if(r <= 2.*hi)then
                call kernel(r,dx,hi,w,dwdr,dwdx,2,j)
                
                mj= mass(j)
                Qj= Q(j)
                aj= a(j)
                
                Qi= Qi + mj*Qj*w/aj
                ai= ai + mj*w
                sumw_out=  sumw_out + mj*w/aj
              end if
            end if
          end do
        end if
      end if
      
      Qi= Qi/sumw_out
      ai= ai/sumw_out   ![PhD thesis, chap. 4, eq. (4.3)]
      
      
      if(slopm <= 1.e-6)then
        doutf= ai/Bwid(2) 
      else
        doutf= 0.5*( -Bwid(2) + dsqrt(Bwid(2)*Bwid(2) + 4.*slopm*ai) )
     &         /slopm
      end if
      !*
      
      
      
      !~ inflow & outflow
      do i= nvirp1,np
        if(typep(i) == 0)then
          if(slopm <= 1.e-6)then
            d(i)= a(i)/Bwid(i) 
          else
            d(i)= 0.5*( -Bwid(i) + 
     &             dsqrt(Bwid(i)*Bwid(i) + 4.*slopm*a(i)) )
     &             /slopm
          end if
        else if(typep(i) == 1)then
          a(i)= AP(doutf,Bwid(2),slopm,1)
          d(i)= doutf
          h(i)= hini(i)*aini(i)/a(i)
        else if(typep(i) == -1)then
          a(i)= AP(dinf,Bwid(1),slopm,1)
          d(i)= dinf
          h(i)= hini(i)*aini(i)/a(i)
        end if
      end do
      !*
      
      
      write(24,'(2(i6,1x))')itime,nn
      
      if(nn > 20) write(24,'(2(i6,1x))')itime,nn
                  
      
      end