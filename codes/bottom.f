      !~2012/09
      subroutine bottom(istart,iend) ![paper, 2011, eq. (21)(22)]
      
c	use AP_comput
c      implicit none
      
      include 'common_global.dat'
      
      real*8:: dx2
      
      open(40,file= 'bottom.out')
      
      
      dx2= 1.e-6
      
      do i= istart,iend
        Bwid(i)= 0.    !be a constant in pipes
        s0(i)= 0.
        nMa(i)= 0.
        sumw_b= dx2
        
        if( typep(i) == 0 .or. typep(i) == -3 )then
          ii= iipar(i) !self
          if( nc(2,ii) > 0 )then
            do k= 1,nc(2,ii)
              j= ibox(2,ii,k)
                      
              dx= x(i) - xvirb(j)
              r= abs(dx)
              hj= hvirb(j) 
              
              if( r <= 2.*hj )then
                call kernel(r,dx,hj,w,dwdr,dwdx,i,j)
                
                mj= massvirb(j)
                
                Bwid(i)= Bwid(i) + mj*Bwidvirb(j)*w !be a constant in pipes
                s0(i)= s0(i) + mj*s0virb(j)*w
                nMa(i)= nMa(i) + mj*nMavirb(j)*w
                sumw_b= sumw_b + mj*w
              end if
            end do
          end if    
                     
          iil= ii - 1 !neighbor-left
          if( iil > 0 .and. nc(2,iil) > 0 )then
            do k= 1,nc(2,iil)
              j= ibox(2,iil,k)
              
              dx= x(i) - xvirb(j)
              r= abs(dx)
              hj= hvirb(j) 
              
              if(r <= 2.*hj)then
                call kernel(r,dx,hj,w,dwdr,dwdx,i,j)
                
                mj= massvirb(j)
                
                Bwid(i)= Bwid(i) + mj*Bwidvirb(j)*w  !be a constant in pipes
                s0(i)= s0(i) + mj*s0virb(j)*w
                nMa(i)= nMa(i) + mj*nMavirb(j)*w
                sumw_b= sumw_b + mj*w
              end if
            end do
          end if
                     
          iir= ii+1 !neighbor-right
          if(iir <= nct .and. nc(2,iir) > 0)then
            do k= 1,nc(2,iir)
              j= ibox(2,iir,k)
                            
              dx= x(i) - xvirb(j)
              r= abs(dx)
              hj= hvirb(j) 
              
              if(r <= 2.*hj)then
                call kernel(r,dx,hj,w,dwdr,dwdx,i,j)
                
                mj= massvirb(j)
                
                Bwid(i)= Bwid(i) + mj*Bwidvirb(j)*w   !be a constant in pipes
                s0(i)= s0(i) + mj*s0virb(j)*w
                nMa(i)= nMa(i) + mj*nMavirb(j)*w
                sumw_b= sumw_b + mj*w
              end if
            end do
          end if
        end if
        
        Bwid(i)= Bwid(i)/sumw_b !be a constant in pipes
        s0(i)= s0(i)/sumw_b
        nMa(i)= nMa(i)/sumw_b
      end do
      
                    
      end