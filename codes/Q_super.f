      !~ 2012/05/15
      subroutine Q_super   ![paper, 2005 2010 2011]
            
      implicit none
      
      include 'common_global.dat'
      
      real*8:: dx2
      
      open(24,file= 'depth.out')
                 
      
      nn= 0
      dx2= 1.e-6
      
      !~ Characteristic method ¡V Qoutf
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
      ai= ai/sumw_out
      Q(2)= Qi 
      a(2)= ai 
      
      Qoutf= Q(2)
      !*
      
      
      end