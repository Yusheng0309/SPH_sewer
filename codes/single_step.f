      !~ 2012/09
      subroutine single_step(icell,jcell,kind_cell)
      
	use AP_comput      
      implicit none
      
      include 'common_global.dat'
      
      real*8:: dx2,du,dot,dpress
      
      open(25,file= 'step.out')
      
      
      dx2= 1.e-6
      
      if(kind_cell == 1)then !neighbor
        do k1= 1,nc(1,icell)
          i= ibox(1,icell,k1)
          do k2= 1,nc(1,jcell)
            j= ibox(1,jcell,k2)
            
            if(typep(i) >= -1 .and. typep(j) >= -1)then
              dx= x(i) - x(j)
              r= abs(dx)
              hlar= dmax1(h(i),h(j))
              
              wbar= 0.
              dwdxbar= 0.
              
              if(r <= 2.*hlar)then        
                call kernel(r,dx,h(i),w,dwdr,dwdx,i,j)
                wbar= wbar + w
                dwdxbar= dwdxbar + dwdx
                call kernel(r,dx,h(j),w,dwdr,dwdx,i,j)
                wbar= wbar + w
                dwdxbar= dwdxbar + dwdx
                                
                wbar= 0.5*wbar
                dwdxbar= 0.5*dwdxbar
                
                mi= mass(i)
                mj= mass(j)
                
                Qi= Q(i) 
                Qj= Q(j)
                
                     !if flow pressurised then use Apipe in M. eqn.
                ai= min(a(i), Amax) 
                aj= min(a(j), Amax) 
                abar= 0.5*(ai+aj)
                                                
                    ! Hs + dia or  Hs + dia/2 ?               
                di= d(i) + h_pressure(i)
                dj= d(j) + h_pressure(j)
                cbar= 0.5*(sqrt(g*di) + sqrt(g*dj))

                ui= Qi/ai
                uj= Qj/aj
                du= ui - uj 
                
                !~ advection
                dQdt(i)= dQdt(i) + ui*mj*du*dwdxbar
                dQdt(j)= dQdt(j) + uj*mi*du*dwdxbar
                !*
                
			                                                                                                                            
                !~ pressure   [paper, 2005, eq. (24)]
                dpress= di/ai/ai + dj/aj/aj
                dQdt(i)= dQdt(i) - g*ai*ai*mj*dpress*dwdxbar
                dQdt(j)= dQdt(j) + g*aj*aj*mi*dpress*dwdxbar
                !*
                                                                                                          
                !~ artificial viscosity [paper, 2005, eq. (26)]
                dot= du*dx
                dQdt(i)= dQdt(i) + abar*
     &            mj*cbar*(dot)*dwdxbar/sqrt(r*r+dx2)/aj
                dQdt(j)= dQdt(j) - abar*
     &            mi*cbar*(dot)*dwdxbar/sqrt(r*r+dx2)/ai
                !*
              end if
            end if
          end do
        end do
      else if(kind_cell == 2)then !self
        do k1= 1,nc(1,icell)-1
          i= ibox(1,icell,k1)
          do k2= k1+1,nc(1,icell)
            j= ibox(1,icell,k2)
            
            if(typep(i) >= -1 .and. typep(j) >= -1)then
              dx= x(i) - x(j)
              r= abs(dx)
              hlar= dmax1(h(i),h(j))
              wbar= 0.
              dwdxbar= 0.
              
              if(r <= 2.*hlar)then        
                call kernel(r,dx,h(i),w,dwdr,dwdx,i,j)
                wbar= wbar + w
                dwdxbar= dwdxbar + dwdx
                call kernel(r,dx,h(j),w,dwdr,dwdx,i,j)
                wbar= wbar + w
                dwdxbar= dwdxbar + dwdx
                                
                wbar= 0.5*wbar
                dwdxbar= 0.5*dwdxbar
                
                mi= mass(i)
                mj= mass(j)
                
                Qi= Q(i) 
                Qj= Q(j)
                             !if flow pressurised then use Apipe in M. eqn.
                ai= min(a(i), Amax) 
                aj= min(a(j), Amax)  
                abar= 0.5*(ai+aj)
                              ! Hs + dia or  Hs + dia/2 ?               
                di= d(i) + h_pressure(i)
                dj= d(j) + h_pressure(j)
                cbar= 0.5*(sqrt(g*d(i)) + sqrt(g*d(j)))

                ui= Qi/ai
                uj= Qj/aj
                du= ui - uj 
                                
                !~ advection
                dQdt(i)= dQdt(i) + ui*mj*du*dwdxbar
                dQdt(j)= dQdt(j) + uj*mi*du*dwdxbar
                !*
                                                                                                                                         
                !~ pressure   [paper, 2005, eq. (24)]
                dpress= di/ai/ai + dj/aj/aj
                dQdt(i)= dQdt(i) - g*ai*ai*mj*dpress*dwdxbar
                dQdt(j)= dQdt(j) + g*aj*aj*mi*dpress*dwdxbar
                !*
                                                                                                          
                !~ artificial viscosity [paper, 2005, eq. (26)]
                dot= du*dx
                dQdt(i)= dQdt(i) + abar*
     &            mj*cbar*(dot)*dwdxbar/sqrt(r*r+dx2)/aj
                dQdt(j)= dQdt(j) - abar*
     &            mi*cbar*(dot)*dwdxbar/sqrt(r*r+dx2)/ai
                !*
              end if
            end if
          end do
        end do
      end if

            
      end