      subroutine calc_vi_filament(xp,yp,zp,vix,viy,viz)
      use wake_variables
      use mparam
      use flowfp
      use input
      use xtra      
      implicit none
      doubleprecision :: xp,yp,zp,vix,viy,viz
      
c  local
      doubleprecision :: dx1,dy1,dz1,dx2,dy2,dz2,num,dem,r1s,r2s,rcs,
     &  ls,rl,es
      doubleprecision :: rx1,rx2,r1r2,r1xr2,r1yr2,r1zr2,gcm,dvix,dviy,
     & dviz,rxrs
      integer j,ll,nm,m

      es  = 0.000001    
      ! to desingularize the kernal for 0 segment length
     
      vix = 0.0
      viy = 0.0
      viz = 0.0
      
!!!!!!Precompute the position vectors for both trailing and shed wake
        
      do m = 1,nrotor       
        do ll = 1,nblds
          do nm = 1,ntrailers
            do j=1,nfil(nm)

              dx1 = xp - xm(j,nm,ll,m)
              dy1 = yp - ym(j,nm,ll,m)
              dz1 = zp - zm(j,nm,ll,m)
          
              dxarr(j,nm) = dx1
              dyarr(j,nm) = dy1
              dzarr(j,nm) = dz1

              r1s = dx1 * dx1 + dy1 * dy1 + dz1 * dz1
      
              if( r1s .le. 0.0) then
              rx1 = 1.0 
        ! trick to deal with / by 0 later since r1 x r2 will be zero!
            else
                rx1  = dsqrt( r1s )
            endif
      
              rsarr(j,nm) = r1s
              rarr(j,nm) = rx1

            enddo
          enddo      
         
!!!!!!!!!! ----- Effect of Trailing wake 
            
          do nm = 1,ntrailers

          dx1 = dxarr(1,nm)
          dy1 = dyarr(1,nm)
          dz1 = dzarr(1,nm)
          r1s = rsarr(1,nm)
          rx1 = rarr(1,nm)
      
           do j=1,nfil(nm)-1

           dx2  = dxarr(j+1,nm)
           dy2  = dyarr(j+1,nm)
           dz2  = dzarr(j+1,nm)
           r2s  = rsarr(j+1,nm)
           rx2  = rarr(j+1,nm)
        
           rcs   = rc(j,nm,ll,m) * rc(j,nm,ll,m)
           
!c  dot product
             r1r2  = dx1 * dx2 + dy1 * dy2 + dz1 * dz2
!c  cross product
             r1xr2 = - ( dy2 * dz1 - dy1 * dz2 )
             r1yr2 = - ( dx1 * dz2 - dx2 * dz1 )
             r1zr2 = - ( dx2 * dy1 - dx1 * dy2 )

             ls    = r1s + r2s - 2.0 * r1r2
             rl    = rc(j,nm,ll,m) * rc(j,nm,ll,m) * ( ls + es )
             rxrs  = r1s * r2s - r1r2 * r1r2 

             num   =  ( rx1 + rx2 ) * ( rx1 * rx2 - r1r2 )
             dem   = rx1 * rx2 * dsqrt( rxrs * rxrs + rl * rl )
          
             gcm = gamu(j,nm,ll,m) * float(bld_inside(j,nm,ll,m)) 
     &                               * num / dem

             dvix  = gcm * r1xr2
             dviy  = gcm * r1yr2
             dviz  = gcm * r1zr2
          
             vix   = vix + dvix
             viy   = viy + dviy
             viz   = viz + dviz
         
!c save current node #2 geometry as node #1 geometry

           dx1   = dx2
           dy1   = dy2
           dz1   = dz2
           r1s   = r2s
           rx1    = rx2
              
           enddo
          enddo ! end of multiple trailers      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!Effect of Shed wake 
          if(rshedwakerev.gt.0) then      
          do j = 2,nmarksh

          dx1 = dxarr(j,1)
          dy1 = dyarr(j,1)
          dz1 = dzarr(j,1)
          r1s = rsarr(j,1)
          rx1  = rarr(j,1)
         
            do nm = 1,ntrailers-1

            dx2   = dxarr(j,nm+1)
            dy2   = dyarr(j,nm+1)
            dz2   = dzarr(j,nm+1)
            r2s   = rsarr(j,nm+1)
            rx2    = rarr(j,nm+1)
        
!c  dot product
              r1r2  = dx1 * dx2 + dy1 * dy2 + dz1 * dz2
!c  cross product
              r1xr2 = - ( dy2 * dz1 - dy1 * dz2 )
              r1yr2 = - ( dx1 * dz2 - dx2 * dz1 )
              r1zr2 = - ( dx2 * dy1 - dx1 * dy2 )

              ls    = r1s + r2s - 2.0 * r1r2
              rl    = rcshed(j,nm,ll,m) * rcshed(j,nm,ll,m)
     &                      * ( ls + es )
              rxrs  = r1s * r2s - r1r2 * r1r2 

              num   =  ( rx1 + rx2 ) * ( rx1 * rx2 - r1r2 )
              dem   = rx1 * rx2 * dsqrt( rxrs * rxrs + rl * rl )
          
          
      gcm   = gamushed(j,nm,ll,m) * float(bld_inside(j,nm,ll,m))
     &                                    * num / dem

              dvix  = gcm * r1xr2
              dviy  = gcm * r1yr2
              dviz  = gcm * r1zr2
          
              vix   = vix + dvix
              viy   = viy + dviy
              viz   = viz + dviz
         
c save current node #2 geometry as node #1 geometry

           dx1   = dx2
           dy1   = dy2
           dz1   = dz2
           r1s   = r2s
             rx1    = rx2
              
           enddo
          enddo 
        endif ! end of shed wake
        
        enddo ! end of blade loop
      enddo ! end of multiple rotor loop
      
      
      return
      end
