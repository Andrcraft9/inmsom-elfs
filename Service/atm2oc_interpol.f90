!=====================================================================
subroutine weight_matrix_intrp(xin,           & !array(1-D) of input x-grid values (input) 
                               yin,           & !array(1-D) of input y-grid values (input)
                              nxin,           & !number of input x-grid points(input)
                              nyin,           & !number of input y-grid points(input)
                              xout,           & !array(2D) of output x-grid values in input coordinate system (input)
                              yout,           & !array(2D) of output y-grid values in input coordinate system (input)
                             nxout1,nxout2,   & !boundaries of arrays for output x-grid coordinates (input)
                             nyout1,nyout2,   & !boundaries of arrays for output y-grid coordinates (input)
                             maskin,          & !input sea-land mask (input)
                             maskout,         & !output sea-land mask (input)
                           i_input,           & !x-grid numbers of input grid for output grid (output)
                           j_input,           & !y-grid numbers of input grid for output grid (output)
                           matrx_elmnt_intrp, & !nonzero matrix elements of interpolation (output)
                           indper,            & !index of input grid periodicity:=0 -nonperiodic,=1 -periodic case
                           fillmiss,          & !filling missed values (0-no, 1-yes)
                                   mmm_in,    & !first significant point in x-direction for input grid (input)
                                    mm_in,    & ! last significant point in x-direction for input grid  (input)
                                   nnn_in,    & !first significant point in y-direction for input grid  (input)
                                    nn_in,    & ! last significant point in y-direction for input grid  (input)
                                   mmm_out,   & !first significant point in x-direction for output grid  (input)
                                    mm_out,   & ! last significant point in x-direction for output grid  (input)
                                   nnn_out,   & !first significant point in y-direction for output grid  (input)
                                    nn_out)     ! last significant point in y-direction for output grid  (input)
   	
!  definition of data parameters
implicit none
include 'crdfnc.fi'

integer indper !index of periodicity:=0 -nonperiodic,=1 -periodic case

integer nxin,nyin,nxout1,nxout2,nyout1,nyout2,fillmiss

integer maskin(nxin,nyin)
real(4) maskout(nxout1:nxout2,nyout1:nyout2)

integer i_input(nxout1:nxout2,nyout1:nyout2,4),      &  !x-grid numbers of input grid for output grid
        j_input(nxout1:nxout2,nyout1:nyout2,4)          !y-grid numbers of input grid for output grid

real(8) matrx_elmnt_intrp(nxout1:nxout2,nyout1:nyout2,4)

real(8) sum_matr_elmnt

real(8) xin(nxin), yin(nyin)                	!input grid

real(8) xout(nxout1:nxout2,nyout1:nyout2),  &
        yout(nxout1:nxout2,nyout1:nyout2)       !output grid

real(8) delta_left,delta_right

integer i_nes,j_nes                             !auxilary indexes

integer i_left,i_right,j_down,j_up              !auxilary indexes

integer i,j,k,m,n

real(8) bilin_denom

real(8) x_left,x_right,y_down,y_up, ret_lon, ret_lat

real(8), allocatable:: xinp(:),yinp(:),  temp(:)        !recalculated grid for periodic case                

integer mmm_in,  mm_in,  nnn_in,  nn_in, mmm_out, mm_out, nnn_out, nn_out



!---------input grid analysis on periodicity----------------------

      if(indper/=1.and.indper/=0) then
	   write (*,*) 'wrong index of periodicity in subroutine weight_matrix_intrp'
	   stop 
	endif
	
		allocate (xinp(nxin+indper),yinp(nyin),temp(4))
! new lon and lat grid definition

	     xinp(mmm_in:mm_in)=xin(mmm_in:mm_in)
           if(indper==1) xinp(mm_in+1)=xin(mmm_in)+360.0d0 !completing of longitude
	     
		 yinp=yin

! interpolation from input to output grid
	do j=nnn_out,nn_out
	 do i=mmm_out,mm_out
           
	  if (maskout(i,j)>0.5) then
            
            ret_lon=xout(i,j)
            ret_lat=yout(i,j)

!          corrections

            if(indper==1) then
	        if(ret_lon > xinp(mm_in+1)) ret_lon=ret_lon-360.0d0
	        if(ret_lon < xinp(mmm_in) ) ret_lon=ret_lon+360.0d0
            else
              
              if(ret_lon > xinp(mm_in)) then
	         delta_right=ret_lon-xinp(mm_in)
	         ret_lon=ret_lon-360.0d0
               delta_left=xin(mmm_in)-ret_lon
	         if(delta_left > delta_right) ret_lon=ret_lon+360.0d0
	         go to 111
	        end if

              if(ret_lon < xinp(mmm_in)) then
	         delta_left=xinp(mmm_in)-ret_lon
	         ret_lon=ret_lon+360.0d0
 	         delta_right=ret_lon-xinp(mm_in)
	         if(delta_right > delta_left) ret_lon=ret_lon-360.0d0
	         go to 111
	        end if
               	       
            end if

111	 continue
!       identification of reversly rotated point on geografic grid  
       
        i_nes=mmm_in
	 
	   n=mm_in+indper

1000     m=(i_nes+n)/2	 
	     
           if (ret_lon < xinp(m)) n=m
	     if (ret_lon >=xinp(m)) i_nes=m
	     if (n > i_nes+1) goto 1000
	 
         i_left =i_nes
	   i_right=i_nes+1
         if((mm_in-mmm_in+1)==1) i_right=i_nes

	   j_nes=nnn_in
	   n=nn_in
1001     m=(j_nes+n)/2	 
	     
           if (ret_lat < yinp(m)) n=m
	     if (ret_lat >=yinp(m)) j_nes=m
	     if (n > j_nes+1) goto 1001
        
	   j_down = j_nes
         j_up   = j_nes+1
         if((nn_in-nnn_in+1)==1)j_up=j_nes

!       main bilinear interpolation

!--------------definition of output i-indexes and preparing of matrix calculating----

1982     if((i_left < mmm_in)   .and.(indper==1)) then
          
          i_input(i,j,1)= i_left + (mm_in-mmm_in+1)
	    i_input(i,j,3)= i_left + (mm_in-mmm_in+1)

	    x_left = xinp(i_input(i,j,1))-360.

         else
          
	    i_input(i,j,1)=i_left
          i_input(i,j,3)=i_left

          x_left =xinp(i_input(i,j,1))

	   end if
	            
	   
	   if((i_right > mm_in).and.(indper==1)) then
          
		i_input(i,j,2)=i_right-(mm_in-mmm_in+1)
		i_input(i,j,4)=i_right-(mm_in-mmm_in+1)

		x_right=xinp(i_input(i,j,2))+360.
	   
	   else
          
		i_input(i,j,2)=i_right
		i_input(i,j,4)=i_right

          x_right=xinp(i_input(i,j,2))
	   
	   end if

!--------------definition of output j-indexes--------------
	   j_input(i,j,1)=j_down
	   j_input(i,j,2)=j_down
	   j_input(i,j,3)=j_up
	   j_input(i,j,4)=j_up

	   y_down =yinp(j_input(i,j,1))
	   y_up   =yinp(j_input(i,j,3))

!---------if all 4 points are undefined moving bounds of rectangle	   
	   if((maskin(i_input(i,j,1),j_input(i,j,1))*        &
             maskin(i_input(i,j,2),j_input(i,j,2))*        &
             maskin(i_input(i,j,3),j_input(i,j,3))*        &
             maskin(i_input(i,j,4),j_input(i,j,4)))/=0) then

           if(fillmiss==0) then
           matrx_elmnt_intrp(i,j,1:4)=0.25d0

	     go to 501
	     end if

           if (indper==0) then

              if ((i_left==mmm_in).and.(i_right==mm_in).and.      & 
                  (j_down==nnn_in).and.(j_up   ==nn_in)) go to 1941

	        if ((i_left==mmm_in).and.(i_right==mm_in)) then
	          if (j_down >= (nnn_in+1)) j_down=j_down-1
	          if (j_up   <= ( nn_in-1)) j_up  =j_up  +1
	          i_left =i_nes+1
	          i_right=i_nes
	        end if

		  if (i_left >=(mmm_in+1)) i_left =i_left -1
	        if (i_right<=( mm_in-1)) i_right=i_right+1
		 
		 else
              
			if (((i_right-i_left) >= (mm_in-mmm_in+1)).and.     &
                       (j_down==nnn_in).and.(j_up   ==nn_in)) go to 1941		 
		 				        
              if((i_right-i_left)>=(mm_in-mmm_in+1)) then
	          if (j_down >=(nnn_in+1)) j_down=j_down-1
	          if (j_up   <=( nn_in-1)) j_up  =j_up  +1
	          i_left =i_nes+1
	          i_right=i_nes
	        end if  	        
			
              if((i_right-i_left) < (mm_in-mmm_in+1)) then
                  i_left =i_left -1
                  i_right=i_right+1
	        end if

		 end if
	     
		 go to 1982
	   
	   end if
!---------end of moving bounds of rectangle--------------------------

	   bilin_denom= 1.

	   if (x_right /= x_left) bilin_denom=bilin_denom / (fnclon(x_right)-fnclon(x_left))
	   if (y_up    /= y_down) bilin_denom=bilin_denom / (fnclat(y_up)   -fnclat(y_down))
	 	 
	   temp=bilin_denom


	    if (x_right /=x_left) then

           temp(1)= temp(1)*(fnclon(x_right) - fnclon(ret_lon))
           temp(2)= temp(2)*(fnclon(ret_lon) - fnclon(x_left) )
	     temp(3)= temp(3)*(fnclon(x_right) - fnclon(ret_lon)) 
	     temp(4)= temp(4)*(fnclon(ret_lon) - fnclon(x_left) )

          end if

	    if (y_up    /=y_down) then

           temp(1)= temp(1)*(fnclat(y_up)    - fnclat(ret_lat))
           temp(2)= temp(2)*(fnclat(y_up)    - fnclat(ret_lat))
	     temp(3)= temp(3)*(fnclat(ret_lat) - fnclat(y_down) ) 
	     temp(4)= temp(4)*(fnclat(ret_lat) - fnclat(y_down) )

          end if

!---------- procedure if 4 points are defined---------	  
	    if((maskin(i_input(i,j,1),j_input(i,j,1))+       &
              maskin(i_input(i,j,2),j_input(i,j,2))+       &
              maskin(i_input(i,j,3),j_input(i,j,3))+       &
              maskin(i_input(i,j,4),j_input(i,j,4)))==0) then
	     	
			do k=1,4
		     matrx_elmnt_intrp(i,j,k)=temp(k)
	        end do
!----------end of procedure if 4 points are defined---------

	    else
!---------- procedure if 3 points are defined---------
	       if((maskin(i_input(i,j,1),j_input(i,j,1))+     &
                 maskin(i_input(i,j,2),j_input(i,j,2))+     &
                 maskin(i_input(i,j,3),j_input(i,j,3))+     &
                 maskin(i_input(i,j,4),j_input(i,j,4)))==1) then

	         !the 1-st point is undefined
		     if(maskin(i_input(i,j,1),j_input(i,j,1))==1) then
	            matrx_elmnt_intrp(i,j,1)=0.0d0
	            matrx_elmnt_intrp(i,j,2)=temp(2)+temp(1)
	            matrx_elmnt_intrp(i,j,3)=temp(3)+temp(1)
	            matrx_elmnt_intrp(i,j,4)=temp(4)-temp(1)
               end if
	 
	         !the 2-nd point is undefined
	         if(maskin(i_input(i,j,2),j_input(i,j,2))==1) then
	            matrx_elmnt_intrp(i,j,1)=temp(1)+temp(2)
	            matrx_elmnt_intrp(i,j,2)=0.0d0
	            matrx_elmnt_intrp(i,j,3)=temp(3)-temp(2)
	            matrx_elmnt_intrp(i,j,4)=temp(4)+temp(2)
	         end if
	
	         !the 3-rd point is undefined
	         if(maskin(i_input(i,j,3),j_input(i,j,3))==1) then
	            matrx_elmnt_intrp(i,j,1)=temp(1)+temp(3)
	            matrx_elmnt_intrp(i,j,2)=temp(2)-temp(3)
	            matrx_elmnt_intrp(i,j,3)=0.0d0
	            matrx_elmnt_intrp(i,j,4)=temp(4)+temp(3)
	         end if

	         !the 4-th point is undefined	       
		     if(maskin(i_input(i,j,4),j_input(i,j,4))==1) then
	            matrx_elmnt_intrp(i,j,1)=temp(1)-temp(4)
	            matrx_elmnt_intrp(i,j,2)=temp(2)+temp(4)
	            matrx_elmnt_intrp(i,j,3)=temp(3)+temp(4)
	            matrx_elmnt_intrp(i,j,4)=0.0d0
		     end if 
!----------end of procedure if 3 points are defined---------
		 
		   else
!---------- procedure if 2 points are defined---------
	          if((maskin(i_input(i,j,1),j_input(i,j,1))+          &
                    maskin(i_input(i,j,2),j_input(i,j,2))+          &
                    maskin(i_input(i,j,3),j_input(i,j,3))+          &
                    maskin(i_input(i,j,4),j_input(i,j,4)))==2) then
	
	         !the 1-st and the 2-nd points are undefined		   
		     if((maskin(i_input(i,j,1),j_input(i,j,1))==1).and.     &
                    (maskin(i_input(i,j,2),j_input(i,j,2))==1)) then
	             matrx_elmnt_intrp(i,j,1)=0.0d0
	             matrx_elmnt_intrp(i,j,2)=0.0d0
	             matrx_elmnt_intrp(i,j,3)=temp(1)+temp(3)
	             matrx_elmnt_intrp(i,j,4)=temp(2)+temp(4)
	         end if

	         !the 1-st and the 3-rd points are undefined
		     if((maskin(i_input(i,j,1),j_input(i,j,1))==1).and.     &
                    (maskin(i_input(i,j,3),j_input(i,j,3))==1)) then
	             matrx_elmnt_intrp(i,j,1)=0.0d0
	             matrx_elmnt_intrp(i,j,2)=temp(1)+temp(2)
	             matrx_elmnt_intrp(i,j,3)=0.0d0
	             matrx_elmnt_intrp(i,j,4)=temp(3)+temp(4)
	         end if

	         !the 3-rd and the 4-th points are undefined		   
		     if((maskin(i_input(i,j,3),j_input(i,j,3))==1).and.     &
                    (maskin(i_input(i,j,4),j_input(i,j,4))==1)) then
	             matrx_elmnt_intrp(i,j,1)=temp(1)+temp(3)
	             matrx_elmnt_intrp(i,j,2)=temp(2)+temp(4)
	             matrx_elmnt_intrp(i,j,3)=0.0d0
	             matrx_elmnt_intrp(i,j,4)=0.0d0
	         end if

	         !the 4-th and the 2-nd points are undefined
		     if((maskin(i_input(i,j,4),j_input(i,j,4))==1).and.     & 
                    (maskin(i_input(i,j,2),j_input(i,j,2))==1)) then
	             matrx_elmnt_intrp(i,j,1)=temp(1)+temp(2)
	             matrx_elmnt_intrp(i,j,2)=0.0d0
	             matrx_elmnt_intrp(i,j,3)=temp(3)+temp(4)
	             matrx_elmnt_intrp(i,j,4)=0.0d0
	         end if

	         !the 1-st and the 4-th points are undefined
		     if((maskin(i_input(i,j,1),j_input(i,j,1))==1).and.     &
                    (maskin(i_input(i,j,4),j_input(i,j,4))==1)) then
	             matrx_elmnt_intrp(i,j,1)=0.0d0
	             matrx_elmnt_intrp(i,j,2)=temp(1)+temp(2)
	             matrx_elmnt_intrp(i,j,3)=temp(3)+temp(4)
	             matrx_elmnt_intrp(i,j,4)=0.0d0
	         end if

	         !the 2-nd and the 3-rd points are undefined
		     if((maskin(i_input(i,j,3),j_input(i,j,3))==1).and.     &
                    (maskin(i_input(i,j,2),j_input(i,j,2))==1)) then
                   matrx_elmnt_intrp(i,j,1)=temp(1)+temp(2)
	             matrx_elmnt_intrp(i,j,2)=0.0d0
	             matrx_elmnt_intrp(i,j,3)=0.0d0
	             matrx_elmnt_intrp(i,j,4)=temp(3)+temp(4)
	         end if

!----------end of procedure if 2 points are defined---------	         
			 
			  else
!---------- procedure if 1 points is defined---------	
		 				
	          do k=1,4
		       if(maskin(i_input(i,j,k),j_input(i,j,k))==0) then
	             matrx_elmnt_intrp(i,j,k)=1.0d0
	           else
	             matrx_elmnt_intrp(i,j,k)=0.0d0
	           end if
	          end do

!----------end of procedure if 2 points are defined---------				 
			  end if
             end if
	    end if

	  else
         do k=1,4
          matrx_elmnt_intrp(i,j,k)=0.0d0
          i_input(i,j,k)=0
          j_input(i,j,k)=0
         end do
	  end if
      
      matrx_elmnt_intrp(I,J,:)=min(matrx_elmnt_intrp(i,j,:),1.0d0)
      matrx_elmnt_intrp(I,J,:)=max(matrx_elmnt_intrp(i,j,:),0.0d0)  
!--------remorming of matrix elements-------------------------------
500       sum_matr_elmnt=matrx_elmnt_intrp(i,j,1)  &
                        +matrx_elmnt_intrp(i,j,2)  &
                        +matrx_elmnt_intrp(i,j,3)  &
                        +matrx_elmnt_intrp(i,j,4)

	  matrx_elmnt_intrp(i,j,:)= matrx_elmnt_intrp(i,j,:)/sum_matr_elmnt

501    continue     	  	  
	 end do
	end do
		
	deallocate (temp,yinp,xinp)
   
      return

1941  write(*,*) 'entire grid undefined!!!'
      stop
endsubroutine weight_matrix_intrp

!===================================================================================
! scalar field interpolation from atm to ocean
subroutine interpolrot_scal(     nxin,          & !number of x-grid points(input)
                                 nyin,          & !number of y-grid points(input)
                                nxout1,nxout2,  & !number of x-grid points(output)
                                nyout1,nyout2,  & !number of y-grid points(output)
                               funcin,          & !input data array
                              funcout,          & !output data array
                              i_input,          & !x-grid numbers of input grid for output grid
                              j_input,          & !y-grid numbers of input grid for output grid
                              coefficient,      & !nonzero matrix elements of interpolation
                              maskout,          & !sea-land mask for output data
                              undefout,         &
                                 mmm_out,       & !first significant point in x-direction (output)
                                 mm_out,        & ! last significant point in x-direction (output)
                                 nnn_out,       & !first significant point in y-direction (output)
                                 nn_out)          ! last significant point in y-direction (output)
implicit none

integer nxin,nyin,nxout1,nxout2,nyout1,nyout2
integer mmm_out, mm_out, nnn_out, nn_out

real(4) funcin(nxin,nyin)
real(8) funcout(nxout1:nxout2,nyout1:nyout2)

real(8) coefficient(nxout1:nxout2,nyout1:nyout2,4),undefout

integer i_input(nxout1:nxout2,nyout1:nyout2,4),     &     !x-grid numbers of input grid for output grid
        j_input(nxout1:nxout2,nyout1:nyout2,4)            !y-grid numbers of input grid for output grid
real(4) maskout(nxout1:nxout2,nyout1:nyout2)

integer i,j

!$omp parallel do private(i,j)
	 do j=nnn_out,nn_out
	  do i=mmm_out,mm_out
         if (maskout(i,j)<=0.5) then
          funcout(i,j)=undefout
	   else

	    funcout(i,j)= coefficient(i,j,1)*dble(funcin(i_input(i,j,1),j_input(i,j,1)))+   &
                        coefficient(i,j,2)*dble(funcin(i_input(i,j,2),j_input(i,j,2)))+   &
                        coefficient(i,j,3)*dble(funcin(i_input(i,j,3),j_input(i,j,3)))+   &
                        coefficient(i,j,4)*dble(funcin(i_input(i,j,4),j_input(i,j,4)))
         end if
	  end do
	 end do
!$omp end parallel do

endsubroutine interpolrot_scal
!===================================================================
! vector field interpolation from atm to ocean
subroutine interpolrot_vec(nxin,         &  !number of x-grid points(input)
                          nyin,          &  !number of y-grid points(input)
                         nxout1,nxout2,  &  !number of x-grid points(output)
                         nyout1,nyout2,  &  !number of y-grid points(output)
                        funcin_z,        &  !input zonal data array
                        funcin_m,        &  !input meridional data array
                       funcout_z,        &  !output zonal data array
                       funcout_m,        &  !output meridional data array
                       i_input,          &  !x-grid numbers of input grid for output grid
                       j_input,          &  !y-grid numbers of input grid for output grid
                       coefficient,      &  !nonzero matrix elements of interpolation
                       rotvec_coeff,     &  !angles between parallels
                       maskout,          &  !sea-land mask for output data
                       undefout,         &
                          mmm_out,       &  !first significant point in x-direction (output)
                          mm_out,        &  ! last significant point in x-direction (output)
                          nnn_out,       &  !first significant point in y-direction (output)
                          nn_out)           ! last significant point in y-direction (output)
implicit none

integer nxin,nyin,nxout1,nxout2,nyout1,nyout2
integer mmm_out, mm_out, nnn_out, nn_out

real(4) funcin_z(nxin,nyin),  funcin_m(nxin,nyin)
real(8) funcout_z(nxout1:nxout2,nyout1:nyout2),funcout_m(nxout1:nxout2,nyout1:nyout2),     &
        func_zon_unrot,func_mer_unrot

real(8) coefficient(nxout1:nxout2,nyout1:nyout2,4),     &
        rotvec_coeff(nxout1:nxout2,nyout1:nyout2,4), undefout

integer i_input(nxout1:nxout2,nyout1:nyout2,4),       &  !x-grid numbers of input grid for output grid
        j_input(nxout1:nxout2,nyout1:nyout2,4)           !y-grid numbers of input grid for output grid
real(4) maskout(nxout1:nxout2,nyout1:nyout2)

integer i,j

!$omp parallel do private(i,j,func_zon_unrot,func_mer_unrot)
	 do j=nnn_out,nn_out
	  do i=mmm_out,mm_out
         if (maskout(i,j).le.0.5) then
          funcout_z(i,j)=undefout
	    funcout_m(i,j)=undefout
	   else
          func_zon_unrot= coefficient(i,j,1)*dble(funcin_z(i_input(i,j,1),j_input(i,j,1)))+      &
                          coefficient(i,j,2)*dble(funcin_z(i_input(i,j,2),j_input(i,j,2)))+      &
                          coefficient(i,j,3)*dble(funcin_z(i_input(i,j,3),j_input(i,j,3)))+      &
                          coefficient(i,j,4)*dble(funcin_z(i_input(i,j,4),j_input(i,j,4)))
	    
	    func_mer_unrot= coefficient(i,j,1)*dble(funcin_m(i_input(i,j,1),j_input(i,j,1)))+    &
                          coefficient(i,j,2)*dble(funcin_m(i_input(i,j,2),j_input(i,j,2)))+    &
                          coefficient(i,j,3)*dble(funcin_m(i_input(i,j,3),j_input(i,j,3)))+    &
                          coefficient(i,j,4)*dble(funcin_m(i_input(i,j,4),j_input(i,j,4)))

      	funcout_z(i,j)=rotvec_coeff(i,j,1)*func_zon_unrot+     &
                           rotvec_coeff(i,j,2)*func_mer_unrot

      	funcout_m(i,j)=rotvec_coeff(i,j,3)*func_zon_unrot+     & 
                           rotvec_coeff(i,j,4)*func_mer_unrot

         end if
	  end do
	 end do
!$omp end parallel do

endsubroutine interpolrot_vec
!=====================================================================
subroutine weight_matrix_intrp_next(xin,           & !array(1-D) of input x-grid values (input) 
                                    yin,           & !array(1-D) of input y-grid values (input)
                                   nxin,           & !number of input x-grid points(input)
                                   nyin,           & !number of input y-grid points(input)
                                   xout,           & !array(2D) of output x-grid values in input coordinate system (input)
                                   yout,           & !array(2D) of output y-grid values in input coordinate system (input)
                                  nxout1,nxout2,   & !boundaries of arrays for output x-grid coordinates (input)
                                  nyout1,nyout2,   & !boundaries of arrays for output y-grid coordinates (input)
                                  maskin,          & !input sea-land mask (input)
                                  maskout,         & !output sea-land mask (input)
                                i_input,           & !x-grid numbers of input grid for output grid (output)
                                j_input,           & !y-grid numbers of input grid for output grid (output)
                                matrx_elmnt_intrp, & !nonzero matrix elements of interpolation (output)
                                indper,            & !index of input grid periodicity:=0 -nonperiodic,=1 -periodic case
                                fillmiss,          & !filling missed values (0-no, 1-yes)
                                        mmm_in,    & !first significant point in x-direction for input grid (input)
                                         mm_in,    & ! last significant point in x-direction for input grid  (input)
                                        nnn_in,    & !first significant point in y-direction for input grid  (input)
                                         nn_in,    & ! last significant point in y-direction for input grid  (input)
                                        mmm_out,   & !first significant point in x-direction for output grid  (input)
                                         mm_out,   & ! last significant point in x-direction for output grid  (input)
                                        nnn_out,   & !first significant point in y-direction for output grid  (input)
                                         nn_out)     ! last significant point in y-direction for output grid  (input)
   	
!  definition of data parameters
implicit none
real(8), parameter:: pi=3.1415926535897d0, pip180=pi/180.0d0
integer indper !index of periodicity:=0 -nonperiodic,=1 -periodic case

integer nxin,nyin,nxout1,nxout2,nyout1,nyout2,fillmiss

integer maskin(nxin,nyin)
real(4) maskout(nxout1:nxout2,nyout1:nyout2)

integer i_input(nxout1:nxout2,nyout1:nyout2,4),      &  !x-grid numbers of input grid for output grid
        j_input(nxout1:nxout2,nyout1:nyout2,4)          !y-grid numbers of input grid for output grid

real(8) matrx_elmnt_intrp(nxout1:nxout2,nyout1:nyout2,4)

real(8) sum_matr_elmnt

real(8) xin(nxin), yin(nyin)                	!input grid

real(8) xout(nxout1:nxout2,nyout1:nyout2),  &
        yout(nxout1:nxout2,nyout1:nyout2)       !output grid

real(8) delta_left,delta_right

integer i_nes,j_nes                             !auxilary indexes

integer i_left,i_right,j_down,j_up              !auxilary indexes

integer inxt, jnxt

integer i,j,k,m,n

real(8) bilin_denom

real(8) x_left,x_right,y_down,y_up, ret_lon, ret_lat
real(8) xl, yl, zl, xo, yo, zo, dist, dist_min
real(8) fnclon, fnclat
real(8), allocatable:: xinp(:),yinp(:),  temp(:)        !recalculated grid for periodic case                

integer mmm_in,  mm_in,  nnn_in,  nn_in, mmm_out, mm_out, nnn_out, nn_out
integer, allocatable:: i_next(:,:), j_next(:,:)
real(8), allocatable:: xd(:,:),yd(:,:),zd(:,:)

allocate (i_next(nxin,nyin), j_next(nxin,nyin))
allocate(xd(nxin,nyin),yd(nxin,nyin),zd(nxin,nyin))


do j=nnn_in, nn_in
 do i=mmm_in, mm_in
     xd(i,j)=dcos(xin(i)*pip180)*dcos(yin(j)*pip180)
     yd(i,j)=dsin(xin(i)*pip180)*dcos(yin(j)*pip180)
     zd(i,j)=dsin(yin(j)*pip180)
 enddo
enddo

do j=nnn_in, nn_in
 do i=mmm_in, mm_in

   i_next(i,j)=i
   j_next(i,j)=j
   
   if(maskin(i,j)==1) then   
     xl=xd(i,j)
     yl=yd(i,j)
     zl=zd(i,j)
     dist_min=10.0d0
     
     do n=nnn_in, nn_in
      do m=mmm_in, mm_in
       
       if(maskin(m,n)==0) then
        xo=xd(m,n)
        yo=yd(m,n)
        zo=zd(m,n)

        dist=(xo-xl)**2+(yo-yl)**2+(zo-zl)**2
        if(dist<dist_min) then
         dist_min=dist
         i_next(i,j)=m
         j_next(i,j)=n
        endif
       
       endif

      enddo
     enddo
   
   endif
 
 enddo
enddo


!      open (95,file='next.txt')
!      write(95,*) 'Number on x   ;   Number on y ;   next i  ;  next j'
!	do j=nnn_in,nn_in
!	 do i=mmm_in,mm_in
!        if(maskin(i,j)==1) then 
!         write(95,'(4i6,3f11.7,2f9.3)') i, j, i_next(i,j), j_next(i,j), xd(i,j), yd(i,j), zd(i,j), xin(i), yin(j)
!        endif
!       enddo
!      enddo
!      close(95)



!---------input grid analysis on periodicity----------------------

      if(indper/=1.and.indper/=0) then
	   write (*,*) 'wrong index of periodicity in subroutine weight_matrix_intrp'
	   stop 
	endif
	
		allocate (xinp(nxin+indper),yinp(nyin),temp(4))
! new lon and lat grid definition

	     xinp(mmm_in:mm_in)=xin(mmm_in:mm_in)
            if(indper==1) xinp(mm_in+1)=xin(mmm_in)+360.0d0 !completing of longitude
           yinp=yin

! interpolation matrix construction

	do j=nnn_out,nn_out
	 do i=mmm_out,mm_out
           
	  if (maskout(i,j)>0.5) then
            ret_lon=xout(i,j)
            ret_lat=yout(i,j)

!          corrections

            if(indper==1) then
	        if(ret_lon > xinp(mm_in+1)) ret_lon=ret_lon-360.0d0
	        if(ret_lon < xinp(mmm_in) ) ret_lon=ret_lon+360.0d0
            else
              
              if(ret_lon > xinp(mm_in)) then
	         delta_right=ret_lon-xinp(mm_in)
	         ret_lon=ret_lon-360.0d0
               delta_left=xin(mmm_in)-ret_lon
	         if(delta_left > delta_right) ret_lon=ret_lon+360.0d0
	         go to 111
	        end if

              if(ret_lon < xinp(mmm_in)) then
	         delta_left=xinp(mmm_in)-ret_lon
	         ret_lon=ret_lon+360.0d0
 	         delta_right=ret_lon-xinp(mm_in)
	         if(delta_right > delta_left) ret_lon=ret_lon-360.0d0
	         go to 111
	        end if
               	       
            end if
111	 continue

!       identification of reversly rotated point on geografic grid  
       
        i_nes=mmm_in
	 
	   n=mm_in+indper

1000     m=(i_nes+n)/2	 
	     
           if (ret_lon < xinp(m)) n=m
	     if (ret_lon >=xinp(m)) i_nes=m
	     if (n > i_nes+1) goto 1000
	 
         i_left =i_nes
	   i_right=i_nes+1
         if((mm_in-mmm_in+1)==1) i_right=i_nes

	   j_nes=nnn_in
	   n=nn_in
1001     m=(j_nes+n)/2	 
	     
           if (ret_lat < yinp(m)) n=m
	     if (ret_lat >=yinp(m)) j_nes=m
	     if (n > j_nes+1) goto 1001
        
	   j_down = j_nes
         j_up   = j_nes+1
         if((nn_in-nnn_in+1)==1) j_up=j_nes

!       main bilinear interpolation

!--------------definition of output i-indexes and preparing of matrix calculating----

1982     if((i_left < mmm_in)   .and.(indper==1)) then
          
          i_input(i,j,1)= i_left + (mm_in-mmm_in+1)
	    i_input(i,j,3)= i_left + (mm_in-mmm_in+1)

	    x_left = xinp(i_input(i,j,1))-360.

         else
          
	    i_input(i,j,1)=i_left
          i_input(i,j,3)=i_left

          x_left =xinp(i_input(i,j,1))

	   end if
	            
	   
	   if((i_right > mm_in).and.(indper==1)) then
          
		i_input(i,j,2)=i_right-(mm_in-mmm_in+1)
		i_input(i,j,4)=i_right-(mm_in-mmm_in+1)

		x_right=xinp(i_input(i,j,2))+360.
	   
	   else
          
		i_input(i,j,2)=i_right
		i_input(i,j,4)=i_right

          x_right=xinp(i_input(i,j,2))
	   
	   end if

!--------------definition of output j-indexes--------------
	   j_input(i,j,1)=j_down
	   j_input(i,j,2)=j_down
	   j_input(i,j,3)=j_up
	   j_input(i,j,4)=j_up

	   y_down =yinp(j_input(i,j,1))
	   y_up   =yinp(j_input(i,j,3))

!---------if all 4 points are undefined checking fillmiss condition	   
	   if((maskin(i_input(i,j,1),j_input(i,j,1))*        &
             maskin(i_input(i,j,2),j_input(i,j,2))*        &
             maskin(i_input(i,j,3),j_input(i,j,3))*        &
             maskin(i_input(i,j,4),j_input(i,j,4)))/=0) then

           if(fillmiss==0) then
            matrx_elmnt_intrp(i,j,1:4)=0.25d0
            go to 501
	     end if
         
         endif

	   bilin_denom= 1.0d0

	   if (x_right /= x_left) bilin_denom=bilin_denom / (fnclon(x_right)-fnclon(x_left))
	   if (y_up    /= y_down) bilin_denom=bilin_denom / (fnclat(y_up)   -fnclat(y_down))
	 	 
	   temp=bilin_denom


	    if (x_right /=x_left) then

           temp(1)= temp(1)*(fnclon(x_right) - fnclon(ret_lon))
           temp(2)= temp(2)*(fnclon(ret_lon) - fnclon(x_left) )
	     temp(3)= temp(3)*(fnclon(x_right) - fnclon(ret_lon)) 
	     temp(4)= temp(4)*(fnclon(ret_lon) - fnclon(x_left) )

          end if

	    if (y_up    /=y_down) then

           temp(1)= temp(1)*(fnclat(y_up)    - fnclat(ret_lat))
           temp(2)= temp(2)*(fnclat(y_up)    - fnclat(ret_lat))
	     temp(3)= temp(3)*(fnclat(ret_lat) - fnclat(y_down) ) 
	     temp(4)= temp(4)*(fnclat(ret_lat) - fnclat(y_down) )

          end if
          
          do k=1,4
           matrx_elmnt_intrp(i,j,k)=temp(k)
           inxt=i_next(i_input(i,j,k),j_input(i,j,k))
           jnxt=j_next(i_input(i,j,k),j_input(i,j,k))
           i_input(i,j,k)=inxt
           j_input(i,j,k)=jnxt
          enddo
	  else
          matrx_elmnt_intrp(i,j,:)=0.0d0
                    i_input(i,j,:)=0
                    j_input(i,j,:)=0
        end if
      
      matrx_elmnt_intrp(I,J,:)=min(matrx_elmnt_intrp(i,j,:),1.0d0)
      matrx_elmnt_intrp(I,J,:)=max(matrx_elmnt_intrp(i,j,:),0.0d0)  
!--------remorming of matrix elements-------------------------------
500       sum_matr_elmnt=matrx_elmnt_intrp(i,j,1)  &
                        +matrx_elmnt_intrp(i,j,2)  &
                        +matrx_elmnt_intrp(i,j,3)  &
                        +matrx_elmnt_intrp(i,j,4)

	  matrx_elmnt_intrp(i,j,:)= matrx_elmnt_intrp(i,j,:)/sum_matr_elmnt

501    continue     	  	  
	 end do
	end do


deallocate (temp,yinp,xinp)
deallocate(zd,yd,xd)
deallocate(j_next,i_next)

endsubroutine weight_matrix_intrp_next

!================================================
function fnclat(latitude)
implicit none
!  interpolation weight functions for grid coordinates
      real(8) latitude,fnclat,lat
      real(8), parameter:: dpip180=3.1415926535897/180.0d0, lat_extr=89.999999d0
!   bilinear interpolation
!     fnclat(latitude) =latitude
      lat=max(min(latitude,lat_extr),-lat_extr)
!   harmonic interpolation
      fnclat=dlog(   (1d0+dsin(dpip180*lat))  &
                   / (1d0-dsin(dpip180*lat))   )
endfunction fnclat

!================================================
function fnclon(longitude)
implicit none
!  interpolation weight functions for grid coordinates
      real(8) longitude,fnclon
!   bilinear interpolation
      fnclon =longitude
endfunction fnclon
