!========================================================================
subroutine sea_surface_fluxes_simple
    use main_basin_pars
    use mpi_parallel_tools
    use basin_grid
    use ocean_variables
    use atm_forcing
    implicit none

    integer :: m, n
    real*8 :: wnd, wnd_mod
    real*8 :: coeff_surf_fric

    coeff_surf_fric = 3.25d-6

    !$omp parallel do private(m,n,wnd,wnd_mod)
    do n=ny_start,ny_end
        do m=nx_start,nx_end

        if(lcu(m,n)>0.5) then
            wnd = (uwnd(m  ,n)*dx(m  ,n)*dy(m  ,n)                              &
                    + uwnd(m+1,n)*dx(m+1,n)*dy(m+1,n))/2.0d0/dxt(m,n)/dyh(m,n)
            wnd_mod = (wind(m  ,n)*dx(m  ,n)*dy(m  ,n)                          &
                        + wind(m+1,n)*dx(m+1,n)*dy(m+1,n))/2.0d0/dxt(m,n)/dyh(m,n)

            surf_stress_x(m,n) = wnd * wnd_mod * coeff_surf_fric
        endif

        if(lcv(m,n)>0.5) then
            wnd = (vwnd(m,n  )*dx(m,n  )*dy(m,n  )                              &
                    + vwnd(m,n+1)*dx(m,n+1)*dy(m,n+1))/2.0d0/dxh(m,n)/dyt(m,n)
            wnd_mod = (wind(m,n  )*dx(m,n  )*dy(m,n  )                          &
                        + wind(m,n+1)*dx(m,n+1)*dy(m,n+1))/2.0d0/dxh(m,n)/dyt(m,n)

            surf_stress_y(m,n) = wnd * wnd_mod * coeff_surf_fric
        endif

        enddo
    enddo
    !$omp end parallel do

endsubroutine sea_surface_fluxes_simple

!========================================================================
subroutine sea_surface_fluxes
    use main_basin_pars
    use mpi_parallel_tools
    use basin_grid
    use ocean_variables
    use atm_forcing
    implicit none

    integer m,n,k,ierr

    real(8) evap_rate
    real(8) wf_ave, sf_ave, m_calc, wf, wnd
    real(8) tmp_real8

    if(ksw_ssbc==1 .or. ksw_ssbc==2) then
        !$omp parallel do
        do n=ny_start,ny_end
            do m=nx_start,nx_end

            if(lcu(m,n)>0.5) then
            surf_stress_x(m,n)= (taux(m  ,n)*dx(m  ,n)*dy(m  ,n)   &
                                +taux(m+1,n)*dx(m+1,n)*dy(m+1,n))/RefDen/2.0d0/dxt(m,n)/dyh(m,n)
            endif

            if(lcv(m,n)>0.5) then
            surf_stress_y(m,n)= (tauy(m,n  )*dx(m,n  )*dy(m,n  )   &
                                +tauy(m,n+1)*dx(m,n+1)*dy(m,n+1))/RefDen/2.0d0/dxh(m,n)/dyt(m,n)
            endif

            enddo
        enddo
        !$omp end parallel do
    endif

    if(ksw_ssbc==3) then
        !$omp parallel do private(m, n, evap_rate)
        do n=ny_start,ny_end
            do m=nx_start,nx_end
                if(lu(m,n)>0.5) then

                    call air_sea_turbulent_fluxes(  wind(m,n),         &   ! wind modulo, m/s
                                                    slpr(m,n),         &   ! sea level pressure, Pa
                                                    tatm(m,n),         &   ! air temperature,  �C
                                                    tatm(m,n),         &   ! sea surface temp, �C
                                                    qatm(m,n),         &   ! air specific humidity, kg/kg
                                                     u_height,         &   ! height of wind datasets, m
                                                     t_height,         &   ! height of tair datasets, m
                                                     q_height,         &   ! height of qair datasets, m
                                                    taux(m,n),         &   !      zonal wind stress, Pa
                                                    tauy(m,n)          )   ! meridional wind stress, Pa

                    wf_tot(m,n)=rain(m,n)+snow(m,n)

                endif
            enddo
        enddo
        !$omp end parallel do

        call syncborder_real8(taux, 1)
        call syncborder_real8(tauy, 1)
        if(periodicity_x/=0) then
            call cyclize8_x(taux,nx,ny,1,mmm,mm)
        endif

        if(periodicity_y/=0) then
            call cyclize8_y(tauy,nx,ny,1,nnn,nn)
        endif

        wf_ave=0.0d0
        sf_ave=0.0d0
        m_calc=0.0d0
        if(ksw_wflux>0) then
            !$omp parallel do private(m,n) reduction(+:wf_ave, m_calc)
            do n=ny_start,ny_end
                do m=nx_start,nx_end
                    if(lu(m,n)>0.5) then
                        wf_ave=wf_ave+wf_tot(m,n)*dx(m,n)*dy(m,n)
                        m_calc=m_calc+dx(m,n)*dy(m,n)
                    endif
                enddo
            enddo
            !$omp end parallel do
            tmp_real8 = wf_ave
            call mpi_allreduce(tmp_real8, wf_ave, 1, mpi_real8, mpi_sum, cart_comm, ierr)
            tmp_real8 = m_calc
            call mpi_allreduce(tmp_real8, m_calc, 1, mpi_real8, mpi_sum, cart_comm, ierr)

            wf_tot = wf_tot - wf_ave/m_calc
        endif

        call syncborder_real8(wf_tot, 1)
        if(periodicity_x/=0) then
            call cyclize8_x(wf_tot,nx,ny,1,mmm,mm)
        endif
        if(periodicity_y/=0) then
            call cyclize8_y(wf_tot,nx,ny,1,nnn,nn)
        endif

        !$omp parallel do private(m,n,wf,wnd)
        do n=ny_start,ny_end
            do m=nx_start,nx_end
                if(lcu(m,n)>0.5) then
                    wnd = (uwnd(m  ,n)*dx(m  ,n)*dy(m  ,n)         &
                        + uwnd(m+1,n)*dx(m+1,n)*dy(m+1,n))/2.0d0/dxt(m,n)/dyh(m,n)

                    wf = (wf_tot(m  ,n)*uwnd(m  ,n)*dx(m  ,n)*dy(m  ,n)         &
                        + wf_tot(m+1,n)*uwnd(m+1,n)*dx(m+1,n)*dy(m+1,n))/2.0d0/dxt(m,n)/dyh(m,n)

                    surf_stress_x(m,n) = (taux(m,n) + taux(m+1,n))/RefDen/2.0d0     &
                        *( wnd ) + wf/RefDen
                endif

                if(lcv(m,n)>0.5) then
                    wnd = (vwnd(m,n  )*dx(m,n  )*dy(m,n  )         &
                    + vwnd(m,n+1)*dx(m,n+1)*dy(m,n+1))/2.0d0/dxh(m,n)/dyt(m,n)

                    wf = (wf_tot(m,n  )*vwnd(m,n  )*dx(m,n  )*dy(m,n  )         &
                    +wf_tot(m,n+1)*vwnd(m,n+1)*dx(m,n+1)*dy(m,n+1))/2.0d0/dxh(m,n)/dyt(m,n)

                    surf_stress_y(m,n) = (tauy(m,n) + tauy(m,n+1))/RefDen/2.0d0     &
                           *( wnd ) + wf/RefDen
                endif

            enddo
        enddo
        !$omp end parallel do

    endif

endsubroutine sea_surface_fluxes
