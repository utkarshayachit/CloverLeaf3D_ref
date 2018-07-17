!Crown Copyright 2012 AWE.
!
! This file is part of CloverLeaf.
!
! CloverLeaf is free software: you can redistribute it and/or modify it under 
! the terms of the GNU General Public License as published by the 
! Free Software Foundation, either version 3 of the License, or (at your option) 
! any later version.
!
! CloverLeaf is distributed in the hope that it will be useful, but 
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
! FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
! details.
!
! You should have received a copy of the GNU General Public License along with 
! CloverLeaf. If not, see http://www.gnu.org/licenses/.

!>  @brief SENSEI support
!>  @author Utkarsh Ayachit
!>  @details Invokes SENSEI in situ library

SUBROUTINE sensei

  USE clover_module
  USE update_halo_module
  USE viscosity_module
  USE ideal_gas_module

  IMPLICIT NONE

  INTEGER :: j,k,l,tile,err,get_unit,u,dummy
  INTEGER :: nxc,nyc,nzc,nxv,nyv,nzv,nblocks,xmin,ymin,zmin
  REAL(KIND=8)    :: temp_var

  CHARACTER(len=10)           :: chunk_name,step_name
  CHARACTER(len=90)           :: filename

  LOGICAL, SAVE :: first_call=.TRUE.

  INTEGER :: fields(NUM_FIELDS)

  REAL(KIND=8) :: kernel_time,timer


  IF(first_call) THEN
    CALL cloverleaf3d_sensei_bridge_init(parallel%max_task, tiles_per_chunk)
    first_call=.FALSE.
  ENDIF

  IF(profiler_on) kernel_time=timer()

  DO tile=1,tiles_per_chunk
    CALL ideal_gas(tile,.FALSE.)
  ENDDO

  IF(profiler_on) profiler%ideal_gas=profiler%ideal_gas+(timer()-kernel_time)

  fields=0
  fields(FIELD_PRESSURE)=1
  fields(FIELD_XVEL0)=1
  fields(FIELD_YVEL0)=1
  fields(FIELD_ZVEL0)=1

  CALL update_halo(fields,1)


  IF(profiler_on) kernel_time=timer()
  CALL viscosity()
  IF(profiler_on) profiler%viscosity=profiler%viscosity+(timer()-kernel_time)

  IF(profiler_on) kernel_time=timer()
  DO tile = 1, tiles_per_chunk
    IF(chunk%task.EQ.parallel%task) THEN
      nxc=chunk%tiles(tile)%t_xmax-chunk%tiles(tile)%t_xmin+1
      nyc=chunk%tiles(tile)%t_ymax-chunk%tiles(tile)%t_ymin+1
      nzc=chunk%tiles(tile)%t_zmax-chunk%tiles(tile)%t_zmin+1
      nxv=nxc+1
      nyv=nyc+1
      nzv=nzc+1

      CALL cloverleaf3d_sensei_bridge_update(                       &
        chunk%task, tile, time, step,                               &
        chunk%tiles(tile)%t_xmin, chunk%tiles(tile)%t_xmax,         &
        chunk%tiles(tile)%t_ymin, chunk%tiles(tile)%t_ymax,         &
        chunk%tiles(tile)%t_zmin, chunk%tiles(tile)%t_zmax,         &
        chunk%tiles(tile)%field%vertexx,                            &
        chunk%tiles(tile)%field%vertexy,                            &
        chunk%tiles(tile)%field%vertexz,                            &
        chunk%tiles(tile)%field%density0,                           &
        chunk%tiles(tile)%field%energy0,                            &
        chunk%tiles(tile)%field%pressure,                           &
        chunk%tiles(tile)%field%viscosity)
    ENDIF
  ENDDO
  IF(profiler_on) profiler%visit=profiler%visit+(timer()-kernel_time)

END SUBROUTINE sensei
