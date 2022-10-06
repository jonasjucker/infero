program ecrad_ml
  use inferof
  use iso_c_binding, only : c_float, c_null_char
  implicit none

  ! Command line arguments
  character(1024) :: model_path
  character(1024) :: model_type
  character(1024) :: yaml_config

  ! model of infero model
  type(infero_model) :: model
  type(infero_tensor_set) :: iset
  type(infero_tensor_set) :: oset

  ! parameters
  REAL (c_float), PARAMETER ::  pi = 3.14159265358979323846264338327950288
  REAL (c_float), PARAMETER ::  rad2deg   = 180.0/pi

  ! define gridsize (icon grid indices)
  !integer,parameter :: batch_size = 10
  !integer,parameter :: batch_size = 100
  !integer,parameter :: batch_size = 1000
  integer,parameter :: batch_size = 10000
  !integer,parameter :: batch_size = 81919

  ! input and output tensors
  real(c_float) :: input_3d(batch_size,60,1,6)
  real(c_float) :: input_2d(batch_size,1,8)
  real(c_float) ::  pred_flx(batch_size,60,4)

  ! netcdf
  character(1024) :: netcdf_data_file,varname,icon_grid
  character(50) :: dim_name(6),grid_dim_name(14)
  character(19) :: timestamp(4)
  integer :: dim_len(6),grid_dim_len(14)

  ! indices
  integer :: step, nsteps, s_idx,e_idx, nc_time_idx(4),i,k,id_d, counter(4)

  ! data fields
  real(c_float), allocatable :: from_netcdf_3d(:,:,:,:), from_netcdf_2d(:,:,:), neighbor_cell_index(:,:), abs_diff(:,:,:,:)
  real(c_float), allocatable :: swflx(:,:,:,:), lwflx(:,:,:,:)
  real(c_float), allocatable :: clat(:), clon(:),lon(:),lat(:)
  real(c_float):: mean_absolute_error(4,4)

  ! scalars
  real(c_float) :: percentage, lat_min,lat_max,lon_min,lon_max

  ! Get CL arguments
  CALL get_command_argument(1, model_path)
  CALL get_command_argument(2, model_type)

  ! SETUP INFERO AND ASSIGN IN/OUT TENSORS

  ! init infero
  call infero_check(infero_initialise())

  ! YAML config string -> LW/SW
  yaml_config = "---"//NEW_LINE('A') &
    //"  path: "//TRIM(model_path)//NEW_LINE('A') &
    //"  type: "//TRIM(model_type)//c_null_char

  ! get a infero model for LW/SW
  call infero_check(model%initialise_from_yaml_string(yaml_config))


  ! init tensor sets
  call infero_check(iset%initialise())
  call infero_check(oset%initialise())

  ! input is identical for both models
  call infero_check(iset%push_tensor(input_3d, "serving_default_input_22"))
  call infero_check(iset%push_tensor(input_2d, "serving_default_input_23"))

  call infero_check(oset%push_tensor(pred_flx, "StatefulPartitionedCall"))



  ! DEFINE PERIOD (FULL SOLAR CYCLE)
  nsteps = 4

  ! correspond to file nr. 5
  timestamp(1) = "2000-01-29 00:00:00"
  timestamp(2) = "2000-01-29 06:00:00"
  timestamp(3) = "2000-01-29 12:00:00"
  timestamp(4) = "2000-01-29 18:00:00"

  ! store NetCDF indices of timestamps defined above
  nc_time_idx(1) = 1
  nc_time_idx(2) = 2
  nc_time_idx(3) = 3
  nc_time_idx(4) = 4

  ! fields to store model output
  ALLOCATE(swflx(batch_size,60,2,nsteps))
  ALLOCATE(lwflx(batch_size,60,2,nsteps))

  ! READ-IN ICON GRID INFORMATION AND COMPUTE DOMAIN EXTENT
  icon_grid='icon_grid.nc'

  ! read icon grid file dimensions
  call get_nc_dims(icon_grid,grid_dim_name,grid_dim_len,14)

  ! fields with icon-grid information
  ALLOCATE(neighbor_cell_index(grid_dim_len(1),grid_dim_len(5)))
  ALLOCATE(clat(grid_dim_len(1)))
  ALLOCATE(clon(grid_dim_len(1)))
  ALLOCATE(lat(batch_size + 1))
  ALLOCATE(lon(batch_size + 1))

  varname="neighbor_cell_index"
  call read_nc_2d(icon_grid,varname,neighbor_cell_index(:,:),grid_dim_len(1),grid_dim_len(5))

  varname="clat"
  call read_nc_1d(icon_grid,varname,clat,grid_dim_len(1))
  varname="clon"
  call read_nc_1d(icon_grid,varname,clon,grid_dim_len(1))

  ! compute domain extent
  DO i=1,batch_size+1
    lon(i) = rad2deg * clon(i)
    lat(i) = rad2deg * clat(i)
  ENDDO

  lon_max = MAXVAL(lon(:))
  lon_min = MINVAL(lon(:))
  lat_max = MAXVAL(lat(:))
  lat_min = MINVAL(lat(:))


  ! READ-IN INPUT-DATA FOR MODEL
  netcdf_data_file='input_data.nc'

  ! read data dimensions
  call get_nc_dims(netcdf_data_file,dim_name,dim_len,6)

  ! fields for NetCDF data
  ALLOCATE(from_netcdf_3d(dim_len(1),dim_len(3),dim_len(6),10))
  ALLOCATE(from_netcdf_2d(dim_len(1),dim_len(6),8))

  ! 2d fields
  varname="pres_sfc"
  call read_nc_2d(netcdf_data_file,varname,from_netcdf_2d(:,:,1),dim_len(1),dim_len(6))

  varname="cosmu0"
  call read_nc_2d(netcdf_data_file,varname,from_netcdf_2d(:,:,2),dim_len(1),dim_len(6))

  varname="qv_s"
  call read_nc_2d(netcdf_data_file,varname,from_netcdf_2d(:,:,3),dim_len(1),dim_len(6))

  varname="albvisdir"
  call read_nc_2d(netcdf_data_file,varname,from_netcdf_2d(:,:,4),dim_len(1),dim_len(6))

  varname="albnirdir"
  call read_nc_2d(netcdf_data_file,varname,from_netcdf_2d(:,:,5),dim_len(1),dim_len(6))

  varname="tsfctrad"
  call read_nc_2d(netcdf_data_file,varname,from_netcdf_2d(:,:,6),dim_len(1),dim_len(6))

  varname="albvisdif"
  call read_nc_2d(netcdf_data_file,varname,from_netcdf_2d(:,:,7),dim_len(1),dim_len(6))

  varname="albnirdif"
  call read_nc_2d(netcdf_data_file,varname,from_netcdf_2d(:,:,8),dim_len(1),dim_len(6))

  ! 3d fields
  varname="clc"
  call read_nc_3d(netcdf_data_file,varname,from_netcdf_3d(:,:,:,1),dim_len(1),dim_len(3),dim_len(6))

  varname="temp"
  call read_nc_3d(netcdf_data_file,varname,from_netcdf_3d(:,:,:,2),dim_len(1),dim_len(3),dim_len(6))

  varname="pres"
  call read_nc_3d(netcdf_data_file,varname,from_netcdf_3d(:,:,:,3),dim_len(1),dim_len(3),dim_len(6))

  varname="qc"
  call read_nc_3d(netcdf_data_file,varname,from_netcdf_3d(:,:,:,4),dim_len(1),dim_len(3),dim_len(6))

  varname="qi"
  call read_nc_3d(netcdf_data_file,varname,from_netcdf_3d(:,:,:,5),dim_len(1),dim_len(3),dim_len(6))

  varname="qv"
  call read_nc_3d(netcdf_data_file,varname,from_netcdf_3d(:,:,:,6),dim_len(1),dim_len(3),dim_len(6))

  varname="lwflx_up"
  call read_nc_3d(netcdf_data_file,varname,from_netcdf_3d(:,:,:,7),dim_len(1),dim_len(3),dim_len(6))

  varname="lwflx_dn"
  call read_nc_3d(netcdf_data_file,varname,from_netcdf_3d(:,:,:,8),dim_len(1),dim_len(3),dim_len(6))

  varname="swflx_up"
  call read_nc_3d(netcdf_data_file,varname,from_netcdf_3d(:,:,:,9),dim_len(1),dim_len(3),dim_len(6))

  varname="swflx_dn"
  call read_nc_3d(netcdf_data_file,varname,from_netcdf_3d(:,:,:,10),dim_len(1),dim_len(3),dim_len(6))


  ! TIMESTEP
  s_idx = 1
  e_idx = 1 + batch_size

  DO step=1,nsteps

    ! update 2D input-tensor

    ! center cell
    input_2d(:,1,:) = from_netcdf_2d(s_idx:e_idx,nc_time_idx(step),:)

    ! update 3D input-tensor

    ! center cell
    input_3d(:,:,1,:) = from_netcdf_3d(s_idx:e_idx,:,nc_time_idx(step),:)


    ! apply model
    call infero_check(model%infer(iset, oset ))

    !CALL infero_check(iset%print())
    !CALL infero_check(oset%print())


    ! store results in permanent fields
    lwflx(:,:,:,step) = pred_flx(:,:,1:2)
    swflx(:,:,:,step) = pred_flx(:,:,3:4)

    ! input
    varname="input_clc"
    CALL stats_2d(varname,input_3d(:,:,1,1))
    varname="input_temp"
    CALL stats_2d(varname,input_3d(:,:,1,2))
    varname="input_pres"
    CALL stats_2d(varname,input_3d(:,:,1,3))
    varname="input_qc"
    CALL stats_2d(varname,input_3d(:,:,1,4))
    varname="input_qi"
    CALL stats_2d('input_qi',input_3d(:,:,1,5))
    varname="input_qv"
    CALL stats_2d(varname,input_3d(:,:,1,6))

    ! output
    varname="lwflx_up"
    CALL stats_2d(varname,pred_flx(:,:,1))
    varname="lwflx_dn"
    CALL stats_2d(varname,pred_flx(:,:,2))
    varname="swflx_up"
    CALL stats_2d(varname,pred_flx(:,:,3))
    varname="swflx_dn"
    CALL stats_2d(varname,pred_flx(:,:,4))
  ENDDO


  ! DOMAIN EXTENT
  write(*,'(A,I7,A)') 'Domain extent for batch_size', batch_size,':'
  write(*,'(A,F8.2,A,F8.2)') '  latmax:', lat_max, ' latmin', lat_min
  write(*,'(A,F8.2,A,F8.2)') '  lonmax:', lon_max, ' lonmin', lon_min


  ! STATISTICS
  write(*,'(A)') ''
  write(*,'(A)') 'STATISTICS'
  write(*,'(A)') ''

  ! values below zero for SW
  counter(:) = 0
  DO step=1,nsteps
    DO i=1,batch_size
      DO k=1,60
        DO id_d=1,2
          IF(swflx(i,k,id_d,step) < 0.0) THEN
            counter(step) = counter(step) + 1
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  write(*,'(I7,A)') batch_size*60*2,' datapoints scanned for each step of SW-flux'
  DO step=1,nsteps
    percentage = MAX(0.0,100.0 * REAL(counter(step))/REAL((batch_size*60*2)))
    write(*,'(F6.1,A,A)') percentage, '% values below 0.0 for ',timestamp(step)
  ENDDO

  ! absolute difference 
  ALLOCATE(abs_diff(batch_size+1,dim_len(3),dim_len(6),4))

  write(*,'(A)') ''
  write(*,'(A)') '  Mean Absolute Error (MAE):'
  DO step=1,nsteps
    ! lwflux_up
    abs_diff(:,:,1,step) = ABS(lwflx(:,:,1,step) - from_netcdf_3d(s_idx:e_idx,:,nc_time_idx(step),7))
    ! lwflux_dn
    abs_diff(:,:,2,step) = ABS(lwflx(:,:,2,step) - from_netcdf_3d(s_idx:e_idx,:,nc_time_idx(step),8))
    ! swflux_up
    abs_diff(:,:,3,step) = ABS(swflx(:,:,1,step) - from_netcdf_3d(s_idx:e_idx,:,nc_time_idx(step),9))
    ! swflux_dn
    abs_diff(:,:,4,step) = ABS(swflx(:,:,2,step) - from_netcdf_3d(s_idx:e_idx,:,nc_time_idx(step),10))

    ! mean values
    CALL mean_2d(abs_diff(:,:,1,step), mean_absolute_error(1,step),batch_size,60)
    CALL mean_2d(abs_diff(:,:,2,step), mean_absolute_error(2,step),batch_size,60)
    CALL mean_2d(abs_diff(:,:,3,step), mean_absolute_error(3,step),batch_size,60)
    CALL mean_2d(abs_diff(:,:,4,step), mean_absolute_error(4,step),batch_size,60)

    write(*,'(A)') ''
    write(*,'(A,A)') '    ',timestamp(step)
    write(*,'(A)') '      Upwards:'
    write(*,'(A,F8.2)') '        SW (all levels): ',mean_absolute_error(3,step)
    write(*,'(A,F8.2)') '        LW (all levels): ',mean_absolute_error(1,step)
    write(*,'(A)') '      Downwards:'
    write(*,'(A,F8.2)') '        SW (all levels): ',mean_absolute_error(4,step)
    write(*,'(A,F8.2)') '        LW (all levels): ',mean_absolute_error(2,step)
  ENDDO

  write(*,'(A)') ''
  write(*,'(A)') '  MIN/MAX values for different levels'
  DO step=1,nsteps
    write(*,'(A)') ''
    write(*,'(A,A)') '  ',timestamp(step)
    write(*,'(A)') '    Upwards:'
    write(*,'(A,F8.2,A,F8.2)') '     SW (all levels): MAX ',MAXVAL(swflx(:,:,1,step)),' MIN ',MINVAL(swflx(:,:,1,step))
    write(*,'(A,F8.2,A,F8.2)') '     LW (all levels): MAX ',MAXVAL(lwflx(:,:,1,step)),' MIN ', MINVAL(lwflx(:,:,1,step))
    write(*,'(A,F8.2,A,F8.2)') '     SW (sfc): MAX ',MAXVAL(swflx(:,1,1,step)),' MIN ', MINVAL(swflx(:,1,1,step))
    write(*,'(A,F8.2,A,F8.2)') '     LW (sfc): MAX ',MAXVAL(lwflx(:,1,1,step)),' MIN ', MINVAL(lwflx(:,1,1,step))
    write(*,'(A,F8.2,A,F8.2)') '     SW (toa): MAX ',MAXVAL(swflx(:,60,1,step)),' MIN ', MINVAL(swflx(:,60,1,step))
    write(*,'(A,F8.2,A,F8.2)') '     LW (toa): MAX ',MAXVAL(lwflx(:,60,1,step)),' MIN ', MINVAL(lwflx(:,60,1,step))
    write(*,'(A)') ''
    write(*,'(A)') '    Downwards:'
    write(*,'(A,F8.2,A,F8.2)') '     SW (all levels): MAX ',MAXVAL(swflx(:,:,2,step)),' MIN ', MINVAL(swflx(:,:,2,step))
    write(*,'(A,F8.2,A,F8.2)') '     LW (all levels): MAX ',MAXVAL(lwflx(:,:,2,step)),' MIN ', MINVAL(lwflx(:,:,2,step))
    write(*,'(A,F8.2,A,F8.2)') '     SW (sfc): MAX ',MAXVAL(swflx(:,1,2,step)),' MIN ', MINVAL(swflx(:,1,2,step))
    write(*,'(A,F8.2,A,F8.2)') '     LW (sfc): MAX ',MAXVAL(lwflx(:,1,2,step)),' MIN ', MINVAL(lwflx(:,1,2,step))
    write(*,'(A,F8.2,A,F8.2)') '     SW (toa): MAX ',MAXVAL(swflx(:,60,2,step)),' MIN ', MINVAL(swflx(:,60,2,step))
    write(*,'(A,F8.2,A,F8.2)') '     LW (toa): MAX ',MAXVAL(lwflx(:,60,2,step)),' MIN ', MINVAL(lwflx(:,60,2,step))
  ENDDO


  ! CLEANUP 

  ! free the model
  call infero_check(model%free())

  ! finalise infero library
  call infero_check(infero_finalise())

  ! deallocate
  DEALLOCATE(swflx)
  DEALLOCATE(lwflx)
  DEALLOCATE(neighbor_cell_index)
  DEALLOCATE(clat)
  DEALLOCATE(clon)
  DEALLOCATE(lat)
  DEALLOCATE(lon)
  DEALLOCATE(from_netcdf_3d)
  DEALLOCATE(from_netcdf_2d)
  DEALLOCATE(abs_diff)

CONTAINS

  SUBROUTINE stats_2d(varname,input)
    use iso_c_binding, only : c_float
    REAL(c_float), INTENT(IN) :: input(:,:)
    CHARACTER(*) :: varname

    INTEGER :: i,j
    REAL(c_float):: mean

    mean = 0.0
    DO i=1,SIZE(input,DIM=1)
      DO j=1,SIZE(input,DIM=2)
        mean = mean + input(i,j)
      ENDDO
    ENDDO

    mean = mean / REAL(i*j)

    !JJ: uncomment for shape information
    !write(message_text,*) SHAPE(input)
    !write(message_text,'(A)') TRIM(message_text)
    !CALL message(TRIM(varname),message_text)

    write(*,'(A,A,F8.2,A,F8.2,A,F8.2)') TRIM(varname),': MEAN ',mean,' MAX ',MAXVAL(input(:,:)),' MIN ', MINVAL(input(:,:))

  END SUBROUTINE stats_2d

end program

SUBROUTINE get_nc_dims(infile,dim_name,dim_len,nr_dims)
  use netcdf

  INTEGER(KIND=4), INTENT(OUT) :: dim_len(nr_dims)
  CHARACTER(LEN=50), INTENT(OUT) :: dim_name(nr_dims)
  INTEGER(KIND=4) :: ncid,i
  CHARACTER(LEN=1024), INTENT(IN) :: infile
  !Open netCDF file
  !:-------:-------:-------:-------:-------:-------:-------:
  CALL check(nf90_open(TRIM(infile), nf90_nowrite, ncid))
  !Inquire about the dimensions
  !:-------:-------:-------:-------:-------:-------:-------:

  write(*,'(A,A)') 'Dimension of NetCDF: ',TRIM(infile)
  DO i=1,nr_dims
    CALL check(nf90_inquire_dimension(ncid,i,dim_name(i),dim_len(i)))
    write(*,'(A,A,A,I7)') '  ',TRIM(dim_name(i)), ':', dim_len(i)
  ENDDO
  !Close netCDF file
  !:-------:-------:-------:-------:-------:-------:-------:
  CALL check(nf90_close(ncid))
END SUBROUTINE get_nc_dims

SUBROUTINE read_nc_3d(infile,varname,idata,nx,ny,nz)
  use iso_c_binding, only : c_float
  use netcdf

  REAL(c_float), DIMENSION(nx,ny,nz), INTENT(OUT) :: idata
  INTEGER(KIND=4), INTENT(IN) :: nx,ny,nz
  INTEGER(KIND=4), DIMENSION(3) :: dimids
  INTEGER(KIND=4) :: ncid, ndims, varid
  CHARACTER(LEN=1024), INTENT(IN) :: infile,varname
  !Open netCDF file
  !:-------:-------:-------:-------:-------:-------:-------:-------:
  CALL check(nf90_open(TRIM(infile), nf90_nowrite, ncid))

  !Get the values of the coordinates and put them in xpos & ypos
  !:-------:-------:-------:-------:-------:-------:-------:-------:
  !CALL check(nf90_inquire_variable(ncid,1,vname,xtype,ndims,dimids))
  CALL check(nf90_inq_varid(ncid,TRIM(varname),varid))
  CALL check(nf90_inquire_variable(ncid=ncid,varid=varid,ndims=ndims))
  CALL check(nf90_inquire_variable(ncid=ncid,varid=varid,ndims=ndims,dimids=dimids))

  write(*,'(A)') TRIM(varname)
  write(*,'(A,I2)') '   ndims:',ndims
  write(*,*) '   dimids:',dimids

  CALL check(nf90_get_var(ncid,varid,idata))
  !Close netCDF file
  !:-------:-------:-------:-------:-------:-------:-------:-------:
  CALL check(nf90_close(ncid))

END SUBROUTINE read_nc_3d

SUBROUTINE read_nc_1d(infile,varname,idata,nx)
  use iso_c_binding, only : c_float
  use netcdf

  REAL(c_float), DIMENSION(nx), INTENT(OUT) :: idata
  INTEGER(KIND=4), INTENT(IN) :: nx
  INTEGER(KIND=4), DIMENSION(1) :: dimids
  INTEGER(KIND=4) :: ncid, ndims, varid
  CHARACTER(LEN=1024), INTENT(IN) :: infile,varname
  !Open netCDF file
  !:-------:-------:-------:-------:-------:-------:-------:-------:
  CALL check(nf90_open(TRIM(infile), nf90_nowrite, ncid))
  CALL check(nf90_inq_varid(ncid,TRIM(varname),varid))
  CALL check(nf90_inquire_variable(ncid=ncid,varid=varid,ndims=ndims,dimids=dimids))

  write(*,'(A)') TRIM(varname)
  write(*,'(A,I2)') '   ndims:',ndims
  write(*,*) '   dimids:',dimids

  CALL check(nf90_get_var(ncid,varid,idata))
  !Close netCDF file
  !:-------:-------:-------:-------:-------:-------:-------:-------:
  CALL check(nf90_close(ncid))

END SUBROUTINE read_nc_1d

SUBROUTINE read_nc_2d(infile,varname,idata,nx,ny)
  use iso_c_binding, only : c_float
  use netcdf

  REAL(c_float), DIMENSION(nx,ny), INTENT(OUT) :: idata
  INTEGER(KIND=4), INTENT(IN) :: nx,ny
  INTEGER(KIND=4), DIMENSION(2) :: dimids
  INTEGER(KIND=4) :: ncid, ndims, varid
  CHARACTER(LEN=1024), INTENT(IN) :: infile,varname
  !Open netCDF file
  !:-------:-------:-------:-------:-------:-------:-------:-------:
  CALL check(nf90_open(TRIM(infile), nf90_nowrite, ncid))
  CALL check(nf90_inq_varid(ncid,TRIM(varname),varid))
  CALL check(nf90_inquire_variable(ncid=ncid,varid=varid,ndims=ndims,dimids=dimids))

  write(*,'(A)') TRIM(varname)
  write(*,'(A,I2)') '   ndims:',ndims
  write(*,*) '   dimids:',dimids

  CALL check(nf90_get_var(ncid,varid,idata))
  !Close netCDF file
  !:-------:-------:-------:-------:-------:-------:-------:-------:
  CALL check(nf90_close(ncid))

END SUBROUTINE read_nc_2d

SUBROUTINE mean_2d(input, mean,nx,ny)
  use iso_c_binding, only : c_float
  INTEGER, INTENT(IN) :: nx,ny
  REAL(c_float), INTENT(IN) :: input(nx,ny)
  REAL(c_float), INTENT(OUT):: mean

  INTEGER :: i,j

  mean = 0.0
  DO i=1,nx
    DO j=1,ny
      mean = mean + input(i,j)
    ENDDO
  ENDDO

  mean = mean / REAL(nx*ny)

END SUBROUTINE mean_2d


SUBROUTINE check(istatus)
  use netcdf

  INTEGER, INTENT (IN) :: istatus
  IF (istatus /= nf90_noerr) THEN
  write(*,'(A)') TRIM((nf90_strerror(istatus)))
  END IF
END SUBROUTINE check

