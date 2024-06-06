program use_raphrrr_sfc
        !
  use kinds, only: i_kind,r_kind,r_single,i_byte
  use module_ncio, only : ncio
  use module_map_utils, only : map_util
  use module_surface, only : use_surface 
  use mpi
!
  implicit none
!
  type(ncio)     :: raphrrr,mpas
  type(map_util) :: map
  type(use_surface) :: sfc
!
!  namelist files
!
  character*80 :: rapfile
  character*80 :: hrrrfile
  character*80 :: hrrr_akfile
  character*80 :: mpasfile
  character*80 :: mpasfile_read
  logical :: do_lake_surgery
  namelist/setup/ rapfile,hrrrfile,hrrr_akfile,mpasfile,do_lake_surgery
!
! MPI variables
  integer :: npe, mype, mypeLocal,ierror
!
! define grid
  integer :: ncells
  integer :: nx_rap,ny_rap
  integer :: nz_mpas,nz_raphrrr
  integer :: nz_mpas_lake,nz_raphrrr_lake
  integer :: nz_mpas_snow,nz_raphrrr_snow
!
! define map
  real(r_single),allocatable,target :: rlon2d_mpas(:),rlat2d_mpas(:)
  real(r_single),allocatable,target :: rlon2d_raphrrr(:,:),rlat2d_raphrrr(:,:)
  integer(i_byte),allocatable :: landmask_raphrrr(:,:)
  integer(i_byte),allocatable :: landmask_mpas(:)
  integer(i_byte),allocatable :: lakemask_raphrrr(:,:)
  integer(i_byte),allocatable :: lakemask_mpas(:)
  real(r_single),allocatable,target :: tmp1d4b(:)
  integer(i_kind),allocatable,target :: tmp1d4i(:)
  real(r_single),allocatable,target :: tmp2d4b(:,:)

  real(r_single),parameter:: rad2deg=180.0/3.1415926
  integer :: i,j,k,n
  character*80 :: raphrrrfile
!
!**********************************************************************
!
!            END OF DECLARATIONS....start of program
! MPI setup
  call MPI_INIT(ierror)
  call MPI_COMM_SIZE(mpi_comm_world,npe,ierror)
  call MPI_COMM_RANK(mpi_comm_world,mype,ierror)

!
  if(mype==0) then
!
     rapfile='missing'
     hrrrfile='missing'
     hrrr_akfile='missing'
     mpasfile='missing'
     do_lake_surgery=.false.
     open(15, file='use_raphrrr_sfc.namelist')
        read(15,setup)
     close(15)
     write(*,setup)

!     mpasfile_read=trim(mpasfile)//"_read"
     mpasfile_read=trim(mpasfile)
! read in mpas latlon
     call mpas%open(trim(mpasfile_read),"r",200)
     call mpas%get_dim("nCells",ncells)
     write(*,*) 'ncells=',ncells

     allocate(rlon2d_mpas(ncells))
     allocate(rlat2d_mpas(ncells))
     call mpas%get_var("lonCell",ncells,rlon2d_mpas)
     call mpas%get_var("latCell",ncells,rlat2d_mpas)
     rlon2d_mpas=rlon2d_mpas*rad2deg
     rlat2d_mpas=rlat2d_mpas*rad2deg
!
! read in mpas land mask
     allocate(landmask_mpas(ncells))
     allocate(tmp1d4i(ncells))
     call mpas%get_var("landmask",ncells,tmp1d4i)
     do i=1,ncells
        landmask_mpas(i)=tmp1d4i(i)
        if(landmask_mpas(i) >=2 ) landmask_mpas(i)=1
     enddo

     call mpas%get_dim("nSoilLevels",nz_mpas)

! read in mpas lake mask
     if(do_lake_surgery) then
       allocate(lakemask_mpas(ncells))
       call mpas%get_var("clm_lake_initialized",ncells,tmp1d4i)
         do i=1,ncells
            lakemask_mpas(i)=int(tmp1d4i(i))
         enddo
        call mpas%get_dim("levlake_clm_lake",nz_mpas_lake)
        call mpas%get_dim("levsnowsoil1_clm_lake",nz_mpas_snow)
     else
        nz_mpas_lake=-99  
        nz_mpas_snow=-99
     endif

     deallocate(tmp1d4i)
     call mpas%close()
     write(*,*) "landmask=",maxval(landmask_mpas), minval(landmask_mpas)
     write(*,*) "nz_mpas=",nz_mpas
     if(do_lake_surgery) then
        write(*,*) "lakemask=",maxval(lakemask_mpas), minval(lakemask_mpas)
        write(*,*) "nz_mpas_lake=",nz_mpas_lake
        write(*,*) "nz_mpas_snow=",nz_mpas_snow
     endif
!
! use RAP
! read in rap latlon
     do n=1,3

        raphrrrfile='missing'
        if(n==1 .and. trim(rapfile) /= 'missing' ) then
           raphrrrfile=rapfile
        elseif(n==2 .and. trim(hrrrfile) /= 'missing' ) then
           raphrrrfile=hrrrfile
        elseif(n==3 .and. trim(hrrr_akfile) /= 'missing' ) then
           raphrrrfile=hrrr_akfile
        else
           cycle
        endif
        write(*,*) "===>"
        write(*,*) "===> tansfer surface fields from ",trim(raphrrrfile)
        write(*,*) "===>"
        call raphrrr%open(trim(raphrrrfile),"r",200)
        call raphrrr%get_dim("west_east",nx_rap)
        call raphrrr%get_dim("south_north",ny_rap)
        call raphrrr%get_dim("soil_layers_stag",nz_raphrrr)
        write(*,*) 'nx_rap,ny_rap=',nx_rap,ny_rap,nz_raphrrr
        if(nz_raphrrr /= nz_mpas ) then
           write(*,*) "Error in vertical level =", nz_raphrrr,nz_mpas
           stop 123
        endif

        if(do_lake_surgery) then
           call raphrrr%get_dim("soil_levels_or_lake_levels_stag",nz_raphrrr_lake)
           call raphrrr%get_dim("snow_and_soil_levels_stag",nz_raphrrr_snow)
           write(*,*) 'nz_raphrrr_lake,nz_raphrrr_snow=',nz_raphrrr_lake,nz_raphrrr_snow
           if(nz_raphrrr_lake /= nz_mpas_lake .or. &
              nz_raphrrr_snow /= nz_mpas_snow) then
              write(*,*) "Error in vertical level lake=", nz_raphrrr_lake,nz_mpas_lake
              write(*,*) "Error in vertical level snow=", nz_raphrrr_snow,nz_mpas_snow
              stop 123
           endif
        else
           nz_raphrrr_lake=-99
           nz_raphrrr_snow=-99
        endif
!
        allocate(rlon2d_raphrrr(nx_rap,ny_rap))
        allocate(rlat2d_raphrrr(nx_rap,ny_rap))
        call raphrrr%get_var("XLONG",nx_rap,ny_rap,rlon2d_raphrrr)
        call raphrrr%get_var("XLAT",nx_rap,ny_rap,rlat2d_raphrrr)
        call map%init_general_transform(nx_rap,ny_rap,rlat2d_raphrrr,rlon2d_raphrrr)
        deallocate(rlon2d_raphrrr)
        deallocate(rlat2d_raphrrr)

        allocate(landmask_raphrrr(nx_rap,ny_rap))
        allocate(tmp2d4b(nx_rap,ny_rap))
        call raphrrr%get_var("LANDMASK",nx_rap,ny_rap,tmp2d4b)
        do j=1,ny_rap
          do i=1,nx_rap
             landmask_raphrrr(i,j)=int(tmp2d4b(i,j))
          enddo
        enddo

        if(do_lake_surgery) then
          allocate(lakemask_raphrrr(nx_rap,ny_rap))
          call raphrrr%get_var("LAKEMASK",nx_rap,ny_rap,tmp2d4b)
          do j=1,ny_rap
            do i=1,nx_rap
               lakemask_raphrrr(i,j)=int(tmp2d4b(i,j))
            enddo
          enddo
        endif

        call raphrrr%get_var("SNOW",nx_rap,ny_rap,tmp2d4b)
        do j=1,ny_rap
          do i=1,nx_rap
             if(tmp2d4b(i,j) > 0.01 ) landmask_raphrrr(i,j)=2  ! snow coverage
          enddo
        enddo
        deallocate(tmp2d4b)
        call raphrrr%close()
!
! initial sfc and map index
        call sfc%init(ncells,nz_mpas,nz_mpas_lake,nz_mpas_snow,nx_rap,ny_rap,4)
        call sfc%build_mapindex(map,rlon2d_mpas,rlat2d_mpas,landmask_mpas,landmask_raphrrr)
        if(do_lake_surgery) &
           call sfc%build_lakeindex(map,rlon2d_mpas,rlat2d_mpas,lakemask_mpas,lakemask_raphrrr)
        call sfc%set_varname()
        call sfc%use_sfc(raphrrrfile,mpasfile,mpasfile_read)
        if(do_lake_surgery) &
           call sfc%use_lake(raphrrrfile,mpasfile,mpasfile_read)

        if(n==2) call sfc%remove_snow(raphrrrfile,mpasfile,rlat2d_mpas)
        call sfc%close()
! release memory
        deallocate(landmask_raphrrr)
        if(do_lake_surgery) deallocate(lakemask_raphrrr)
        call map%destory_general_transform()

     enddo ! n
!
     deallocate(landmask_mpas)
     if(do_lake_surgery)  deallocate(lakemask_mpas)
!
     write(6,*) "=== USE_RAPHRRR_SFC REPROCCESS SUCCESS ==="

  endif ! mype==0

  call MPI_FINALIZE(ierror)
!
end program use_raphrrr_sfc
