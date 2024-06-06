!-------------------------------------------------------------------------
!    NOAA/NCEP, National Centers for Environmental Prediction GSI        !
!-------------------------------------------------------------------------
!BOP
!
module module_surface

  use kinds, only: r_kind,r_single,i_kind,i_byte,i_short
  implicit none

  public :: use_surface
!
!-------------------------------------------------------------------------

! set default to private
  private
  type :: use_surface
      integer :: nvar,nvarlake
      integer :: nvar3d,nvar3dlake,nvar3dsnow
      character(len=20),allocatable :: var_mpas(:),var_rap(:)
      character(len=20),allocatable :: var_mpas_lake(:),var_rap_lake(:)
      integer :: halo
      integer :: ncells
      integer :: nlev,nlev_snow,nlev_lake
      integer :: nlat_target,nlon_target
      integer(i_short), allocatable :: index_x(:)
      integer(i_short), allocatable :: index_y(:)
      integer(i_short), allocatable :: index_x_nomatch(:)
      integer(i_short), allocatable :: index_y_nomatch(:)
      integer(i_short), allocatable :: lake_index_x(:)
      integer(i_short), allocatable :: lake_index_y(:)
    contains
      procedure :: init
      procedure :: build_mapindex
      procedure :: build_lakeindex
      procedure :: set_varname
      procedure :: use_sfc
      procedure :: use_lake
      procedure :: remove_snow
      procedure :: close
  end type use_surface
!
! constants
!
contains
   
  subroutine set_varname(this)
!                .      .    .                                       .
! subprogram:   build_mapindex
!   prgmmr:
!
! abstract:
!
! program history log:
!
!   input argument list:
!   output argument list:
!
!
! 3D fields
! 1. soil moisture SMOIS (9 levels) to smois (9 levels)
! 2. Liquid soil water SH2O (9 levels) to sh2o (9 levels). 
! 3. Soil temperature TSLB (9 levels) to tslb(top 8 levels) and tiice (top 8 levels)
!
! 2D fields
! 4. TSK to tsfc, and tsfcl,tsea,tref,tisfc
! 5, VEGFRA to vegfra
! 8. SNOW to snow
! 9. SNOWH to snowh
! 10. SNOWC to snowc
!
    implicit none

    class(use_surface) :: this
!
    this%nvar=10
    this%nvar3d=3
!
    if(allocated(this%var_mpas)) deallocate(this%var_mpas)
    allocate(this%var_mpas(this%nvar))
    if(allocated(this%var_rap)) deallocate(this%var_rap)
    allocate(this%var_rap(this%nvar))
!
! Soil moisture: 3D float
    this%var_rap(1)="SMOIS"
    this%var_mpas(1)="smois"

! Liquid soil water: 3D float
    this%var_rap(2)="SH2O"
    this%var_mpas(2)="sh2o"

! Soil temperature: 3D float
    this%var_rap(3)="TSLB"
    this%var_mpas(3)="tslb"

! SURFACE SKIN TEMPERATURE: 2D float
    this%var_rap(4)="TSK"
    this%var_mpas(4)="skintemp" 

! VEGETATION FRACTION: RAP/HRRR has real-time VIIRS GVF, MPAS - climatology
    this%var_rap(5)="VEGFRA"
    this%var_mpas(5)="vegfra"

    this%var_rap(6)="SHDMAX"
    this%var_mpas(6)="shdmax"

    this%var_rap(7)="SHDMIN"
    this%var_mpas(7)="shdmin"

! SNOW WATER EQUIVALENT [mm]
    this%var_rap(8)="SNOW"
    this%var_mpas(8)="snow"

! PHYSICAL SNOW DEPTH [m]:
    this%var_rap(9)="SNOWH"
    this%var_mpas(9)="snowh"

    this%var_rap(10)="SNOWC"
    this%var_mpas(10)="snowc"
!
!    Add lake model fields
! We have to read the variable clm_lake_initialized. If it is =1, then it is a lake in the RRFS. 
! If this point is not a lake in RAP/HRRR we keep the values from cold-start run.
! The RAP/HRRR variables go to:
! 3dsoil   1  float T_LAKE3D --> lake_t_lake3d
! 3dsoil   2  float LAKE_ICEFRAC3D --> lake_icefrac3d, fice (2d) 
! 3dsnow   3  float T_SOISNO3D   --> lake_t_soisno3d
! 3dsnow   4  float H2OSOI_ICE3D --> lake_h2osoi_ice3d
! 3dsnow   5  float H2OSOI_LIQ3D --> lake_h2osoi_liq3d
! 3dsnow   6  float H2OSOI_VOL3D --> lake_h2osoi_vol3d
! 3dsnow   7  float DZ3D --> lake_snow_dz3d
! 3dsnow   8  float Z3D --> lake_snow_z3d
! 3dsnow+1 9  float ZI3D --> lake_snow_zi3d
!         10  float SAVEDTKE12D  --> lake_savedtke12d
!         11  float SNOWDP2D --> lake_sndpth2d
!         12  float SNOWDP2D --> snodi
!         13  float H2OSNO2D --> lake_h2osno2d
!         14  float H2OSNO2D --> weasdi
!         15  float T_GRND2D --> lake_tsfc
!         16  float T_GRND2D --> tisfc
!         17  float T_GRND2D --> tsea
!         18  float T_GRND2D --> T_snow
!         19  float T_GRND2D --> tsfc
!         20  float T_GRND2D --> tsnow_ice
!         21  float SNL2D    --> lake_snl2d
!
    this%nvarlake=21
    this%nvar3dlake=2
    this%nvar3dsnow=7
    if(allocated(this%var_mpas_lake)) deallocate(this%var_mpas_lake)
    allocate(this%var_mpas_lake(this%nvarlake))
    if(allocated(this%var_rap_lake)) deallocate(this%var_rap_lake)
    allocate(this%var_rap_lake(this%nvarlake))
!
    this%var_rap_lake(1)="T_LAKE3D"
    this%var_mpas_lake(1)="lake_t_lake3d"

    this%var_rap_lake(2)="LAKE_ICEFRAC3D"
    this%var_mpas_lake(2)="lake_icefrac3d"

    this%var_rap_lake(3)="T_SOISNO3D"
    this%var_mpas_lake(3)="lake_t_soisno3d"

    this%var_rap_lake(4)="H2OSOI_ICE3D"
    this%var_mpas_lake(4)="lake_h2osoi_ice3d"

    this%var_rap_lake(5)="H2OSOI_LIQ3D"
    this%var_mpas_lake(5)="lake_h2osoi_liq3d"

    this%var_rap_lake(6)="H2OSOI_VOL3D"
    this%var_mpas_lake(6)="lake_h2osoi_vol3d"

    this%var_rap_lake(7)="DZ3D"
    this%var_mpas_lake(7)="lake_snow_dz3d"

    this%var_rap_lake(8)="Z3D"
    this%var_mpas_lake(8)="lake_snow_z3d"

    this%var_rap_lake(9)="ZI3D"
    this%var_mpas_lake(9)="lake_snow_zi3d"

    this%var_rap_lake(10)="SAVEDTKE12D"
    this%var_mpas_lake(10)="lake_savedtke12d"

    this%var_rap_lake(11)="SNOWDP2D"
    this%var_mpas_lake(11)="lake_sndpth2d"
    this%var_rap_lake(12)="SNOWDP2D"
    this%var_mpas_lake(12)="snodi"

    this%var_rap_lake(13)="H2OSNO2D"
    this%var_mpas_lake(13)="lake_h2osno2d"
    this%var_rap_lake(14)="H2OSNO2D"
    this%var_mpas_lake(14)="weasdi"

    this%var_rap_lake(15)="T_GRND2D"
    this%var_mpas_lake(15)="lake_tsfc"
    this%var_rap_lake(16)="T_GRND2D"
    this%var_mpas_lake(16)="tisfc"
    this%var_rap_lake(17)="T_GRND2D"
    this%var_mpas_lake(17)="tsea"
    this%var_rap_lake(18)="T_GRND2D"
    this%var_mpas_lake(18)="T_snow"
    this%var_rap_lake(19)="T_GRND2D"
    this%var_mpas_lake(19)="tsfc"
    this%var_rap_lake(20)="T_GRND2D"
    this%var_mpas_lake(20)="tsnow_ice"

    this%var_rap_lake(21)="SNL2D"
    this%var_mpas_lake(21)="lake_snl2d"
!
  end subroutine set_varname

  subroutine use_lake(this,rapfile,mpasfile,mpasfile_read)
!                .      .    .                                       .
! subprogram:   build_mapindex
!   prgmmr:
!
! abstract:
!
! program history log:
!
!   input argument list:
!
!   output argument list:
!
    use module_ncio, only : ncio
  
    implicit none

    character*80,intent(in) :: rapfile
    character*80,intent(in) :: mpasfile
    character*80,intent(in) :: mpasfile_read

    class(use_surface) :: this
    type(ncio)     :: raphrrr,mpas,mpasr
!
    real(r_single),allocatable,target :: fld3d4b(:,:,:)
!
    real(r_single),allocatable,target :: tmp2d4b(:,:)
    real(r_single),allocatable,target :: tmp3d4b(:,:,:)
    real(r_single),  allocatable,target :: tmp2d4br(:)
    real(r_single),  allocatable,target :: tmp3d4br(:,:)
!
    integer  :: nx_rap,ny_rap
    integer  :: ncells_mpas,nz_mpas
    integer  :: i,j,ix,jx
    character(len=20) :: thisvar_mpas,thisvar_rap
    integer  :: k,kx 
!

    ncells_mpas=this%ncells
    nx_rap=this%nlon_target
    ny_rap=this%nlat_target
!
    call raphrrr%open(trim(rapfile),"r",200)
    call mpasr%open(trim(mpasfile_read),"r",200)
    call mpas%open(trim(mpasfile),"w",200)
!
    do k=1,this%nvarlake
       thisvar_mpas=this%var_mpas_lake(k)
       thisvar_rap=this%var_rap_lake(k)
       write(*,*) "=============================================="
       write(*,'(I4,4a)') k," Working to replace mpas lake variable ",trim(thisvar_mpas), &
                  " from RAP/HRRR ",trim(thisvar_rap)
       if(k <= (this%nvar3dlake+this%nvar3dsnow)) then
          if(k<=this%nvar3dlake) then
             nz_mpas=this%nlev_lake
          else
             nz_mpas=this%nlev_snow
             if(k==9) nz_mpas=this%nlev_snow+1
          endif
          allocate(tmp3d4b(nx_rap,ny_rap,nz_mpas))
          call raphrrr%get_var(trim(thisvar_rap),nx_rap,ny_rap,nz_mpas,tmp3d4b)

          allocate(tmp3d4br(nz_mpas,ncells_mpas))
          call mpasr%get_var(trim(thisvar_mpas),nz_mpas,ncells_mpas,tmp3d4br)

            do i=1,ncells_mpas
               ix=this%lake_index_x(i)
               jx=this%lake_index_y(i)
               if( (ix >= 1 .and. ix <= nx_rap) .and. &
                   (jx >= 1 .and. jx <= ny_rap) ) then
                  do kx=1,nz_mpas
                     tmp3d4br(kx,i)=tmp3d4b(ix,jx,kx)
                  enddo
               endif
            enddo
          if(trim(thisvar_mpas) == "lake_icefrac3d") then
             allocate(tmp2d4br(ncells_mpas))
             do i=1,ncells_mpas
                ix=this%lake_index_x(i)
                jx=this%lake_index_y(i)
                if( (ix >= 1 .and. ix <= nx_rap) .and. &
                    (jx >= 1 .and. jx <= ny_rap) ) then
                   tmp2d4br(i)=tmp3d4b(ix,jx,1)
                endif
             enddo
          endif
          deallocate(tmp3d4b)

          do kx = 1,nz_mpas
            write(*,*) k,maxval(tmp3d4br(:,kx)),minval(tmp3d4br(:,kx))
          enddo
          call mpas%replace_var(trim(thisvar_mpas),nz_mpas,ncells_mpas,tmp3d4br)
          if(trim(thisvar_mpas) == "lake_icefrac3d") then
             call mpas%replace_var("fice",ncells_mpas,tmp2d4br(:))
             deallocate(tmp2d4br)
          endif
          deallocate(tmp3d4br)
       else
          allocate(tmp2d4b(nx_rap,ny_rap))
          call raphrrr%get_var(trim(thisvar_rap),nx_rap,ny_rap,tmp2d4b)

          allocate(tmp2d4br(ncells_mpas))
          call mpasr%get_var(trim(thisvar_mpas),ncells_mpas,tmp2d4br)
        
             do i=1,ncells_mpas
                ix=this%lake_index_x(i)
                jx=this%lake_index_y(i)
                if( (ix >= 1 .and. ix <= nx_rap) .and. &
                    (jx >= 1 .and. jx <= ny_rap) ) then
                  tmp2d4br(i)=tmp2d4b(ix,jx)
                endif
             enddo
          deallocate(tmp2d4b)

          write(*,*) k,maxval(tmp2d4br),minval(tmp2d4br)
          call mpas%replace_var(trim(thisvar_mpas),ncells_mpas,tmp2d4br)
          deallocate(tmp2d4br)
!
       endif
    enddo
!
    call mpas%close()
    call mpasr%close()
    call raphrrr%close()

  end subroutine use_lake

  subroutine use_sfc(this,rapfile,mpasfile,mpasfile_read)
!                .      .    .                                       .
! subprogram:   build_mapindex
!   prgmmr:
!
! abstract:
!
! program history log:
!
!   input argument list:
!
!   output argument list:
!
    use module_ncio, only : ncio
  
    implicit none

    character*80,intent(in) :: rapfile
    character*80,intent(in) :: mpasfile
    character*80,intent(in) :: mpasfile_read

    class(use_surface) :: this
    type(ncio)     :: raphrrr,mpas,mpasr
!
!    real(r_single),allocatable,target :: fld2d4b(:,:)
    real(r_single),allocatable,target :: fld3d4b(:,:,:)
!
    real(r_single),allocatable,target :: tmp2d4b(:,:)
    real(r_single),allocatable,target :: tmp3d4b(:,:,:)
    real(r_single),  allocatable,target :: tmp2d4br(:)
    real(r_single),  allocatable,target :: tmp3d4br(:,:)
!    integer(i_byte),allocatable :: landmask_mpas(:,:)
!
    integer  :: nx_rap,ny_rap
    integer  :: ncells_mpas,nz_mpas
    integer  :: i,j,ix,jx
    character(len=20) :: thisvar_mpas,thisvar_rap
    integer  :: k,kx 
!

    ncells_mpas=this%ncells
    nz_mpas=this%nlev
    nx_rap=this%nlon_target
    ny_rap=this%nlat_target
!
    call raphrrr%open(trim(rapfile),"r",200)
    call mpasr%open(trim(mpasfile_read),"r",200)
    call mpas%open(trim(mpasfile),"w",200)
!
    do k=1,this%nvar
       thisvar_mpas=this%var_mpas(k)
       thisvar_rap=this%var_rap(k)
       write(*,*) "=============================================="
       write(*,'(I4,4a)') k," Working to replace mpas variable ",trim(thisvar_mpas), &
                  " from RAP/HRRR ",trim(thisvar_rap)
       if(k <= this%nvar3d) then
          allocate(tmp3d4b(nx_rap,ny_rap,nz_mpas))
          call raphrrr%get_var(trim(thisvar_rap),nx_rap,ny_rap,nz_mpas,tmp3d4b)

          allocate(tmp3d4br(nz_mpas,ncells_mpas))
          call mpasr%get_var(trim(thisvar_mpas),nz_mpas,ncells_mpas,tmp3d4br)

          do i=1,ncells_mpas
               ix=this%index_x(i)
               jx=this%index_y(i)
               if( (ix >= 1 .and. ix <= nx_rap) .and. &
                   (jx >= 1 .and. jx <= ny_rap) ) then
!                  if(trim(thisvar_mpas) == "tslb") then
!                     do kx=1,nz_mpas-1
!                        tmp3d4br(ks,i)=tmp3d4b(ix,jx,kx)
!                     enddo
!                  else
                     do kx=1,nz_mpas
                        tmp3d4br(kx,i)=tmp3d4b(ix,jx,kx)
                     enddo
!                  endif
               endif
          enddo
          deallocate(tmp3d4b)

          call mpas%replace_var(trim(thisvar_mpas),nz_mpas,ncells_mpas,tmp3d4br)
          deallocate(tmp3d4br)
       else
          allocate(tmp2d4b(nx_rap,ny_rap))
          call raphrrr%get_var(trim(thisvar_rap),nx_rap,ny_rap,tmp2d4b)

          allocate(tmp2d4br(ncells_mpas))
          call mpasr%get_var(trim(thisvar_mpas),ncells_mpas,tmp2d4br)
        
          if((maxval(tmp2d4br) - minval(tmp2d4br) )< 1.0e-10) then
              do i=1,ncells_mpas
                ix=this%index_x_nomatch(i)
                jx=this%index_y_nomatch(i)
                if( (ix >= 1 .and. ix <= nx_rap) .and. &
                    (jx >= 1 .and. jx <= ny_rap) ) then
                  tmp2d4br(i)=tmp2d4b(ix,jx)
                endif
              enddo
          else
              do i=1,ncells_mpas
                ix=this%index_x(i)
                jx=this%index_y(i)
                if( (ix >= 1 .and. ix <= nx_rap) .and. &
                    (jx >= 1 .and. jx <= ny_rap) ) then
                  tmp2d4br(i)=tmp2d4b(ix,jx)
                endif
              enddo
          endif
          deallocate(tmp2d4b)

          if(trim(thisvar_mpas) == "snowc") then
             do i=1,ncells_mpas
                tmp2d4br(i)=min(max(tmp2d4br(i),0.0),1.0)
             enddo
          elseif(trim(thisvar_mpas) == "snow") then
             do i=1,ncells_mpas
                tmp2d4br(i)=min(tmp2d4br(i),3000.0)
             enddo
          elseif(trim(thisvar_mpas) == "showh") then
             do i=1,ncells_mpas
                tmp2d4br(i)=min(tmp2d4br(i),7.5)
             enddo
          endif

          call mpas%replace_var(trim(thisvar_mpas),ncells_mpas,tmp2d4br)
          deallocate(tmp2d4br)
!
       endif
    enddo
!
    call mpas%close()
    call mpasr%close()
    call raphrrr%close()

  end subroutine use_sfc

  subroutine remove_snow(this,rapfile,mpasfile,rlat)
!                .      .    .                                       .
! subprogram:   remove snow
!
! To remove RRFS snow in the southern regions:
! if HRRR snow = 0 south of 37N, then set in RRFS sheleg=weasdl=snwdph=snodl = 0.

!   prgmmr:
!
! abstract:
!
! program history log:
!
!   input argument list:
!
!   output argument list:
!
    use module_ncio, only : ncio
  
    implicit none

    class(use_surface) :: this

    character*80,intent(in) :: rapfile
    character*80,intent(in) :: mpasfile
    real(r_single),  intent(in) :: rlat(this%ncells)
!
    type(ncio)     :: raphrrr,mpas
!
    real(r_single),allocatable,target :: tmp2d4b(:,:)
    real(r_single),allocatable,target :: tmp2d4br(:)
!
    integer  :: nx_rap,ny_rap
    integer  :: ncells_mpas,nz_mpas
    integer  :: i,j,ix,jx,nn
    character(len=20) :: thisvar_mpas
!

    ncells_mpas=this%ncells
    nz_mpas=this%nlev
    nx_rap=this%nlon_target
    ny_rap=this%nlat_target
!
    allocate(tmp2d4b(nx_rap,ny_rap))
    allocate(tmp2d4br(ncells_mpas))

    call raphrrr%open(trim(rapfile),"r",200)
    call raphrrr%get_var("SNOW",nx_rap,ny_rap,tmp2d4b)
    call raphrrr%close()
!
    do nn=1,3
       if(nn==1) thisvar_mpas="snow"
       if(nn==2) thisvar_mpas="snowh"
       if(nn==3) thisvar_mpas="snowc"

       call mpas%open(trim(mpasfile),"r",200)
       call mpas%get_var(trim(thisvar_mpas),ncells_mpas,tmp2d4br)
       call mpas%close()

       do i=1,ncells_mpas
          ix=this%index_x(i)
          jx=this%index_y(i)
          if( (ix >= 1 .and. ix <= nx_rap) .and. &
              (jx >= 1 .and. jx <= ny_rap) ) then
             if(tmp2d4b(ix,jx) < 1.0e-10_r_single .and. rlat(i) < 37.0 ) then
                tmp2d4br(i)=0.0
             endif
          endif
       enddo

       call mpas%open(trim(mpasfile),"w",200)
       call mpas%replace_var(trim(thisvar_mpas),ncells_mpas,tmp2d4br)
       call mpas%close()

    end do ! nn
    deallocate(tmp2d4br)
    deallocate(tmp2d4b)
  
  end subroutine remove_snow

  subroutine build_lakeindex(this,map,rlon,rlat,lakemask_this,lakemask_target)
!                .      .    .                                       .
! subprogram:   build_mapindex
!   prgmmr:
!
! abstract:
!
! program history log:
!
!   input argument list:
!    map: grid convert class
!    lakemask_this: lake mask for this grid
!    lakemask_target: lake mask for target grid
!
!   output argument list:
!
    use module_map_utils, only : map_util
    implicit none

    class(use_surface) :: this
    type(map_util),  intent(in) :: map
    integer(i_byte), intent(in) :: lakemask_this(this%ncells)
    integer(i_byte), intent(in) :: lakemask_target(this%nlon_target,this%nlat_target)
    real(r_single),  intent(in) :: rlon(this%ncells)
    real(r_single),  intent(in) :: rlat(this%ncells)
!
    real(r_single) :: xc,yc
    integer  :: ixc,jyc
    integer  :: i,j,ix,jx
    integer  :: ii,jj,ixx,jxx
    real(r_single) :: dist,thisdist
!
    do i=1,this%ncells
       call map%tll2xy(rlon(i),rlat(i),xc,yc)
       ixc=int(xc+0.5)
       jyc=int(yc+0.5)
       if( (ixc > 0 .and. ixc < this%nlon_target + 1) .and. &
           (jyc > 0 .and. jyc < this%nlat_target + 1) ) then
          ix=min(max(ixc,1),this%nlon_target)
          jx=min(max(jyc,1),this%nlat_target)
! lake in RRFS
          if( lakemask_this (i) == 1 ) then
             if( lakemask_target(ix,jx) == 1 ) then
! lake in RAP HRRR
                this%lake_index_x(i)=ix
                this%lake_index_y(i)=jx 
             else
! deal with not matched lake: if there is a lake within a grid point, count as
!                             match
                dist=99999.0
                do jj=max(1,jx-1),min(jx+1,this%nlat_target)
                do ii=max(1,ix-1),min(ix+1,this%nlon_target)
                   if(lakemask_this(i)==lakemask_target(ii,jj)) then
                      thisdist=(jj-jx)*(jj-jx)+(ii-ix)*(ii-ix)
                      thisdist=sqrt(thisdist)
                      if(thisdist <= dist) then
                         dist=thisdist
                         ixx=ii
                         jxx=jj
                      endif
                   endif
                enddo
                enddo
                if(dist < 88888.0) then
                   this%lake_index_x(i)=ixx
                   this%lake_index_y(i)=jxx
                  !write(*,*) i,this%lake_index_x(i),this%lake_index_y(i)
                endif
             endif
          endif  ! if RRFS lake
       endif  ! if inside the traget domain
    enddo

  end subroutine build_lakeindex

  subroutine build_mapindex(this,map,rlon,rlat,landmask_this,landmask_target)
!                .      .    .                                       .
! subprogram:   build_mapindex
!   prgmmr:
!
! abstract:
!
! program history log:
!
!   input argument list:
!    map: grid convert class
!    landmask_this: land mask for this grid
!    landmask_target: land mask for target grid
!
!   output argument list:
!
    use module_map_utils, only : map_util
    implicit none

    class(use_surface) :: this
    type(map_util),  intent(in) :: map
    integer(i_byte), intent(in) :: landmask_this(this%ncells)
    integer(i_byte), intent(in) :: landmask_target(this%nlon_target,this%nlat_target)
    real(r_single),  intent(in) :: rlon(this%ncells)
    real(r_single),  intent(in) :: rlat(this%ncells)
!
    real(r_single) :: xc,yc
    integer  :: ixc,jyc
    integer  :: i,j,ix,jx
    integer  :: ii,jj,ixx,jxx
    real(r_single) :: dist,thisdist
!
    do i=1,this%ncells
       call map%tll2xy(rlon(i),rlat(i),xc,yc)
       ixc=int(xc+0.5)
       jyc=int(yc+0.5)
       if( (ixc > 0 .and. ixc < this%nlon_target + 1) .and. &
           (jyc > 0 .and. jyc < this%nlat_target + 1) ) then
          ix=min(max(ixc,1),this%nlon_target)
          jx=min(max(jyc,1),this%nlat_target)
          this%index_x_nomatch(i)=ix
          this%index_y_nomatch(i)=jx 
          if( landmask_target(ix,jx) == 2 ) then
! snow/ice in RAP HRRR
             this%index_x(i)=ix
             this%index_y(i)=jx 
          elseif( landmask_this(i)==landmask_target(ix,jx) ) then
! match water and land
             this%index_x(i)=ix
             this%index_y(i)=jx 
          else
! deal with not matched water and land
             dist=99999.0
             do jj=max(1,jx-this%halo),min(jx+this%halo,this%nlat_target)
               do ii=max(1,ix-this%halo),min(ix+this%halo,this%nlon_target)
                 if(landmask_this(i)==landmask_target(ii,jj)) then
                    thisdist=(jj-jx)*(jj-jx)+(ii-ix)*(ii-ix)
                    thisdist=sqrt(thisdist)
                    if(thisdist <= dist) then
                       dist=thisdist
                       ixx=ii
                       jxx=jj
                    endif
                 endif
               enddo
             enddo
             if(dist < 88888.0) then
               this%index_x(i)=ixx
               this%index_y(i)=jxx
               !write(*,*) i,this%index_x(i),this%index_y(i)
             endif
          endif
       endif
    enddo
!     
!
  end subroutine build_mapindex

  subroutine init(this,ncells,nlev,nlev_lake,nlev_snow,nlon_target,nlat_target,halo)
!                .      .    .                                       .
! subprogram:    init
!   prgmmr:
!
! abstract:
!
! program history log:
!
!   input argument list:
!    nlon,nlat - grid dimension
!
!   output argument list:
!
    implicit none

    class(use_surface) :: this
    integer,intent(in) :: halo
    integer,intent(in) :: ncells,nlev,nlev_lake,nlev_snow
    integer,intent(in) :: nlon_target,nlat_target

    this%halo=halo
    this%ncells=ncells
    this%nlev=nlev
    this%nlev_lake=nlev_lake
    this%nlev_snow=nlev_snow
    this%nlon_target=nlon_target
    this%nlat_target=nlat_target
    allocate(this%index_x(this%ncells))
    allocate(this%index_y(this%ncells))
    this%index_x=-999
    this%index_y=-999
    allocate(this%index_x_nomatch(this%ncells))
    allocate(this%index_y_nomatch(this%ncells))
    this%index_x_nomatch=-999
    this%index_y_nomatch=-999
    allocate(this%lake_index_x(this%ncells))
    allocate(this%lake_index_y(this%ncells))
    this%lake_index_x=-999
    this%lake_index_y=-999

  end subroutine init

  subroutine close(this)
!                .      .    .                                       .
! subprogram:    close
!   prgmmr:
!
! abstract:
!
! program history log:
!
!   input argument list:
!
!   output argument list:
!
    implicit none

    class(use_surface) :: this

    this%nvar=0
    this%nvar3d=0
    this%nvarlake=0
    this%nvar3dlake=0
    this%nvar3dsnow=0
    this%halo=0
    this%ncells=0
    this%nlev=0
    this%nlev_lake=0
    this%nlev_snow=0
    this%nlon_target=0
    this%nlat_target=0

    deallocate(this%index_x)
    deallocate(this%index_y)
    deallocate(this%index_x_nomatch)
    deallocate(this%index_y_nomatch)
    deallocate(this%lake_index_x)
    deallocate(this%lake_index_y)

    deallocate(this%var_mpas)
    deallocate(this%var_rap)
    deallocate(this%var_mpas_lake)
    deallocate(this%var_rap_lake)

  end subroutine close

end module module_surface
