module mpic_mod
  ! This is essentially the same as cold tank_experiments, but we also fix the environmental
  ! relative humidity by adjusting the vapour content.
       
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH
  use monc_component_mod, only : component_descriptor_type
  use state_mod, only : model_state_type
  use grids_mod, only : X_INDEX, Y_INDEX, Z_INDEX
  use science_constants_mod, only : pi, r_over_cp
  use optionsdatabase_mod, only :  options_get_real_array, options_get_real, options_get_logical, &
     options_get_array_size, options_get_string_array
  use q_indices_mod, only: get_q_index, standard_q_names
  implicit none
#ifndef TEST_MODE
  private
#endif
 
  integer, parameter :: UNSET_REAL=-999.0
  integer, parameter :: MAXBUBBLES=10 !maximum number of bubbles
  integer, parameter :: MAXQIN=10 !maximum number of bubbles
  real(kind=DEFAULT_PRECISION) :: x_cen, y_cen, z_cen   ! coordinates of tank_experiments centre
  real(kind=DEFAULT_PRECISION) :: rad   ! radial parameters 
  real(kind=DEFAULT_PRECISION) :: zb,zd,zm,zc,rhb,mu
  real(kind=DEFAULT_PRECISION) :: e1,e2,e3,r_smooth_frac,r_edge
  real(kind=DEFAULT_PRECISION) :: b_m=12.5
 
  logical :: l_moist   ! Moist bubbles
  ! Use this value of PI here (slightly different from the scientific constants to match the LEM tank_experiments setup)
  real(kind=DEFAULT_PRECISION), parameter :: my_pi = 4.0_DEFAULT_PRECISION*atan(1.0_DEFAULT_PRECISION)
  ! q variables
  real(kind=DEFAULT_PRECISION), dimension(MAXBUBBLES*MAXQIN) :: bubble_q_values_tmp ! initial mixing ratios of bubble q variables
  real(kind=DEFAULT_PRECISION), allocatable :: bubble_q_values(:,:) ! initial mixing ratios of bubble q variables
  character(len=STRING_LENGTH) :: bubble_q_names(MAXQIN)='unset'  ! names of q variables to initialize
  integer :: iqv,iql ! index for vapour
  integer :: number_of_bubbles ! The number of bubbles
  integer :: nq_bubbles ! Number of q variables to initialize in the bubbles
  public mpic_get_descriptor
contains
  type(component_descriptor_type) function mpic_get_descriptor()
    mpic_get_descriptor%name="mpic"
    mpic_get_descriptor%version=0.2
    mpic_get_descriptor%initialisation=>initialisation_callback
    mpic_get_descriptor%timestep=>timestep_callback
  end function mpic_get_descriptor
  !! Note that this is not the most efficient way to iterate through theta (j heavy), but it is the same as the LEM set up
  !! so directly comparable and probably doesn't matter too much as it is just called onec in the initialisation
  subroutine initialisation_callback(current_state)
    type(model_state_type), intent(inout), target :: current_state
    integer :: i ! loop counter
    ! Read in parameters from options database
    l_moist=options_get_logical(current_state%options_database, "mpic_lmoist") 
    x_cen=options_get_real(current_state%options_database, "mpic_x_cen") 
    y_cen=options_get_real(current_state%options_database, "mpic_y_cen") 
    z_cen=options_get_real(current_state%options_database, "mpic_z_cen") 
    rad=options_get_real(current_state%options_database, "mpic_rad") 
    rhb=options_get_real(current_state%options_database, "mpic_rhb") 
    zc=options_get_real(current_state%options_database, "mpic_zc") 
    mu=options_get_real(current_state%options_database, "mpic_mu") 
    zd=options_get_real(current_state%options_database, "mpic_zd") 
    zm=options_get_real(current_state%options_database, "mpic_zm") 
    e1=options_get_real(current_state%options_database, "mpic_e1") 
    e2=options_get_real(current_state%options_database, "mpic_e2") 
    e3=options_get_real(current_state%options_database, "mpic_e3") 
    r_smooth_frac=options_get_real(current_state%options_database, "mpic_r_smooth_frac") 
    nq_bubbles=options_get_array_size(current_state%options_database, "bubble_q_names")
    if (nq_bubbles >  0)then
      call options_get_string_array(current_state%options_database, "bubble_q_names", bubble_q_names)
      call options_get_real_array(current_state%options_database, "bubble_q_values", bubble_q_values_tmp)
      allocate(bubble_q_values(nq_bubbles, number_of_bubbles))
      bubble_q_values=reshape(bubble_q_values_tmp, (/ nq_bubbles, number_of_bubbles/))
    end if
   
    ! Add in the q variable index for vapour
    if (l_moist) iqv=get_q_index(standard_q_names%VAPOUR, 'mpic')
    if (l_moist) iql=get_q_index(standard_q_names%CLOUD_LIQUID_MASS, 'mpic')
    call generate_bubbles(current_state)
  end subroutine initialisation_callback
  subroutine generate_bubbles(current_state)
    type(model_state_type), intent(inout), target :: current_state
    real(kind=DEFAULT_PRECISION) :: x, y              ! x and y coordinate
    real(kind=DEFAULT_PRECISION) :: zb,hpl,hen,dbdz,radsq,bpl      ! auxiliary vars
    real(kind=DEFAULT_PRECISION) :: delx, dely, delz  ! useful things
   
    integer :: i,j,k,iq,n ! loop counters
    hpl=exp(-zc)
    hen=mu*hpl
    zb=log(rhb/hen)
    dbdz=b_m*(hpl-exp(-zm))/(zm-zd)
    bpl=dbdz*(zd-zb)
    radsq=rad**2
    e1=e1/radsq
    e2=e2/radsq
    e3=e3/radsq
   
    do i=current_state%local_grid%local_domain_start_index(X_INDEX), current_state%local_grid%local_domain_end_index(X_INDEX)
      x = (current_state%local_grid%start(X_INDEX) + (i-current_state%local_grid%local_domain_start_index(X_INDEX))) * &
         current_state%global_grid%configuration%horizontal%dx
      do k=current_state%local_grid%local_domain_start_index(Z_INDEX), current_state%local_grid%local_domain_end_index(Z_INDEX)
        do j=current_state%local_grid%local_domain_start_index(Y_INDEX), current_state%local_grid%local_domain_end_index(Y_INDEX)         
          y = (current_state%local_grid%start(Y_INDEX) + (j-current_state%local_grid%local_domain_start_index(Y_INDEX))) * &
               current_state%global_grid%configuration%horizontal%dy       
         
          delx = ( x - x_cen)
          delz = (current_state%global_grid%configuration%vertical%zn(k) - z_cen)   
          dely = ( y - y_cen)
          if ((delx**2+dely**2+delz**2) .gt. radsq) then
            !Outside plume:
            if (current_state%global_grid%configuration%vertical%zn(k) .lt. zb) then
              !Mixed layer:
              current_state%th%data(k,j,i)=0.0_DEFAULT_PRECISION
              if (l_moist) current_state%q(iqv)%data(k,j,i)=hen
            else
              !Upper stratified zone:
              current_state%th%data(k,j,i)=dbdz*(current_state%global_grid%configuration%vertical%zn(k)-zb)
              if (l_moist) current_state%q(iqv)%data(k,j,i)=rhb*exp(-current_state%global_grid%configuration%vertical%zn(k))
            endif
          elseif ((delx*delx+dely*dely+delz*delz) <= radsq*r_smooth_frac*r_smooth_frac) then
            !Inside plume:
            current_state%th%data(k,j,i)=bpl*(1.0_DEFAULT_PRECISION+e1*delx*dely+e2*delx*delz+e3*dely*delz)
            if(l_moist) current_state%q(iqv)%data(k,j,i)=hpl
          else 
            !Inside plume edge
            r_edge=(sqrt(delx*delx+dely*dely+delz*delz)-rad*r_smooth_frac)/(rad*(1.0_DEFAULT_PRECISION-r_smooth_frac))
            current_state%th%data(k,j,i)=bpl*(1.0_DEFAULT_PRECISION+e1*delx*dely+e2*delx*delz+e3*dely*delz)*&
            (1.0_DEFAULT_PRECISION - (6.0_DEFAULT_PRECISION*r_edge**5-15.0_DEFAULT_PRECISION*r_edge**4+&
            10.0_DEFAULT_PRECISION*r_edge**3))
            if(l_moist) then 
              current_state%q(iqv)%data(k,j,i)=hpl*(1.0_DEFAULT_PRECISION+(mu-1.0_DEFAULT_PRECISION)*&
              (6.0_DEFAULT_PRECISION*r_edge**5-15.0_DEFAULT_PRECISION*r_edge**4+&
              10.0_DEFAULT_PRECISION*r_edge**3))
            end if
          endif
        end do
      end do
    end do
    if (nq_bubbles >  0)deallocate(bubble_q_values)
  end subroutine generate_bubbles
  subroutine timestep_callback(current_state)
    type(model_state_type), intent(inout), target :: current_state
    integer :: i,j,k ! look counters
    do i=current_state%local_grid%local_domain_start_index(X_INDEX), current_state%local_grid%local_domain_end_index(X_INDEX)
      do k=current_state%local_grid%local_domain_start_index(Z_INDEX), current_state%local_grid%local_domain_end_index(Z_INDEX)
        do j=current_state%local_grid%local_domain_start_index(Y_INDEX), current_state%local_grid%local_domain_end_index(Y_INDEX)         
           if (l_moist) current_state%q(iql)%data(k,j,i)=max(current_state%q(iqv)%data(k,j,i)-&
           exp(-current_state%global_grid%configuration%vertical%zn(k)),0.0)
        end do
      end do
    end do
  end subroutine timestep_callback
end module mpic_mod
