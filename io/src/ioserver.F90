!> The main IO server functionality which handles waiting for commands and data both of which are delt with.
!! The lower level details of the communication, configuration parsing etc are all held elsewhere. The server
!! can be thought of similar to a bus, with command and data channels. The command gives context to what is on
!! the data channel and not all commands require data (such as deregistration of MONC process)
module io_server_mod
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH
  use configuration_parser_mod, only : DATA_SIZE_STRIDE, io_configuration_type, io_configuration_data_definition_type, &
       io_configuration_registered_monc_type, configuration_parse, extend_registered_moncs_array, retrieve_data_definition, &
       build_definition_description_type_from_configuration, build_field_description_type_from_configuration
  use mpi_communication_mod, only : build_mpi_datatype, data_receive, test_for_command, register_command_receive, &
       cancel_requests, free_mpi_type, get_number_io_servers, get_my_io_rank, test_for_inter_io
  use diagnostic_federator_mod, only : initialise_diagnostic_federator, finalise_diagnostic_federator, &
       check_diagnostic_federator_for_completion, pass_fields_to_diagnostics_federator
  use writer_federator_mod, only : initialise_writer_federator, finalise_writer_federator, check_writer_for_trigger
  use writer_field_manager_mod, only : initialise_writer_field_manager, finalise_writer_field_manager, &
       provide_monc_data_to_writer_federator
  use collections_mod, only : map_type, c_get, c_put, c_is_empty, c_remove, c_add, c_size, c_value_at
  use conversions_mod, only : conv_to_integer, conv_to_generic, conv_to_string
  use io_server_client_mod, only : REGISTER_COMMAND, DEREGISTER_COMMAND, INTER_IO_COMMUNICATION, DATA_COMMAND_START, DATA_TAG, &
       LOCAL_SIZES_KEY, LOCAL_START_POINTS_KEY, LOCAL_END_POINTS_KEY, data_sizing_description_type, definition_description_type, &
       field_description_type, build_mpi_type_data_sizing_description, get_data_description_from_name, &
       build_mpi_type_field_description, build_mpi_type_definition_description
  use forthread_mod, only : forthread_rwlock_rdlock, forthread_rwlock_wrlock, forthread_rwlock_tryrdlock, &
       forthread_rwlock_unlock, forthread_rwlock_init, forthread_rwlock_destroy, forthread_mutex_init, forthread_mutex_lock, &
       forthread_mutex_unlock, forthread_cond_wait, forthread_cond_signal, forthread_cond_init
  use threadpool_mod, only : threadpool_init, threadpool_finalise, threadpool_start_thread, check_thread_status, &
       threadpool_deactivate, threadpool_is_idle
  use global_callback_inter_io_mod, only : perform_global_callback
  use logging_mod, only : LOG_ERROR, LOG_WARN, log_log, initialise_logging
  use mpi, only : MPI_COMM_WORLD, MPI_STATUSES_IGNORE
  implicit none

#ifndef TEST_MODE
  private
#endif  

  integer :: mpi_type_data_sizing_description, & !< The MPI type for field sizing (i.e. array size etc send when MONCs register)
       mpi_type_definition_description, & !< The MPI data type for data descriptions sent to MONCs
       mpi_type_field_description !< The MPI data type for field descriptions sent to MONCs
  type(io_configuration_type), volatile, save :: io_configuration !< Internal representation of the IO configuration
  logical :: contine_poll_messages !< Whether to continue waiting command messages from any MONC processes
  logical, volatile :: contine_poll_interio_messages, already_registered_finishing_call
  type(field_description_type), dimension(:), allocatable :: registree_field_descriptions
  type(definition_description_type), dimension(:), allocatable :: registree_definition_descriptions

  integer :: monc_registration_lock

  public io_server_run
contains

  !> Called to start the IO server and once this subroutine returns then it indicates that the IO server has finished.
  !! The runtime is spent in here awaiting commands and then dealing with them. Termination occurs when all MONC processes
  !! have deregistered, note that to trigger this then at least one MONC process must first register
  !! @param io_communicator_arg The IO communicator containing just the IO servers
  !! @param io_xml_configuration Textual XML configuration that is used to set up the IO server
  subroutine io_server_run(options_database, io_communicator_arg, io_xml_configuration, &
       provided_threading, total_global_processes)
    type(map_type), intent(inout) :: options_database
    integer, intent(in) :: io_communicator_arg, provided_threading, total_global_processes
    character, dimension(:), allocatable, intent(inout) :: io_xml_configuration

    integer :: recv_type, i, command, source, ierr
    character, dimension(:), allocatable :: data_buffer
    logical :: message_waiting

    call configuration_parse(options_database, io_xml_configuration, io_configuration)
    deallocate(io_xml_configuration)
    call threadpool_init(io_configuration, provided_threading)
    call check_thread_status(forthread_rwlock_init(monc_registration_lock, -1))
    contine_poll_messages=.true.
    contine_poll_interio_messages=.true.
    already_registered_finishing_call=.false.
    io_configuration%io_communicator=io_communicator_arg
    io_configuration%number_of_io_servers=get_number_io_servers(io_communicator_arg)
    io_configuration%number_of_global_moncs=total_global_processes-io_configuration%number_of_io_servers
    io_configuration%my_io_rank=get_my_io_rank(io_communicator_arg)
    call initialise_logging(io_configuration%my_io_rank)
    registree_definition_descriptions=build_definition_description_type_from_configuration(io_configuration)
    registree_field_descriptions=build_field_description_type_from_configuration(io_configuration)
    call initialise_diagnostic_federator(io_configuration)
    call initialise_writer_federator(io_configuration)
    call initialise_writer_field_manager(io_configuration)

    mpi_type_data_sizing_description=build_mpi_type_data_sizing_description()
    mpi_type_definition_description=build_mpi_type_definition_description()
    mpi_type_field_description=build_mpi_type_field_description()

    call register_command_receive()

    do while (await_command(command, source, data_buffer))      
      call handle_command_message(command, source, data_buffer)
    end do
    call threadpool_deactivate()
    call finalise_writer_field_manager()
    call finalise_writer_federator()
    call finalise_diagnostic_federator(io_configuration)
    call check_thread_status(forthread_rwlock_destroy(monc_registration_lock))
    call threadpool_finalise()
    call free_individual_registered_monc_aspects()
    call cancel_requests()
    call free_mpi_type(mpi_type_data_sizing_description)
    call free_mpi_type(mpi_type_definition_description)
    call free_mpi_type(mpi_type_field_description)    
  end subroutine io_server_run

  !> Awaits a command or shutdown from MONC processes and other IO servers
  !! @param command The command received is output
  !! @param source The source process received is output
  !! @returns Whether to continue polling for commands (and whether to process the current output)
  logical function await_command(command, source, data_buffer)
    integer, intent(out) :: command, source
    character, dimension(:), allocatable :: data_buffer

    logical :: completed, inter_io_complete

    completed=.false.
    await_command=.false.
    do while(.not. completed)
      if (.not. contine_poll_messages .and. .not. contine_poll_interio_messages) return
      if (contine_poll_messages) then
        if (test_for_command(command, source)) then
          await_command=.true.
          return
        end if
      end if
      if (contine_poll_interio_messages .and. allocated(io_configuration%inter_io_communications)) then       
        inter_io_complete=test_for_inter_io(io_configuration%inter_io_communications, &
             io_configuration%number_inter_io_communications, io_configuration%io_communicator, command, source, data_buffer) 
        if (inter_io_complete) then
          await_command=.true.
          return
        end if
      end if
      if (.not. contine_poll_messages .and. .not. already_registered_finishing_call) then
        if (check_diagnostic_federator_for_completion(io_configuration) .and. threadpool_is_idle()) then
          already_registered_finishing_call=.true.          
          call perform_global_callback(io_configuration, "termination", 1, termination_callback)          
        end if
      end if      
    end do    
  end function await_command

  !> This is the termination callback which is called once all MONCs have deregistered, no sends are active by inter IO
  !! communications and all threads are idle. This shuts down the inter IO listening and kickstarts finalisation and closure
  !! @param io_configuration The IO server configuration
  !! @param values Values (ignored)
  !! @param field_name Field name identifier
  !! @param timestep Timestep identifier
  subroutine termination_callback(io_configuration, values, field_name, timestep)
    type(io_configuration_type), intent(inout) :: io_configuration
    real(DEFAULT_PRECISION), dimension(:) :: values
    character(len=STRING_LENGTH) :: field_name
    integer :: timestep

    contine_poll_interio_messages=.false.
  end subroutine termination_callback  

  !> Called to handle a specific command that has been recieved
  !! @param command The command which has been received from some process
  !! @param source The PID of the source (MONC) process
  subroutine handle_command_message(command, source, data_buffer)
    integer, intent(in) :: command, source
    character, dimension(:), allocatable, intent(inout) :: data_buffer

    if (command == REGISTER_COMMAND) then
      call threadpool_start_thread(handle_monc_registration, (/ source /))
    else if (command == DEREGISTER_COMMAND) then
      call threadpool_start_thread(handle_deregistration_command, (/ source /))
    else if (command == INTER_IO_COMMUNICATION) then
      call threadpool_start_thread(handle_inter_io_communication_command, (/ source /), data_buffer=data_buffer)      
      deallocate(data_buffer)
    else if (command .ge. DATA_COMMAND_START) then
      call threadpool_start_thread(handle_data_message, (/ source,  command-DATA_COMMAND_START /))
    end if    
  end subroutine handle_command_message

  !> A helper function to get the location of a MONC's configuration in the IO data structure
  !! @param source Source index of the MONC process
  !! @returns Index that that MONC corresponds to
  integer function get_monc_location(source)
    integer, intent(in) :: source

    class(*), pointer :: generic

    generic=>c_get(io_configuration%monc_to_index, conv_to_string(source))
    get_monc_location=conv_to_integer(generic, .false.)
  end function get_monc_location

  !> Handles inter IO server communications
  !! @param arguments The thread based arguments, this is the index of the inter IO server description
  subroutine handle_inter_io_communication_command(arguments, data_buffer)
    integer, dimension(:), intent(in) :: arguments
    character, dimension(:), intent(inout), optional :: data_buffer

    integer :: source

    source=arguments(1)

    call io_configuration%inter_io_communications(source)%handling_procedure(io_configuration, data_buffer, source)
  end subroutine handle_inter_io_communication_command

  !> Frees up the memory associated with individual registered MONCs. This is done at the end for all MONCs as we can't
  !! deallocate dynamically in a threaded environment without excessive ordering and locking in case some data processing
  !! is queued or in progress
  subroutine free_individual_registered_monc_aspects()
    class(*), pointer :: generic
    integer :: i, j, specific_monc_data_type
    
    do i=1, size(io_configuration%registered_moncs)
      do j=1, c_size(io_configuration%registered_moncs(i)%registered_monc_types)
        generic=>c_value_at(io_configuration%registered_moncs(i)%registered_monc_types, j)
        specific_monc_data_type=conv_to_integer(generic, .false.)
        call free_mpi_type(specific_monc_data_type)
      end do
      if (allocated(io_configuration%registered_moncs(i)%field_start_locations)) &
           deallocate(io_configuration%registered_moncs(i)%field_start_locations)
      if (allocated(io_configuration%registered_moncs(i)%field_end_locations)) &
           deallocate(io_configuration%registered_moncs(i)%field_end_locations)
      if (allocated(io_configuration%registered_moncs(i)%definition_names)) &
           deallocate(io_configuration%registered_moncs(i)%definition_names)
      if (allocated(io_configuration%registered_moncs(i)%dimensions)) deallocate(io_configuration%registered_moncs(i)%dimensions)
    end do
  end subroutine free_individual_registered_monc_aspects  

  !> Deregisteres a specific MONC source process
  !! @param source The MONC process PID that we are deregistering
  subroutine handle_deregistration_command(arguments, data_buffer)
    integer, dimension(:), intent(in) :: arguments
    character, dimension(:), intent(inout), optional :: data_buffer

    integer :: monc_location, source

    source=arguments(1)
    monc_location=get_monc_location(source)
    call check_thread_status(forthread_mutex_lock(io_configuration%registered_moncs(monc_location)%active_mutex))
    do while (io_configuration%registered_moncs(monc_location)%active_threads .gt. 0)
      call check_thread_status(forthread_cond_wait(io_configuration%registered_moncs(monc_location)%deactivate_condition_variable,&
             io_configuration%registered_moncs(monc_location)%active_mutex))
    end do
    call check_thread_status(forthread_mutex_unlock(io_configuration%registered_moncs(monc_location)%active_mutex))
    call check_thread_status(forthread_rwlock_wrlock(monc_registration_lock))
    io_configuration%active_moncs=io_configuration%active_moncs-1
    if (io_configuration%active_moncs==0) contine_poll_messages=.false.
    call check_thread_status(forthread_rwlock_unlock(monc_registration_lock))
  end subroutine handle_deregistration_command  

  !> Handles the command for data download from a specific process. This will allocate the receive buffer
  !! and then call to get the data. Once it has been received then the data is run against handling rules
  !! @param data_set The data set ID that is being sent over
  !! @param source The source PID
  subroutine handle_data_message(arguments, input_data_buffer)
    integer, dimension(:), intent(in) :: arguments
    character, dimension(:), intent(inout), optional :: input_data_buffer

    integer :: specific_monc_data_type, specific_monc_buffer_size, recv_count, monc_location, data_set, &
         source, matched_datadefn_index
    character, dimension(:), allocatable :: data_buffer
    class(*), pointer :: generic

    source=arguments(1)
    data_set=arguments(2)
    monc_location=get_monc_location(source) 

    call check_thread_status(forthread_rwlock_rdlock(monc_registration_lock))

    call check_thread_status(forthread_mutex_lock(io_configuration%registered_moncs(monc_location)%active_mutex))
    io_configuration%registered_moncs(monc_location)%active_threads=&
         io_configuration%registered_moncs(monc_location)%active_threads+1
    call check_thread_status(forthread_mutex_unlock(io_configuration%registered_moncs(monc_location)%active_mutex))

    generic=>c_get(io_configuration%registered_moncs(monc_location)%registered_monc_types, conv_to_string(data_set))
    specific_monc_data_type=conv_to_integer(generic, .false.)
    generic=>c_get(io_configuration%registered_moncs(monc_location)%registered_monc_buffer_sizes, conv_to_string(data_set))
    specific_monc_buffer_size=conv_to_integer(generic, .false.)

    allocate(data_buffer(specific_monc_buffer_size))      
    recv_count=data_receive(specific_monc_data_type, 1, source, dump_data=data_buffer, data_dump_id=data_set)
    
    matched_datadefn_index=retrieve_data_definition(io_configuration, &
         io_configuration%registered_moncs(monc_location)%definition_names(data_set))

    if (matched_datadefn_index .gt. 0) then
      call pass_fields_to_diagnostics_federator(io_configuration, source, data_set, data_buffer)
      call provide_monc_data_to_writer_federator(io_configuration, source, data_set, data_buffer)
      call check_writer_for_trigger(io_configuration, source, data_set, data_buffer)
    else
      call log_log(LOG_WARN, "IO server can not find matching data definition with name "&
           //io_configuration%registered_moncs(monc_location)%definition_names(data_set))
    end if    

    call check_thread_status(forthread_mutex_lock(io_configuration%registered_moncs(monc_location)%active_mutex))
    io_configuration%registered_moncs(monc_location)%active_threads=&
         io_configuration%registered_moncs(monc_location)%active_threads-1
    call check_thread_status(forthread_cond_signal(io_configuration%registered_moncs(monc_location)%deactivate_condition_variable))
    call check_thread_status(forthread_mutex_unlock(io_configuration%registered_moncs(monc_location)%active_mutex))
    call check_thread_status(forthread_rwlock_unlock(monc_registration_lock))
    deallocate(data_buffer)
  end subroutine handle_data_message

  !> Handles registration from some MONC process. The source process sends some data description to this IO server which
  !! basically tells the IO server the size of the array datas (which might be different on different processes in the case
  !! of uneven decomposition.) Based upon this a communication (MPI) data type is constructed and the data size in bytes determined
  !! @param source The PID of the MONC process that is registering itself
  subroutine handle_monc_registration(arguments, data_buffer)
    integer, dimension(:), intent(in) :: arguments
    character, dimension(:), intent(inout), optional :: data_buffer

    class(*), pointer :: generic
    integer :: configuration_send_request(2), ierr, number_data_definitions, this_monc_index, source, suspension_mutex

    source=arguments(1)
    configuration_send_request=send_configuration_to_registree(source)
    number_data_definitions=io_configuration%number_of_data_definitions

    call check_thread_status(forthread_rwlock_wrlock(monc_registration_lock))

    io_configuration%number_of_moncs=io_configuration%number_of_moncs+1
    this_monc_index=io_configuration%number_of_moncs
    if (io_configuration%number_of_moncs .gt. size(io_configuration%registered_moncs)) then
      call extend_registered_moncs_array(io_configuration)      
    end if

    io_configuration%active_moncs=io_configuration%active_moncs+1
    call check_thread_status(forthread_rwlock_unlock(monc_registration_lock))

    generic=>conv_to_generic(this_monc_index, .true.)
    call c_put(io_configuration%monc_to_index, conv_to_string(source), generic)

    call check_thread_status(forthread_mutex_init(io_configuration%registered_moncs(this_monc_index)%active_mutex, -1))
    call check_thread_status(forthread_cond_init(&
         io_configuration%registered_moncs(this_monc_index)%deactivate_condition_variable, -1))
    io_configuration%registered_moncs(this_monc_index)%active_threads=0

    allocate(io_configuration%registered_moncs(this_monc_index)%field_start_locations(number_data_definitions), &
         io_configuration%registered_moncs(this_monc_index)%field_end_locations(number_data_definitions), &
         io_configuration%registered_moncs(this_monc_index)%definition_names(number_data_definitions), &
         io_configuration%registered_moncs(this_monc_index)%dimensions(number_data_definitions))

    ! Wait for configuration to have been sent to registree
    call mpi_waitall(2, configuration_send_request, MPI_STATUSES_IGNORE, ierr) 
    call init_data_definition(source, io_configuration%registered_moncs(this_monc_index))
  end subroutine handle_monc_registration

  !> Sends the data and field descriptions to the MONC process that just registered with the IO server
  !! @param source The MPI rank (MPI_COMM_WORLD) of the registree
  !! @returns The nonblocking send request handles which can be waited for completion later (overlap compute and communication)
  function send_configuration_to_registree(source)
    integer, intent(in) :: source
    integer :: send_configuration_to_registree(2)
    
    integer :: ierr, srequest(2)

    call mpi_isend(registree_definition_descriptions, size(registree_definition_descriptions), mpi_type_definition_description, &
         source, DATA_TAG, MPI_COMM_WORLD, srequest(1), ierr)
    call mpi_isend(registree_field_descriptions, size(registree_field_descriptions), mpi_type_field_description, &
         source, DATA_TAG, MPI_COMM_WORLD, srequest(2), ierr)

    send_configuration_to_registree=srequest    
  end function send_configuration_to_registree  

  !> Initialise the sizing of data definitions from a MONC process. The IO server determines, from configuration, the
  !! structure of each data definition but the size of the arrays depends upon the MONC process (due to uneven distribution
  !! of data etc...) This receives the sizing message and then builds the MPI datatype for each data definition that the IO 
  !! server will receive from that specific MONC process. The field sizings are for all fields in every data definition, and
  !! these are applied to each data definition which will simply ignore non matching fields
  !! @param source The source MONC PID
  !! @param monc_defn The corresponding MONC definition data structure
  subroutine init_data_definition(source, monc_defn)
    integer, intent(in) :: source
    type(io_configuration_registered_monc_type), intent(inout) :: monc_defn

    type(io_configuration_data_definition_type) :: matched_datadefn
    type(data_sizing_description_type) :: data_description(io_configuration%number_of_distinct_data_fields+3)
    integer :: created_mpi_type, data_size, recv_count, i
    class(*), pointer :: generic
    logical :: data_defn_found

    recv_count=data_receive(mpi_type_data_sizing_description, io_configuration%number_of_distinct_data_fields+3, &
         source, description_data=data_description)

    call handle_monc_dimension_information(data_description, monc_defn)
    do i=1, io_configuration%number_of_data_definitions
      created_mpi_type=build_mpi_datatype(io_configuration%data_definitions(i), data_description, data_size, &
           monc_defn%field_start_locations(i), monc_defn%field_end_locations(i), monc_defn%dimensions(i))

      generic=>conv_to_generic(created_mpi_type, .true.)
      call c_put(monc_defn%registered_monc_types, conv_to_string(i), generic)
      generic=>conv_to_generic(data_size, .true.)
      call c_put(monc_defn%registered_monc_buffer_sizes, conv_to_string(i), generic)

      monc_defn%definition_names(i)=io_configuration%data_definitions(i)%name
    end do
  end subroutine init_data_definition

  !> Handles the provided local MONC dimension and data layout information
  !! @param data_description The data descriptions sent over from MONC
  !! @param monc_defn The corresponding MONC definition data structure
  subroutine handle_monc_dimension_information(data_description, monc_defn)
    type(io_configuration_registered_monc_type), intent(inout) :: monc_defn
    type(data_sizing_description_type), dimension(:) :: data_description

    type(data_sizing_description_type) :: field_description
    integer :: i
    logical :: field_found

    field_found=get_data_description_from_name(data_description, LOCAL_SIZES_KEY, field_description)
    if (.not. field_found) call log_log(LOG_ERROR, "Malformed MONC registration, no local size information")
    do i=1,3
      monc_defn%local_dim_sizes(i)=field_description%dim_sizes(i)
    end do
    field_found=get_data_description_from_name(data_description, LOCAL_START_POINTS_KEY, field_description)
    if (.not. field_found) call log_log(LOG_ERROR, "Malformed MONC registration, no local start point information")
    do i=1,3
      monc_defn%local_dim_starts(i)=field_description%dim_sizes(i)
    end do
    field_found=get_data_description_from_name(data_description, LOCAL_END_POINTS_KEY, field_description)
    if (.not. field_found) call log_log(LOG_ERROR, "Malformed MONC registration, no local end point information")
    do i=1,3
      monc_defn%local_dim_ends(i)=field_description%dim_sizes(i)
    end do
  end subroutine handle_monc_dimension_information  
end module io_server_mod