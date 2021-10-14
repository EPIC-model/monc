!> Dummy stub when not compiling in socrates
module socrates_couple_mod
  use monc_component_mod, only : component_descriptor_type
  use def_socrates_options, only: str_socrates_options
  implicit none

  type (str_socrates_options) :: socrates_opt ! need this to compile DEPHY, even without radiation

  public socrates_couple_get_descriptor, socrates_opt
contains

  !> Provides the descriptor back to the caller and is used in component registration
  !! @returns The termination check component descriptor
  type(component_descriptor_type) function socrates_couple_get_descriptor()
    socrates_couple_get_descriptor%name="socrates_couple"
    socrates_couple_get_descriptor%version=0.1

  end function socrates_couple_get_descriptor

end module socrates_couple_mod






