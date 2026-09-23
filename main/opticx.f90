
program opticx
  use constants_math
  use parser_input_file
  use parser_wannier90_tb
  use parser_kp
  use parser_optics_xatu_dim
  use exciton_envelopes
  use bands
  use ome_sp
  use ome_ex
  use ome
  use optical_response
  use sigma_first_sp
  use sigma_first_ex
  use sigma_second_sp
  use sigma_second_ex

  
  implicit none
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !read input file for optical calculations
  call get_input_file()
  !get TB/kp Hamiltonian from the wannier90/kp file
  if (inputfile_type == 'Wannier') then
    call wannier90_get(material_name_in)
  else if (inputfile_type == 'kp') then
    call kp_get(material_name_in)
  end if
  !get grid, number of bands and exciton quantities: 
  !either by reading from xatu-output or opticx-input
  call get_optics_xatu_dim()
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call get_energy_bands()
  !pause
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !evaluate single particle optical matrix elements: only VME by now
  call get_ome()
  !evaluate and print requested optical responses
  call get_optical_response()
 
  write(*,*) 'Opticx calculation ended'

end program 