require 'mkmf'
dir_config("vtc")
have_library('vtc', 'Nbody')
have_header('vtc.h')
have_header('vtclocal.h')

create_makefile("Nbody")
