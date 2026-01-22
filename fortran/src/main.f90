program main
use sim_m
use udp_m
implicit none
character(100) :: filename

call udp_initialize()
call get_command_argument(1,filename)
call init(filename)
call run()

call udp_finalize()
end program main




