! Fatal error handler

subroutine fatal(msg)
    character(len=*), intent(in) :: msg

    print *, "ERROR: " // msg
    stop
end subroutine fatal
