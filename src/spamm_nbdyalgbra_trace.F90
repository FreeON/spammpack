  recursive function absmax_qutree (A) result(absmax)

    type(qutree), pointer, intent(in) :: A
    real(kind(0d0)) :: absmax

    absmax = 0

    if(.not. associated(A)) then
       return
    endif

    if(A%i_upper-A%i_lower+1 == SPAMM_BLOCK_SIZE) then
       absmax = maxval(abs(A%blok))
    else
       if(associated(A%quad11)) then
          absmax = max(absmax, absmax_qutree(A%quad11))
       endif
       if(associated(A%quad12)) then
          absmax = max(absmax, absmax_qutree(A%quad12))
       endif
       if(associated(A%quad21)) then
          absmax = max(absmax, absmax_qutree(A%quad21))
       endif
       if(associated(A%quad22)) then
          absmax = max(absmax, absmax_qutree(A%quad22))
       endif
    endif

  end function absmax_qutree
