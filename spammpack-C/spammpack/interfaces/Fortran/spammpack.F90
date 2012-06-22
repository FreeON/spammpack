module spammpack

  use spamm_derived
  use spamm_convert

  interface spamm_version
    function spamm_version ()
    end function spamm_version
  end interface spamm_version

  interface spamm_exit
    subroutine spamm_exit (exitcode)
      integer :: exitcode
    end subroutine spamm_exit
  end interface spamm_exit

end module spammpack
