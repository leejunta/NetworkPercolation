!  XDR Fortran Interface
!  2014 (c) James W. Barnett <jbarnet4@tulane.edu>
!  https://github.com/wesbarnett/

module xtc

    use, intrinsic :: iso_c_binding, only: C_PTR, C_CHAR, C_FLOAT, C_INT

    implicit none
    private
    public xdrfile_open,xdrfile_close, read_xtc_natoms, read_xtc

    ! the data type located in libxdrfile
    type, public, bind(C) :: xdrfile
      type(C_PTR) :: fp, xdr
      character(kind=C_CHAR) :: mode
      integer(C_INT) :: buf1, buf1size, buf2, buf2size
    end type xdrfile

    ! interface with libxdrfile
    interface 

      integer(C_INT) function read_xtc_natoms(filename,NATOMS) bind(C, name='read_xtc_natoms')
        import
        character(kind=C_CHAR), intent(in) :: filename
        integer(C_INT), intent(out) :: NATOMS
      end function
      
      type(C_PTR) function xdrfile_open(filename,mode) bind(C, name='xdrfile_open')
        import
        character(kind=C_CHAR), intent(in) :: filename(*), mode(*)
      end function

      integer(C_INT) function read_xtc(xd,NATOMS,STEP,time,box,x,prec) bind(C, name='read_xtc')
        import
        type(xdrfile), intent(in) :: xd
        integer(C_INT), intent(out) :: NATOMS, STEP
        real(C_FLOAT), intent(out) :: time, prec, box(*), x(*)
      end function

      integer(C_INT) function xdrfile_close(xd) bind(C,name='xdrfile_close')
        import
        type(xdrfile), intent(in) :: xd
      end function

    end interface

end module xtc