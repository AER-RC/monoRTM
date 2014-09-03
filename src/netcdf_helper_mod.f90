#ifdef USENETCDF
module netcdf_helper_mod
! This interface is used to obtain the netcdf-xtype
! It returns NF90_FLOAT for real*4 and NF90_DOUBLE for real*8
    interface nf90_xtype
      module procedure nf90_xtype_float, nf90_xtype_double
    end interface

contains

  function nf90_xtype_float(var)
    USE netcdf
    integer*4 :: nf90_xtype_float
    real*4, intent (in) :: var
    nf90_xtype_float = NF90_FLOAT
    end function nf90_xtype_float

  function nf90_xtype_double(var)
    USE netcdf
    integer*4 :: nf90_xtype_double
    real*8, intent (in) :: var
    nf90_xtype_double = NF90_DOUBLE
    end function nf90_xtype_double

end module netcdf_helper_mod

#endif

