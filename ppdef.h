
! fortran unit of error messages

#define stderr 0

! fortran unit of screen output

#define stdout 6

! initfile global pattern matching (c function gmatch in library -lgmatch -lgen )

!#define GMATCH

#define procnum 1
#define ALLOCATE_LOCAL_VARS

#define HAS_GETPID

! the exit subroutine is not conform with the Fortran 95 standard
! but it is very usefull for shell scripts

#define ERROR_STOP call exit(1)

! if strict Fortran 95 standard is required use, ERROR_STOP should be
! defined as

!#define ERROR_STOP stop


! uncomment this line for MISP compiler
! For g95, pgf90 and ifort flush can be called without a status variable

#define flush(unit,stat) FLUSH(unit)

 
! quick and dirty hack for double precision. FIX ME
#ifdef DOUBLE_PRECISION
#  define ssaupd dsaupd
#  define sseupd dseupd
#  define snrm2 dnrm2
#  define saxpy daxpy
#  define sscal dscal
#  define scopy dcopy
#  define sgemv dgemv
#endif

! enable NetCDF input and output
#define NETCDF

! enable diagnostics
!#define DIAG_CROSSVAL

!filter B by diffusion
#define B_DIFF

!#define ALL_POINTS_CROSS_VALIDATION
