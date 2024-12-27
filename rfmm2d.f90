MODULE glue
  USE, INTRINSIC :: ISO_C_Binding
  IMPLICIT NONE
CONTAINS
  SUBROUTINE rfmm2d(ierr, iprec, nsource, source, dipstr, dipvec, ntarget, target, pottarg) BIND (C, NAME = 'rfmm2d')
    IMPLICIT NONE
    
    ! Déclarations des arguments
    INTEGER(C_INT),  INTENT(OUT)  :: ierr
    INTEGER(C_INT),  INTENT(IN)   :: iprec
    INTEGER(C_INT),  INTENT(IN)   :: nsource
    REAL(C_DOUBLE),   dimension(2, nsource), INTENT(IN)  :: source
    REAL(C_DOUBLE),   dimension(nsource),    INTENT(IN)  :: dipstr
    REAL(C_DOUBLE),   dimension(2, nsource), INTENT(IN)  :: dipvec
    INTEGER(C_INT),  INTENT(IN)  :: ntarget
    REAL(C_DOUBLE),  dimension(ntarget),  INTENT(IN)  :: target
    REAL(C_DOUBLE),  dimension(ntarget),  INTENT(OUT) :: pottarg
    
    ! Déclarations explicites des variables locales
    INTEGER(C_INT) :: ifcharge, ifdipole, ifpot, ifgrad, ifhess
    INTEGER(C_INT) :: ifpottarg, ifgradtarg, ifhesstarg
    REAL(C_DOUBLE), DIMENSION(:), ALLOCATABLE :: charge, pot, grad, hess, gradtarg, hesstarg

    ! Initialisation des variables locales
    ifcharge = 0
    ifdipole = 1
    ifpot = 0
    ifgrad = 0
    ifhess = 0
    ifpottarg = 1
    ifgradtarg = 0
    ifhesstarg = 0
    
    ! Appel de la routine externe
    CALL rfmm2dparttarg(ierr, iprec, nsource, source, &
         ifcharge, charge, ifdipole, dipstr, dipvec, &
         ifpot, pot, ifgrad, grad, ifhess, hess, &
         ntarget, target, ifpottarg, pottarg, ifgradtarg, gradtarg, &
         ifhesstarg, hesstarg)
         
  END SUBROUTINE rfmm2d

END MODULE glue

