#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* nleigs.c */
/* Fortran interface file */

/*
* This file was generated automatically by bfort from the C source
* file.  
 */

#ifdef PETSC_USE_POINTER_CONVERSION
#if defined(__cplusplus)
extern "C" { 
#endif 
extern void *PetscToPointer(void*);
extern int PetscFromPointer(void *);
extern void PetscRmPointer(void*);
#if defined(__cplusplus)
} 
#endif 

#else

#define PetscToPointer(a) (*(PetscFortranAddr *)(a))
#define PetscFromPointer(a) (PetscFortranAddr)(a)
#define PetscRmPointer(a)
#endif

#include "slepcnep.h"
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepnleigssetrestart_ NEPNLEIGSSETRESTART
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepnleigssetrestart_ nepnleigssetrestart
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepnleigsgetrestart_ NEPNLEIGSGETRESTART
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepnleigsgetrestart_ nepnleigsgetrestart
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepnleigssetlocking_ NEPNLEIGSSETLOCKING
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepnleigssetlocking_ nepnleigssetlocking
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepnleigsgetlocking_ NEPNLEIGSGETLOCKING
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepnleigsgetlocking_ nepnleigsgetlocking
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepnleigssetinterpolation_ NEPNLEIGSSETINTERPOLATION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepnleigssetinterpolation_ nepnleigssetinterpolation
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepnleigsgetinterpolation_ NEPNLEIGSGETINTERPOLATION
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepnleigsgetinterpolation_ nepnleigsgetinterpolation
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepnleigssettrueresidual_ NEPNLEIGSSETTRUERESIDUAL
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepnleigssettrueresidual_ nepnleigssettrueresidual
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepnleigsgettrueresidual_ NEPNLEIGSGETTRUERESIDUAL
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepnleigsgettrueresidual_ nepnleigsgettrueresidual
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
PETSC_EXTERN void PETSC_STDCALL  nepnleigssetrestart_(NEP nep,PetscReal *keep, int *__ierr ){
*__ierr = NEPNLEIGSSetRestart(
	(NEP)PetscToPointer((nep) ),*keep);
}
PETSC_EXTERN void PETSC_STDCALL  nepnleigsgetrestart_(NEP nep,PetscReal *keep, int *__ierr ){
*__ierr = NEPNLEIGSGetRestart(
	(NEP)PetscToPointer((nep) ),keep);
}
PETSC_EXTERN void PETSC_STDCALL  nepnleigssetlocking_(NEP nep,PetscBool *lock, int *__ierr ){
*__ierr = NEPNLEIGSSetLocking(
	(NEP)PetscToPointer((nep) ),*lock);
}
PETSC_EXTERN void PETSC_STDCALL  nepnleigsgetlocking_(NEP nep,PetscBool *lock, int *__ierr ){
*__ierr = NEPNLEIGSGetLocking(
	(NEP)PetscToPointer((nep) ),lock);
}
PETSC_EXTERN void PETSC_STDCALL  nepnleigssetinterpolation_(NEP nep,PetscReal *tol,PetscInt *maxits, int *__ierr ){
*__ierr = NEPNLEIGSSetInterpolation(
	(NEP)PetscToPointer((nep) ),*tol,*maxits);
}
PETSC_EXTERN void PETSC_STDCALL  nepnleigsgetinterpolation_(NEP nep,PetscReal *tol,PetscInt *maxits, int *__ierr ){
*__ierr = NEPNLEIGSGetInterpolation(
	(NEP)PetscToPointer((nep) ),tol,maxits);
}
PETSC_EXTERN void PETSC_STDCALL  nepnleigssettrueresidual_(NEP nep,PetscBool *trueres, int *__ierr ){
*__ierr = NEPNLEIGSSetTrueResidual(
	(NEP)PetscToPointer((nep) ),*trueres);
}
PETSC_EXTERN void PETSC_STDCALL  nepnleigsgettrueresidual_(NEP nep,PetscBool *trueres, int *__ierr ){
*__ierr = NEPNLEIGSGetTrueResidual(
	(NEP)PetscToPointer((nep) ),trueres);
}
#if defined(__cplusplus)
}
#endif
