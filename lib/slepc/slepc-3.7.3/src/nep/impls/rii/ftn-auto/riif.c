#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* rii.c */
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
#define nepriisetmaximumiterations_ NEPRIISETMAXIMUMITERATIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepriisetmaximumiterations_ nepriisetmaximumiterations
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepriigetmaximumiterations_ NEPRIIGETMAXIMUMITERATIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepriigetmaximumiterations_ nepriigetmaximumiterations
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepriisetlagpreconditioner_ NEPRIISETLAGPRECONDITIONER
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepriisetlagpreconditioner_ nepriisetlagpreconditioner
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepriigetlagpreconditioner_ NEPRIIGETLAGPRECONDITIONER
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepriigetlagpreconditioner_ nepriigetlagpreconditioner
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepriisetconstcorrectiontol_ NEPRIISETCONSTCORRECTIONTOL
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepriisetconstcorrectiontol_ nepriisetconstcorrectiontol
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepriigetconstcorrectiontol_ NEPRIIGETCONSTCORRECTIONTOL
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepriigetconstcorrectiontol_ nepriigetconstcorrectiontol
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepriisetksp_ NEPRIISETKSP
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepriisetksp_ nepriisetksp
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define nepriigetksp_ NEPRIIGETKSP
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define nepriigetksp_ nepriigetksp
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
PETSC_EXTERN void PETSC_STDCALL  nepriisetmaximumiterations_(NEP nep,PetscInt *its, int *__ierr ){
*__ierr = NEPRIISetMaximumIterations(
	(NEP)PetscToPointer((nep) ),*its);
}
PETSC_EXTERN void PETSC_STDCALL  nepriigetmaximumiterations_(NEP nep,PetscInt *its, int *__ierr ){
*__ierr = NEPRIIGetMaximumIterations(
	(NEP)PetscToPointer((nep) ),its);
}
PETSC_EXTERN void PETSC_STDCALL  nepriisetlagpreconditioner_(NEP nep,PetscInt *lag, int *__ierr ){
*__ierr = NEPRIISetLagPreconditioner(
	(NEP)PetscToPointer((nep) ),*lag);
}
PETSC_EXTERN void PETSC_STDCALL  nepriigetlagpreconditioner_(NEP nep,PetscInt *lag, int *__ierr ){
*__ierr = NEPRIIGetLagPreconditioner(
	(NEP)PetscToPointer((nep) ),lag);
}
PETSC_EXTERN void PETSC_STDCALL  nepriisetconstcorrectiontol_(NEP nep,PetscBool *cct, int *__ierr ){
*__ierr = NEPRIISetConstCorrectionTol(
	(NEP)PetscToPointer((nep) ),*cct);
}
PETSC_EXTERN void PETSC_STDCALL  nepriigetconstcorrectiontol_(NEP nep,PetscBool *cct, int *__ierr ){
*__ierr = NEPRIIGetConstCorrectionTol(
	(NEP)PetscToPointer((nep) ),cct);
}
PETSC_EXTERN void PETSC_STDCALL  nepriisetksp_(NEP nep,KSP ksp, int *__ierr ){
*__ierr = NEPRIISetKSP(
	(NEP)PetscToPointer((nep) ),
	(KSP)PetscToPointer((ksp) ));
}
PETSC_EXTERN void PETSC_STDCALL  nepriigetksp_(NEP nep,KSP *ksp, int *__ierr ){
*__ierr = NEPRIIGetKSP(
	(NEP)PetscToPointer((nep) ),ksp);
}
#if defined(__cplusplus)
}
#endif
