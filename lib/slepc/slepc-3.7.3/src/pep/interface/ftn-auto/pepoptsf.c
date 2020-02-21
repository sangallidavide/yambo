#include "petscsys.h"
#include "petscfix.h"
#include "petsc/private/fortranimpl.h"
/* pepopts.c */
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

#include "slepcpep.h"
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepsetfromoptions_ PEPSETFROMOPTIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepsetfromoptions_ pepsetfromoptions
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepgettolerances_ PEPGETTOLERANCES
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepgettolerances_ pepgettolerances
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepsettolerances_ PEPSETTOLERANCES
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepsettolerances_ pepsettolerances
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepgetdimensions_ PEPGETDIMENSIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepgetdimensions_ pepgetdimensions
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepsetdimensions_ PEPSETDIMENSIONS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepsetdimensions_ pepsetdimensions
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepsetwhicheigenpairs_ PEPSETWHICHEIGENPAIRS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepsetwhicheigenpairs_ pepsetwhicheigenpairs
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepgetwhicheigenpairs_ PEPGETWHICHEIGENPAIRS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepgetwhicheigenpairs_ pepgetwhicheigenpairs
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepsetproblemtype_ PEPSETPROBLEMTYPE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepsetproblemtype_ pepsetproblemtype
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepgetproblemtype_ PEPGETPROBLEMTYPE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepgetproblemtype_ pepgetproblemtype
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepsetbasis_ PEPSETBASIS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepsetbasis_ pepsetbasis
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepgetbasis_ PEPGETBASIS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepgetbasis_ pepgetbasis
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepsettrackall_ PEPSETTRACKALL
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepsettrackall_ pepsettrackall
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepgettrackall_ PEPGETTRACKALL
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepgettrackall_ pepgettrackall
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepsetconvergencetest_ PEPSETCONVERGENCETEST
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepsetconvergencetest_ pepsetconvergencetest
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepgetconvergencetest_ PEPGETCONVERGENCETEST
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepgetconvergencetest_ pepgetconvergencetest
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepsetstoppingtest_ PEPSETSTOPPINGTEST
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepsetstoppingtest_ pepsetstoppingtest
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepgetstoppingtest_ PEPGETSTOPPINGTEST
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepgetstoppingtest_ pepgetstoppingtest
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepsetscale_ PEPSETSCALE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepsetscale_ pepsetscale
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepgetscale_ PEPGETSCALE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepgetscale_ pepgetscale
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepsetextract_ PEPSETEXTRACT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepsetextract_ pepsetextract
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepgetextract_ PEPGETEXTRACT
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepgetextract_ pepgetextract
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepsetrefine_ PEPSETREFINE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepsetrefine_ pepsetrefine
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define pepgetrefine_ PEPGETREFINE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define pepgetrefine_ pepgetrefine
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
PETSC_EXTERN void PETSC_STDCALL  pepsetfromoptions_(PEP pep, int *__ierr ){
*__ierr = PEPSetFromOptions(
	(PEP)PetscToPointer((pep) ));
}
PETSC_EXTERN void PETSC_STDCALL  pepgettolerances_(PEP pep,PetscReal *tol,PetscInt *maxits, int *__ierr ){
*__ierr = PEPGetTolerances(
	(PEP)PetscToPointer((pep) ),tol,maxits);
}
PETSC_EXTERN void PETSC_STDCALL  pepsettolerances_(PEP pep,PetscReal *tol,PetscInt *maxits, int *__ierr ){
*__ierr = PEPSetTolerances(
	(PEP)PetscToPointer((pep) ),*tol,*maxits);
}
PETSC_EXTERN void PETSC_STDCALL  pepgetdimensions_(PEP pep,PetscInt *nev,PetscInt *ncv,PetscInt *mpd, int *__ierr ){
*__ierr = PEPGetDimensions(
	(PEP)PetscToPointer((pep) ),nev,ncv,mpd);
}
PETSC_EXTERN void PETSC_STDCALL  pepsetdimensions_(PEP pep,PetscInt *nev,PetscInt *ncv,PetscInt *mpd, int *__ierr ){
*__ierr = PEPSetDimensions(
	(PEP)PetscToPointer((pep) ),*nev,*ncv,*mpd);
}
PETSC_EXTERN void PETSC_STDCALL  pepsetwhicheigenpairs_(PEP pep,PEPWhich *which, int *__ierr ){
*__ierr = PEPSetWhichEigenpairs(
	(PEP)PetscToPointer((pep) ),*which);
}
PETSC_EXTERN void PETSC_STDCALL  pepgetwhicheigenpairs_(PEP pep,PEPWhich *which, int *__ierr ){
*__ierr = PEPGetWhichEigenpairs(
	(PEP)PetscToPointer((pep) ),which);
}
PETSC_EXTERN void PETSC_STDCALL  pepsetproblemtype_(PEP pep,PEPProblemType *type, int *__ierr ){
*__ierr = PEPSetProblemType(
	(PEP)PetscToPointer((pep) ),*type);
}
PETSC_EXTERN void PETSC_STDCALL  pepgetproblemtype_(PEP pep,PEPProblemType *type, int *__ierr ){
*__ierr = PEPGetProblemType(
	(PEP)PetscToPointer((pep) ),type);
}
PETSC_EXTERN void PETSC_STDCALL  pepsetbasis_(PEP pep,PEPBasis *basis, int *__ierr ){
*__ierr = PEPSetBasis(
	(PEP)PetscToPointer((pep) ),*basis);
}
PETSC_EXTERN void PETSC_STDCALL  pepgetbasis_(PEP pep,PEPBasis *basis, int *__ierr ){
*__ierr = PEPGetBasis(
	(PEP)PetscToPointer((pep) ),basis);
}
PETSC_EXTERN void PETSC_STDCALL  pepsettrackall_(PEP pep,PetscBool *trackall, int *__ierr ){
*__ierr = PEPSetTrackAll(
	(PEP)PetscToPointer((pep) ),*trackall);
}
PETSC_EXTERN void PETSC_STDCALL  pepgettrackall_(PEP pep,PetscBool *trackall, int *__ierr ){
*__ierr = PEPGetTrackAll(
	(PEP)PetscToPointer((pep) ),trackall);
}
PETSC_EXTERN void PETSC_STDCALL  pepsetconvergencetest_(PEP pep,PEPConv *conv, int *__ierr ){
*__ierr = PEPSetConvergenceTest(
	(PEP)PetscToPointer((pep) ),*conv);
}
PETSC_EXTERN void PETSC_STDCALL  pepgetconvergencetest_(PEP pep,PEPConv *conv, int *__ierr ){
*__ierr = PEPGetConvergenceTest(
	(PEP)PetscToPointer((pep) ),conv);
}
PETSC_EXTERN void PETSC_STDCALL  pepsetstoppingtest_(PEP pep,PEPStop *stop, int *__ierr ){
*__ierr = PEPSetStoppingTest(
	(PEP)PetscToPointer((pep) ),*stop);
}
PETSC_EXTERN void PETSC_STDCALL  pepgetstoppingtest_(PEP pep,PEPStop *stop, int *__ierr ){
*__ierr = PEPGetStoppingTest(
	(PEP)PetscToPointer((pep) ),stop);
}
PETSC_EXTERN void PETSC_STDCALL  pepsetscale_(PEP pep,PEPScale *scale,PetscReal *alpha,Vec Dl,Vec Dr,PetscInt *its,PetscReal *lambda, int *__ierr ){
*__ierr = PEPSetScale(
	(PEP)PetscToPointer((pep) ),*scale,*alpha,
	(Vec)PetscToPointer((Dl) ),
	(Vec)PetscToPointer((Dr) ),*its,*lambda);
}
PETSC_EXTERN void PETSC_STDCALL  pepgetscale_(PEP pep,PEPScale *scale,PetscReal *alpha,Vec *Dl,Vec *Dr,PetscInt *its,PetscReal *lambda, int *__ierr ){
*__ierr = PEPGetScale(
	(PEP)PetscToPointer((pep) ),scale,alpha,Dl,Dr,its,lambda);
}
PETSC_EXTERN void PETSC_STDCALL  pepsetextract_(PEP pep,PEPExtract *extract, int *__ierr ){
*__ierr = PEPSetExtract(
	(PEP)PetscToPointer((pep) ),*extract);
}
PETSC_EXTERN void PETSC_STDCALL  pepgetextract_(PEP pep,PEPExtract *extract, int *__ierr ){
*__ierr = PEPGetExtract(
	(PEP)PetscToPointer((pep) ),extract);
}
PETSC_EXTERN void PETSC_STDCALL  pepsetrefine_(PEP pep,PEPRefine *refine,PetscInt *npart,PetscReal *tol,PetscInt *its,PEPRefineScheme *scheme, int *__ierr ){
*__ierr = PEPSetRefine(
	(PEP)PetscToPointer((pep) ),*refine,*npart,*tol,*its,*scheme);
}
PETSC_EXTERN void PETSC_STDCALL  pepgetrefine_(PEP pep,PEPRefine *refine,PetscInt *npart,PetscReal *tol,PetscInt *its,PEPRefineScheme *scheme, int *__ierr ){
*__ierr = PEPGetRefine(
	(PEP)PetscToPointer((pep) ),refine,npart,tol,its,scheme);
}
#if defined(__cplusplus)
}
#endif
