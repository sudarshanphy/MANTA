#include "config.h"

#define MASTER_PROC 0

#define MAXDIM 3
#define IAXIS 1
#define JAXIS 2
#define KAXIS 3

#define LEFT 121
#define RIGHT 122

#define FLOW 221
#define REFLECT 222
#define PERIODIC 223
#define CUSTOM 224

/* Primitive variables */
#define NVAR_BEGIN 1

#define DENS_VAR 1
#define VELX_VAR 2
#define VELY_VAR 3
#define VELZ_VAR 4
#define PRES_VAR 5

#ifdef MHD
#define BMFX_VAR 6
#define BMFY_VAR 7
#define BMFZ_VAR 8
#define BPSI_VAR 9
#define ENER_VAR 10
#else
#define ENER_VAR 6
#endif

#define NVAR_END ENER_VAR
#define NVAR_NUMBER (NVAR_END - NVAR_BEGIN + 1)

/* Conserved variables */
#define NCONSVAR_BEGIN 1

#define DENS_CONS 1
#define MOMX_CONS 2
#define MOMY_CONS 3
#define MOMZ_CONS 4
#define ENER_CONS 5

#ifdef MHD
#define BMFX_CONS 6
#define BMFY_CONS 7
#define BMFZ_CONS 8
#define BPSI_CONS 9
#define NCONSVAR_END BPSI_CONS
#else
#define NCONSVAR_END ENER_CONS
#endif

#define NCONSVAR_NUMBER (NCONSVAR_END - NCONSVAR_BEGIN + 1)
