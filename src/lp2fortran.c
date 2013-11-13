/*-------------------------------------------------------------------*
 * LP2FORTRAN                                                        *
 * ==========                                                        *
 * This C file offers a convenient interface to LP_SOLVE for Fortran *
 * programs.                                                         *
 *                                                                   *
 * Notice the blatant lack of error checks.  Don't blame me if this  *
 * asplodes in your face.                                            *
 *                                                                   *
 * Notice the blatant lack of documentation. This is becaus          *
 *                                                                   *
 * PLR 02.2012                                                       *
 *-------------------------------------------------------------------*/

#ifndef LP_WRAPPER
  #ifdef F90_DOUBLE_UNDERSCORE
    #ifdef F90_CAPITALS
      #define LP_WRAPPER LP_WRAPPER__
    #else
      #define LP_WRAPPER lp_wrapper__
    #endif
  #else
    #ifdef F90_NO_UNDERSCORE
      #ifdef F90_CAPITALS
        #define LP_WRAPPER LP_WRAPPER
      #else
        #define LP_WRAPPER lp_wrapper
      #endif
    #else
      #ifdef F90_CAPITALS
        #define LP_WRAPPER LP_WRAPPER_
      #else
        #define LP_WRAPPER lp_wrapper_
      #endif
    #endif
  #endif
#endif

#include "lp_lib.h"

void LP_WRAPPER(int *neqn, int *nvar,
                int *vartype,
                REAL *lparray, int *ineq_sign, REAL *lprhs,
                REAL *targetcoeff, int *target_opt_sign,
                REAL *solution, int *solution_status) {
 /* Define variables */
 lprec* lp;
 int *colno=NULL;
 REAL *row=NULL;
 int ieqn; int ivar; int icol;

 /* Initialize linear problem */
 lp=make_lp(0,*nvar);
 for (ivar=0;ivar<=*nvar-1;ivar++) {
  if (vartype[ivar]==1) {
   set_int(lp,ivar+1,TRUE);
  } else if (vartype[ivar]==2) {
   set_binary(lp,ivar+1,TRUE);
  }
 }
 set_add_rowmode(lp,TRUE);
 colno=(int*)malloc(*nvar*sizeof(*colno));
 row=(REAL*)malloc(*nvar*sizeof(*row));

 /* Create equations (inequations, really) */
 for (ieqn=0;ieqn<=*neqn-1;ieqn++) {
  icol=0;
  for (ivar=0;ivar<=*nvar-1;ivar++) {
   if (lparray[ieqn*(*nvar)+ivar]!=0.0) {
    colno[icol]=ivar+1; /* colno is 1-based */
    row[icol++]=lparray[ieqn*(*nvar)+ivar];
   }
  }
  if (ineq_sign[ieqn]>0) {
   add_constraintex(lp,icol,row,colno,GE,lprhs[ieqn]);
  } else if (ineq_sign[ieqn]==0) {
   add_constraintex(lp,icol,row,colno,EQ,lprhs[ieqn]);
  } else {
   add_constraintex(lp,icol,row,colno,LE,lprhs[ieqn]);
  }
 }
 set_add_rowmode(lp,FALSE);

 /* Define objective function */
 icol=0;
 for (ivar=0;ivar<=*nvar-1;ivar++) {
  if (targetcoeff[ivar]!=0.0) {
   colno[icol]=ivar+1; /* colno is 1-based */
   row[icol++]=targetcoeff[ivar];
  }
 }
 set_obj_fnex(lp,icol,row,colno);

 /* Solve */
 if (*target_opt_sign>=0) {
  set_maxim(lp);
 } else {
  set_minim(lp);
 }
 set_verbose(lp, IMPORTANT);
 *solution_status=solve(lp);
 get_variables(lp,solution);

 /* Clean up */
 if(row!=NULL)free(row);
 if(colno!=NULL)free(colno);
 if(lp!=NULL)delete_lp(lp);

 return;
}
