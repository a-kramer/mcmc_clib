/*
 *  ODEmodel11S26P4U_cvs.h
 *
 *  CVODES C prototype file for the functions defined in ODEmodel11S26P4U_cvs.c
 *
 *  This file was generated by the program VFGEN (Version:2.4.1)
 *  Generated on 27-Jan-2017 at 17:06
 */

int ODEmodel11S26P4U_vf(realtype, N_Vector, N_Vector, void *);
int ODEmodel11S26P4U_jac(int, realtype,
                N_Vector, N_Vector,
                DlsMat, void *,
                N_Vector, N_Vector, N_Vector);
int ODEmodel11S26P4U_sens(int, realtype, N_Vector, N_Vector,
                int, N_Vector, N_Vector,
                void *, N_Vector, N_Vector);
int ODEmodel11S26P4U_func(realtype, N_Vector, realtype *, void *);
int ODEmodel11S26P4U_func_sens(realtype, N_Vector, N_Vector *, double *, void *);
