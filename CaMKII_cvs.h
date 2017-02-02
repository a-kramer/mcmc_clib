/*
 *  CaMKII_cvs.h
 *
 *  CVODES C prototype file for the functions defined in CaMKII_cvs.c
 *
 *  This file was generated by the program VFGEN (Version:2.4.1)
 *  Generated on 14-Dec-2016 at 19:00
 */

int CaMKII_vf(realtype, N_Vector, N_Vector, void *);
int CaMKII_jac(int, realtype,
                N_Vector, N_Vector,
                DlsMat, void *,
                N_Vector, N_Vector, N_Vector);
int CaMKII_sens(int, realtype, N_Vector, N_Vector,
                int, N_Vector, N_Vector,
                void *, N_Vector, N_Vector);
int CaMKII_func(realtype, N_Vector, realtype *, void *);
int CaMKII_func_sens(realtype, N_Vector, N_Vector *, double *, void *);
