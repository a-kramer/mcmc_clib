#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
/* The error code indicates how to pre-allocate memory
 * for output values such as `f_`. The _vf function returns
 * the number of state variables, if any of the args are `NULL`.
 * GSL_SUCCESS is returned when no error occurred.
 */

/* ode vector field: y'=f(t,y;p), the Activation expression is currently unused */
int DemoModel_vf(double t, const double y_[], double f_[], void *par)
{
	double *p_=par;
	if (!y_ || !f_) return 6;
	double inv_tau=1000;
	double kf_R0=p_[0];
	double kr_R0=p_[1];
	double kf_R1=p_[2];
	double kr_R1=p_[3];
	double kf_R2=p_[4];
	double kr_R2=p_[5];
	double kf_R3=p_[6];
	double kr_R3=p_[7];
	double kf_R4=p_[8];
	double kr_R4=p_[9];
	double kf_R5=p_[10];
	double kr_R5=p_[11];
	double u=p_[12];
	double t_on=p_[13];
	double A=y_[0];
	double B=y_[1];
	double C=y_[2];
	double AB=y_[3];
	double AC=y_[4];
	double ABC=y_[5];
	double Activation=1/(1-exp(-(t-t_on)*inv_tau));
	double ReactionFlux0=u * kf_R0 * A * B - kr_R0 * AB;
	double ReactionFlux1=kf_R1 * A * C - kr_R1 * AC;
	double ReactionFlux2=kf_R2 * AB * C - kr_R2 * ABC;
	double ReactionFlux3=kf_R3 * AC * B - kr_R3 * ABC;
	f_[0] = -ReactionFlux0-ReactionFlux1;
	f_[1] = -ReactionFlux0-ReactionFlux3;
	f_[2] = -ReactionFlux1-ReactionFlux2;
	f_[3] = +ReactionFlux0-ReactionFlux2;
	f_[4] = +ReactionFlux1-ReactionFlux3;
	f_[5] = +ReactionFlux2+ReactionFlux3;
	return GSL_SUCCESS;
}
/* ode Jacobian df(t,y;p)/dy */
int DemoModel_jac(double t, const double y_[], double *jac_, double *dfdt_, void *par)
{
	double *p_=par;
	if (!y_ || !jac_) return 6*6;
	double inv_tau=1000;
	double kf_R0=p_[0];
	double kr_R0=p_[1];
	double kf_R1=p_[2];
	double kr_R1=p_[3];
	double kf_R2=p_[4];
	double kr_R2=p_[5];
	double kf_R3=p_[6];
	double kr_R3=p_[7];
	double kf_R4=p_[8];
	double kr_R4=p_[9];
	double kf_R5=p_[10];
	double kr_R5=p_[11];
	double u=p_[12];
	double t_on=p_[13];
	double A=y_[0];
	double B=y_[1];
	double C=y_[2];
	double AB=y_[3];
	double AC=y_[4];
	double ABC=y_[5];
	double Activation=1/(1-exp(-(t-t_on)*inv_tau));
/* column 1 (df/dy_0) */
	jac_[0] = -((u*kf_R0)*B)-(kf_R1*C);
	jac_[6] = -((u*kf_R0)*B)-0;
	jac_[12] = -(kf_R1*C)-0;
	jac_[18] = +((u*kf_R0)*B)-0;
	jac_[24] = +(kf_R1*C)-0;
	jac_[30] = +0+0;
/* column 2 (df/dy_1) */
	jac_[1] = -(A*(u*kf_R0))-0;
	jac_[7] = -(A*(u*kf_R0))-(kf_R3*AC);
	jac_[13] = -0-0;
	jac_[19] = +(A*(u*kf_R0))-0;
	jac_[25] = +0-(kf_R3*AC);
	jac_[31] = +0+(kf_R3*AC);
/* column 3 (df/dy_2) */
	jac_[2] = -0-(kf_R1*A);
	jac_[8] = -0-0;
	jac_[14] = -(kf_R1*A)-(kf_R2*AB);
	jac_[20] = +0-(kf_R2*AB);
	jac_[26] = +(kf_R1*A)-0;
	jac_[32] = +(kf_R2*AB)+0;
/* column 4 (df/dy_3) */
	jac_[3] = -(0-kr_R0)-0;
	jac_[9] = -(0-kr_R0)-0;
	jac_[15] = -0-(kf_R2*C);
	jac_[21] = +(0-kr_R0)-(kf_R2*C);
	jac_[27] = +0-0;
	jac_[33] = +(kf_R2*C)+0;
/* column 5 (df/dy_4) */
	jac_[4] = -0-(0-kr_R1);
	jac_[10] = -0-(kf_R3*B);
	jac_[16] = -(0-kr_R1)-0;
	jac_[22] = +0-0;
	jac_[28] = +(0-kr_R1)-(kf_R3*B);
	jac_[34] = +0+(kf_R3*B);
/* column 6 (df/dy_5) */
	jac_[5] = -0-0;
	jac_[11] = -0-(0-kr_R3);
	jac_[17] = -0-(0-kr_R2);
	jac_[23] = +0-(0-kr_R2);
	jac_[29] = +0-(0-kr_R3);
	jac_[35] = +(0-kr_R2)+(0-kr_R3);
	return GSL_SUCCESS;
}
/* ode Functions F(t,y;p) */
int DemoModel_func(double t, const double y_[], double *func_, void *par)
{
	double *p_=par;
	if (!y_ || !func_) return 3;
	double inv_tau=1000;
	double kf_R0=p_[0];
	double kr_R0=p_[1];
	double kf_R1=p_[2];
	double kr_R1=p_[3];
	double kf_R2=p_[4];
	double kr_R2=p_[5];
	double kf_R3=p_[6];
	double kr_R3=p_[7];
	double kf_R4=p_[8];
	double kr_R4=p_[9];
	double kf_R5=p_[10];
	double kr_R5=p_[11];
	double u=p_[12];
	double t_on=p_[13];
	double A=y_[0];
	double B=y_[1];
	double C=y_[2];
	double AB=y_[3];
	double AC=y_[4];
	double ABC=y_[5];
	double Activation=1/(1-exp(-(t-t_on)*inv_tau));
	func_[0] = A+AB+AC+ABC;
	func_[1] = B+AB+ABC;
	func_[2] = C+AC+ABC;
	return GSL_SUCCESS;
}
/* ode default parameters */
int DemoModel_default(double t, void *par)
{
	double *p_=par;
	if (!p_) return 14;
	double inv_tau=1000;
	p_[0] = 1.0;
	p_[1] = 1.0;
	p_[2] = 1.0;
	p_[3] = 1.0;
	p_[4] = 1.0;
	p_[5] = 1.0;
	p_[6] = 1.0;
	p_[7] = 1.0;
	p_[8] = 1.0;
	p_[9] = 1.0;
	p_[10] = 1.0;
	p_[11] = 1.0;
	p_[12] = 1;
	p_[13] = 0;
	return GSL_SUCCESS;
}
/* ode initial values */
int DemoModel_init(double t, double *y_, void *par)
{
	double *p_=par;
	if (!y_) return 6;
	double inv_tau=1000;
	double kf_R0=p_[0];
	double kr_R0=p_[1];
	double kf_R1=p_[2];
	double kr_R1=p_[3];
	double kf_R2=p_[4];
	double kr_R2=p_[5];
	double kf_R3=p_[6];
	double kr_R3=p_[7];
	double kf_R4=p_[8];
	double kr_R4=p_[9];
	double kf_R5=p_[10];
	double kr_R5=p_[11];
	double u=p_[12];
	double t_on=p_[13];
	/* the initial value of y may depend on the parameters. */
	y_[0] = 1000;
	y_[1] = 10;
	y_[2] = 10;
	y_[3] = 0;
	y_[4] = 0;
	y_[5] = 0;
	return GSL_SUCCESS;
}
