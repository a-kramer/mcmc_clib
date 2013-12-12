/*
 *  vf_esens.cpp
 *  sde_mcmc
 *
 *  Created by Vassilios Stathopoulos on 25/01/2012.
 *  Copyright 2012 Computing Science. All rights reserved.
 *
 */

#include <fstream>
#include <sstream>
#include <string>
#include <ginac/ginac.h>
#include <ginac/matrix.h>

#include "vf.h"
#include "vf_utils.h"

using namespace std;
using namespace GiNaC;


static matrix getMatrixGrad(const matrix& M, const lst& expr, const lst& params){
	int i, j, k;
	
	int M_rows = M.rows();
	int M_cols = M.cols();
	
	int P = params.nops();

	matrix H = matrix(M_rows*M_cols, P);
	
	for (k = 0; k < P; ++k){
		symbol v = ex_to<symbol>(params[k]);
		int idx = 0;
		for (i = 0; i<M_rows; i++) {
			for (j = 0; j<M_cols; j++){
				/* substitute expressions */
				ex f;
				if (expr == 0)
					f = M(i,j);
				else
					f = iterated_subs(M(i,j),expr);
				
				ex df = f.diff(v);
				if (df != 0)
					H(idx,k) = df;
				
				idx++;
			}
		}
	}
	return H;
}


/*
 * PrintESENS -- The ESENS Code Generator
 */
void VectorField::PrintESENS(map<string,string> options){

	int N = varname_list.nops();
	int P = parname_list.nops();
	
	/* create symbolic matrix for the sensitivity variables */
	matrix S(N, P);
	int i,j;
	for (i = 0; i < N; i++){
		for (j = 0; j < P; j++){
			ostringstream varName;
			varName << "d" << varname_list[i] << "_" << parname_list[j];
			S(i,j) = symbol(varName.str());
		}
	}
	
	/* create Jacobian of the original system */
	matrix Jy = getJacobian(varvecfield_list, expreqn_list, varname_list);
	
	/* create derivatives of the original system */
	matrix Jp = getJacobian(varvecfield_list, expreqn_list, parname_list);
	
	/* calculate matrix with the rhs of the sensitivity equations */
	/* dot(S) =  Jy * S + Jp */
	matrix dV = Jy.mul(S).add(Jp);

	/* add new variables and exressions to the VF object and print */
	/* order by parameter not by variable */
	for (j = 0; j < P; j++){
		for (i = 0; i < N; i++){
			ostringstream varName;
			varName << S(i,j);
			StateVariable *var = new StateVariable( varName.str() );
			ostringstream varExpresion;
			varExpresion << dV(i,j);
			var->Formula( varExpresion.str() );
			var->DefaultInitialCondition("0.0");
			AddStateVariable(var);
		}
	}
	
	/* second order sensitivity equations */
	if (options["second"] == "yes"){
		
		/* create symbolic matrix for the 2nd order sensitivity variables */
		matrix S2(N*P, P);
		int k;
		int idx = 0;
		for (i = 0; i < N; i++){
			for (j = 0; j < P; j++){
				for (k = 0; k < P; k++){
					ostringstream varName;
					varName << "d" << varname_list[i] << "_" << parname_list[k] << "_" << parname_list[j];
					S2(idx,k) = symbol(varName.str());
				}
				idx ++;
			}
		}
		
		/* create the kronecker product matrix [Fy @ Ip] */
		matrix A(N*P,N*P);
		int idx_i, idx_j;
		
		for (i=0; i<N; i++) {
			idx_i = i*P;
			for (j=0; j<N; j++) {
				idx_j = j*P;
				for (k = 0; k<P; k++)
					A(idx_i+k,idx_j+k) = Jy(i,j);
			}	
		}
		
		/* create the kronecker product matrix [In @ S^T] */
		matrix B(N*P,N*N);
		for (k=0; k<N; k++) {
			idx_i = k*P;
			idx_j = k*N;
			for (j=0; j<P; j++)
				for (i = 0; i<N; i++)
					B(idx_i+j,idx_j+i) = S(i,j);
		}
		
		/* create matrices Fyy, Fyp, Fpy, Fpp */
		matrix Fyy = getMatrixGrad(Jy, expreqn_list, varname_list);
		matrix Fyp = getMatrixGrad(Jy, expreqn_list, parname_list);
		matrix Fpy = getMatrixGrad(Jp, expreqn_list, varname_list);
		matrix Fpp = getMatrixGrad(Jp, expreqn_list, parname_list);
		
		/* Second order sensitivities */
		/* dot(S2) = [Fy @ Ip]*S2 + [In @ S^T]*[Fyy*S + Fyp] + [Fpy*S + Fpp] */
		matrix dV2 = A.mul(S2).add(B.mul(Fyy.mul(S).add(Fyp))).add(Fpy.mul(S).add(Fpp));
		
		/* add new variables and exressions to the VF object and print */
		/* order by parameter not by variable */
		for (k = 0; k < P; k++){
			for (j = 0; j < P; j++){
				for (i = 0; i < N; i++){
					ostringstream varName;

					varName << S2(i*P +j ,k);
					
					StateVariable *var = new StateVariable( varName.str() );
					ostringstream varExpresion;
					
					varExpresion << dV2(i*P +j, k);
					
					var->Formula( varExpresion.str() );
					var->DefaultInitialCondition("0.0");
					AddStateVariable(var);
				}
			}
		}
	
	}
	
	Name(Name()+"_esens");
	PrintXML("esens");
}
