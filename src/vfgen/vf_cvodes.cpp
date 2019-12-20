
//
//  vf_cvodes.cpp
//
//  This file defines the VectorField::PrintCVODES method.
//
//
//
//  Copyright (C) 2011 Vassilios Stathopoulos
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License, Version 2, as
//  published by the Free Software Foundation.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with this program.  The file LICENSE-GPL2 in the VFGEN source directory
//  contains a copy. You may also obtain a copy by writing to the Free Software
//  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
//

#include <fstream>
#include <string>
#include <ginac/ginac.h> 
#include <ginac/matrix.h>

#include "vf.h"
#include "vf_utils.h"

using namespace std;
using namespace GiNaC;


/*
 * PrintCVODES -- The CVODES Code Generator
 */
void VectorField::PrintCVODES(map<string,string> options){

    int nc, np, nv, na, nf;

    if (options.count("version") > 0 && options["version"] != "2.6.0"){
        cerr << "vfgen CVODES command: unknown version specified: " << options["version"] << endl;
        cerr << "Versions of CVODES supported by VFGEN are 2.6.0. Default: version=2.6.0" << endl;
        exit(-1);
	} 

    symbol t(IndependentVariable);
    nc = conname_list.nops();
    nv = varname_list.nops();
    np = parname_list.nops();
    na = exprname_list.nops();
    nf = funcname_list.nops();

    string filename = Name()+"_cvs.c";
    ofstream fout;
    fout.open(filename.c_str());
    fout << csrc << left;

    string pfilename = Name()+"_cvs.h";
    ofstream pout;
    pout.open(pfilename.c_str());
    pout << csrc << left;

    /*  Print C file and header information. */
    fout << "/*" << endl;
    fout << " *  " << filename << endl;
    fout << " *" << endl;
    fout << " *  CVODES C file for the vector field named: " << Name() << endl;
    fout << " *" << endl;
    PrintVFGENComment(fout," *  ");
    fout << " */" << endl;
    fout << endl;

    pout << "/*" << endl;
    pout << " *  " << pfilename << endl;
    pout << " *" << endl;
    pout << " *  CVODES C prototype file for the functions defined in " << filename << endl;
    pout << " *" << endl;
    PrintVFGENComment(pout," *  ");
    pout << " */" << endl;
    pout << endl;

    fout << "#include <math.h>" << endl;
	
	/* Default is CVODES v2.6.0 */
    /* Headers for CVODES */
	fout << "/* Include headers for SUNDIALS v2.4.0, CVODES v2.6.0 */" << endl;
	fout << "#include <cvodes/cvodes.h> " << endl;	
	fout << "#include <cvodes/cvodes_dense.h> " << endl;	
	fout << "#include <nvector/nvector_serial.h> " << endl;	
	fout << "#include <sundials/sundials_types.h> " << endl;
	/* Declerations for ode_model.h	*/
	fout << endl;
	fout << "typedef int (*rhs_f)(realtype, N_Vector, N_Vector, void *);" << endl;
	fout << "typedef int (*jac_f)(long int, realtype," << endl;
	fout << "                     N_Vector, N_Vector," << endl;
	fout << "                     DlsMat, void*," << endl;
	fout << "                     N_Vector, N_Vector, N_Vector);" << endl;
	fout << endl;
	fout << "typedef int (*jacp_f)(int, realtype," << endl;
	fout << "                     N_Vector, N_Vector," << endl;
	fout << "                     DlsMat, void*," << endl;
	fout << "                     N_Vector, N_Vector, N_Vector);" << endl;
	fout << endl;
	fout << "typedef int (*rhs_sens)(int, realtype, N_Vector, N_Vector," << endl;
	fout << "                        int, N_Vector, N_Vector," << endl;
	fout << "                        void *, N_Vector, N_Vector);" << endl;
	fout << endl;
	fout << "typedef int (*func)(realtype, N_Vector, realtype *, void *);" << endl;
	fout << "typedef int (*func_sens)(realtype, N_Vector, N_Vector *, double *, void *);" << endl;
	fout << endl;
	fout << "typedef struct{" << endl;
	fout << "const int N;" << endl;
	fout << "const int P;" << endl;
	fout << "const int F;" << endl;
	fout << "const double* v;" << endl;
	fout << "const double* p;" << endl;
	fout << "rhs_f vf_eval;" << endl;
	fout << "jac_f vf_jac;" << endl;
	fout << "jacp_f vf_jacp;" << endl;
	fout << "rhs_sens vf_sens;" << endl;
	fout << "func vf_func;" << endl;
	fout << "func_sens vf_func_sens;" << endl;
	fout << "const char** v_names;" << endl;
	fout << "const char** p_names;" << endl;
	fout << "const char** f_names;" << endl;
	fout << "void*	dylib;"<< endl;
	fout << "} ode_model;" << endl;
	fout << endl;
	
	if (HasPi)
		fout << "#define Pi  M_PI\n";
	
	/* ================================
	 * Print the vector field function.
	 * ================================  
     */
    fout << "/*" << endl;
    fout << " *  The vector field." << endl;
    fout << " */" << endl;
    fout << endl;
    string func_return_type = "int";
    
    fout << func_return_type << " " << Name() << "_vf(realtype t, N_Vector y_, N_Vector f_, void *params)" << endl;
    pout << func_return_type << " " << Name() << "_vf(realtype, N_Vector, N_Vector, void *);" << endl;
    fout << "    {" << endl;

	/* constants */
	for (int i = 0; i < nc; ++i){
		fout << "    const realtype " << conname_list[i] << " = RCONST(" << convalue_list[i] << ");" << endl;
	}
	
	/* declerations for variables, parameters and expressions */
    CDeclare(fout,"realtype",varname_list);
    CDeclare(fout,"realtype",parname_list);
    CDeclare(fout,"realtype",exprname_list);
    fout << "    realtype *p_;" << endl;
    fout << endl;
    fout << "    p_ = (realtype *) params;" << endl;
    fout << endl;
    
	/* definitions for variables */
    for (int i = 0; i < nv; ++i){
        fout << "    ";
        fout.width(10);
        fout << varname_list[i];
        fout.width(0);
        fout << " = NV_Ith_S(y_," << i << ");" << endl;
	}
	fout << endl;
	
	/* definitions for parameters*/
    GetFromVector(fout,"    ",parname_list,"p_","[]",0,";");
    fout << endl;
	
	/* definitions for expressions */
    for (int i = 0; i < na; ++i){
		fout << "    " << exprname_list[i] << " = " << exprformula_list[i] << ";" << endl;
	}
    if (na > 0)
        fout << endl;

    /* RHS of the vector field */
	for (int i = 0; i < nv; ++i){
        fout << "    NV_Ith_S(f_," << i << ") = " << varvecfield_list[i] << ";" << endl;
	}
	
    /* end of vector field function */
	fout << "    return 0;\n";
    fout << "    }" << endl;
    fout << endl;
    
		
	/* ============================
	 * Print the Jacobian function.
	 * ============================  
     */
    fout << "/*" << endl;
    fout << " *  The Jacobian." << endl;
    fout << " */" << endl;
    fout << endl;
    fout << func_return_type << " " << Name() << "_jac(long N_, realtype t," << endl;
    fout << "                N_Vector y_, N_Vector fy_," << endl;
	fout << "                DlsMat jac_, void *params," << endl;
    fout << "                N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)" << endl;

    pout << func_return_type << " " << Name() << "_jac(long, realtype," << endl;
    pout << "                N_Vector, N_Vector," << endl;
	pout << "                DlsMat, void *," << endl;
    pout << "                N_Vector, N_Vector, N_Vector);" << endl;
    fout << "    {" << endl;

	/* constants */
	for (int i = 0; i < nc; ++i){
		fout << "    const realtype " << conname_list[i] << " = RCONST(" << convalue_list[i] << ");" << endl;
	}
	
	/* declarations for variables and parameters */
    CDeclare(fout,"realtype",varname_list);
    CDeclare(fout,"realtype",parname_list);
    fout << "    realtype *p_;" << endl;
    fout << endl;
    fout << "    p_ = (realtype *) params;" << endl;
    fout << endl;
	
    /* definitions for variables */
    for (int i = 0; i < nv; ++i){
        fout << "    ";
        fout.width(10);
        fout << varname_list[i];
        fout.width(0);
        fout << " = NV_Ith_S(y_," << i << ");" << endl;
	}
	fout << endl;

	/* definitions for paramters */
    GetFromVector(fout,"    ",parname_list,"p_","[]",0,";");
    fout << endl;
		
	/* calculate Jacobian symbolically */
	matrix Jac_y = getJacobian(varvecfield_list, expreqn_list, varname_list);
	/* RHS of the Jacobian function */
    for (int i = 0; i < nv; ++i){
        for (int j = 0; j < nv; ++j){            
            // Skip zero elements.  CVODE initializes jac_ to zero before calling the Jacobian function.
			if (Jac_y(i,j) != 0)
                fout << "    DENSE_ELEM(jac_, " << i << ", " << j << ") = " << Jac_y(i,j) << ";" << endl;
		}
	}
	/* end of jacobian function */
	fout << "    return 0;\n";
    fout << "    }" << endl;

    // print Parameter Derivative, dvf/dp:
    fout << "/*" << endl;
    fout << " *  The Parameter Derivative." << endl;
    fout << " */" << endl;
    fout << endl;
    fout << func_return_type << " " << Name() << "_jacp(int N_, realtype t," << endl;
    fout << "                N_Vector y_, N_Vector fy_," << endl;
	fout << "                DlsMat jacp_, void *params," << endl;
    fout << "                N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)" << endl;

    pout << func_return_type << " " << Name() << "_jacp(int, realtype," << endl;
    pout << "                N_Vector, N_Vector," << endl;
	pout << "                DlsMat, void *," << endl;
    pout << "                N_Vector, N_Vector, N_Vector);" << endl;
    fout << "    {" << endl;

    for (int i = 0; i < nc; ++i)
        {
        fout << "    const realtype " << conname_list[i] << " = RCONST(" << convalue_list[i] << ");" << endl;
        }
    CDeclare(fout,"realtype",varname_list);
    CDeclare(fout,"realtype",parname_list);
    fout << "    realtype *p_;" << endl;
    fout << endl;
    fout << "    p_ = (realtype *) params;" << endl;
    fout << endl;
    // GetFromVector(fout,"    ",varname_list,"y_","[]",0,";");
    for (int i = 0; i < nv; ++i)
        {
        fout << "    ";
        fout.width(10);
        fout << varname_list[i];
        fout.width(0);
        fout << " = NV_Ith_S(y_," << i << ");" << endl;
        }
    fout << endl;
    GetFromVector(fout,"    ",parname_list,"p_","[]",0,";");
    fout << endl;
    for (int i = 0; i < nv; ++i)
        {
        ex f = iterated_subs(varvecfield_list[i],expreqn_list);
        for (int j = 0; j < np; ++j)
            {
	    symbol p = ex_to<symbol>(parname_list[j]);
            ex df = f.diff(p);
            // Skip zero elements.  CVODE initializes jac_ to zero before calling the Jacobian function.
            if (df != 0)
                fout << "    DENSE_ELEM(jacp_, " << i << ", " << j << ") = " << df << ";" << endl;
            }
        }
    if (options["version"] != "2.3.0")
        fout << "    return 0;\n";
    fout << "    }" << endl;
    

	if (options["sens"] == "yes"){
		
		/* ==================================
		 * Print the Sensitivities functions.
		 * ==================================  
		 */
		fout << endl;
        fout << "/*" << endl;
        fout << " *  Sensitivity functions " << endl;
        fout << " */" << endl;
        fout << endl;
		
		fout << func_return_type << " " << Name() << "_sens(int Ns_, realtype t, N_Vector y_, N_Vector ydot_," << endl;
		fout << "                int iS_, N_Vector yS_, N_Vector ySdot_," << endl;
		fout << "                void *params, N_Vector tmp1, N_Vector tmp2)" << endl;
		
		pout << func_return_type << " " << Name() << "_sens(int, realtype, N_Vector, N_Vector," << endl;
		pout << "                int, N_Vector, N_Vector," << endl;
		pout << "                void *, N_Vector, N_Vector);" << endl;
		fout << "    {" << endl;
		
		/* constants */
		for (int i = 0; i < nc; ++i){
			fout << "    const realtype " << conname_list[i] << " = RCONST(" << convalue_list[i] << ");" << endl;
        }
		
		/* declarations for variables and parameters */
		CDeclare(fout,"realtype", varname_list);
		CDeclare(fout,"realtype", parname_list);

		/* create symbolic sensitivity variables */
		matrix S;
		S = ex_to<matrix>(symbolic_matrix(nv, 1, "S"));
		
		/* declarations for sensitivity variables */
		fout << "    realtype ";
		for (int i = 0; i < nv; ++i) {
			if (i > 0)
				fout << ", ";
			fout << S(i,0);	
		}
		fout << ";" << endl;
		
		/* type cast parameters to double array*/
		fout << "    realtype *p_;" << endl;
		fout << endl;
		fout << "    p_ = (realtype *) params;" << endl;
		fout << endl;
		
		/* definitions for variables */
		for (int i = 0; i < nv; ++i){
			fout << "    ";
			fout.width(10);
			fout << varname_list[i];
			fout.width(0);
			fout << " = NV_Ith_S(y_," << i << ");" << endl;
        }
		
		/* definitions for sensitivities */
		for (int i = 0; i < nv; ++i){ 
			fout << "    ";
			fout.width(10);
			fout << S(i,0);
			fout.width(0);
			fout << " = NV_Ith_S(yS_," << i << ");" << endl;
        }
		fout << endl;
		
		/* definitions for parameters */
		GetFromVector(fout,"    ",parname_list,"p_","[]",0,";");
		fout << endl;
		
		/* print RHS Jac_y*yS */
		matrix JyS = Jac_y.mul(S);
		for (int i = 0; i < nv; ++i){
			if (JyS(i,0) != 0)
				fout << "    NV_Ith_S(ySdot_," << i << ") = " << JyS(i,0) << ";" << endl;
        }
		fout << endl;
		
		// print RHS + df/dp
		matrix J_p = getJacobian(varvecfield_list, expreqn_list, parname_list);
		
		fout << "	switch (iS_) {" << endl;
		for (int i = 0; i < np; ++i) {
			fout << "        ";
			fout << "case " << i << ":" << endl;
			for (int j = 0; j < nv; ++j) {
				if (J_p(j,i) != 0){
					fout << "	        ";
					fout << "    NV_Ith_S(ySdot_," << j << ") += " << J_p(j,i) << ";" << endl;
				}
			}
			fout << "	        break;" << endl;
		}
		fout << "	 }" << endl;
		
		/* end of sensitivities function */
		fout << "    return 0;\n";
		fout << "    }" << endl;
		
	}
	
    if (options["func"] == "yes" & nf > 0){
        /* ================================================
		 * Print the user-defined functions.
		 * A single function is created that puts all the
		 * user-defined function values in an array.  This
		 * function is defined so that it can be used with
		 * the CVODE rootfinding code.
		 * ================================================
		 */
        fout << endl;
        fout << "/*" << endl;
        fout << " *  User-defined functions. " << endl;
        fout << " */" << endl;
        fout << endl;
        fout << func_return_type << " " << Name() << "_func(realtype t, N_Vector y_, realtype *func_, void *params)" << endl;
        pout << func_return_type << " " << Name() << "_func(realtype, N_Vector, realtype *, void *);" << endl;
        fout << "    {" << endl;
        
		/* constants */
        for (int i = 0; i < nc; ++i){
			fout << "    const realtype " << conname_list[i] << " = RCONST(" << convalue_list[i] << ");" << endl;
		}
			
		/* declarations for variables, parameters and expressions */
        CDeclare(fout,"realtype",varname_list);
        CDeclare(fout,"realtype",parname_list);
        CDeclare(fout,"realtype",exprname_list);
        fout << "    realtype *p_;" << endl;
        fout << endl;
        fout << "    p_ = (realtype *) params;" << endl;
        fout << endl;
			
		/* declarations for variables */
        for (int i = 0; i < nv; ++i){
			fout << "    ";
            fout.width(10);
            fout << varname_list[i];
            fout.width(0);
            fout << " = NV_Ith_S(y_," << i << ");" << endl;
		}
        fout << endl;

		/* declarations for parameters */
        GetFromVector(fout,"    ",parname_list,"p_","[]",0,";");
        fout << endl;
		
		/* declarations for expressions */
        for (int i = 0; i < na; ++i){
            fout << "    " << exprname_list[i] << " = " << exprformula_list[i] << ";" << endl;
		}
        if (na > 0)
            fout << endl;
		
		/* print RHS for functions */
        for (int n = 0; n < nf; ++n){
            fout << "    /* " << funcname_list[n] << ":  */" << endl;
            fout << "    func_[" << n << "] = " << funcformula_list[n] << ";" << endl;
		}
		
		/* end of the user-defined functions*/
        fout << "    return 0;\n";
        fout << "    }" << endl;
	}
		
	if (options["func"] == "yes" & nf > 0 & options["sens"] == "yes"){
		/* =================================================
		 * Print derivatives for the user-defined functions.
		 * A single function is created that puts all the
		 * sensitivity values in an array.
		 * =================================================
		 */
		fout << endl;
		fout << "/*" << endl;
		fout << " *  Sensitivities of User-defined functions. " << endl;
		fout << " */" << endl;
		fout << endl;
		fout << func_return_type << " " << Name() << "_func_sens(realtype t, N_Vector y_, N_Vector * yS_, double* ret_, void *params)" << endl;
		pout << func_return_type << " " << Name() << "_func_sens(realtype, N_Vector, N_Vector *, double *, void *);" << endl;
		fout << "    {" << endl;
		
		/* constants */
		for (int i = 0; i < nc; ++i){
			fout << "    const realtype " << conname_list[i] << " = RCONST(" << convalue_list[i] << ");" << endl;
		}
		
		/* declarations for variables, parameters and expressions */
		CDeclare(fout,"realtype",varname_list);
		CDeclare(fout,"realtype",parname_list);
		CDeclare(fout,"realtype",exprname_list);
		fout << "    realtype *p_;" << endl;
		fout << endl;
		fout << "    p_ = (realtype *) params;" << endl;
		fout << endl;
		
		/* definitions for variables */
		for (int i = 0; i < nv; ++i){
			fout << "    ";
			fout.width(10);
			fout << varname_list[i];
			fout.width(0);
			fout << " = NV_Ith_S(y_," << i << ");" << endl;
		}
		fout << endl;
		
		/* definitions for parameters */
		GetFromVector(fout,"    ",parname_list,"p_","[]",0,";");
		fout << endl;
		
		/* definitions for expressions */
		for (int i = 0; i < na; ++i){
			fout << "    " << exprname_list[i] << " = " << exprformula_list[i] << ";" << endl;
		}
		if (na > 0)
			fout << endl;
		
		/* declarations of sensitivity variables */
		matrix yS;
		yS = ex_to<matrix>(symbolic_matrix(np,nv,"yS_"));										
		fout << "    realtype ";
		int k = 0;
		for (int j = 0; j < np; ++j) {
			for (int i = 0; i < nv; ++i){
				if (k > 0)
					fout << ", ";
				fout << yS(j,i);
				k++;
			}
		}
		fout << ";" << endl;
		
		// definitions of the sensitivity variables */
		for (int j = 0; j < np; ++j){
			for (int i = 0; i < nv; ++i){
				fout <<"    "<< yS(j,i) << " = NV_Ith_S(yS_[" << j << " ]," << i << ");\n";
			}
		}
		fout << endl;
		
		/* compute function derivatives symbolically */
		/* get dG/dy * dy/dp^T + dG/dp */
		matrix dG_y = getJacobian(funcformula_list, expreqn_list, varname_list);
		matrix dG_p = getJacobian(funcformula_list, expreqn_list, parname_list);
		matrix func_result = dG_y.mul(yS.transpose()).add(dG_p);
		
		/* print RHS for function derivatives */
		k = 0;
		for (int i = 0; i < np; ++i){
			fout << "    /* " << parname_list[i] << ":  */" << endl;
			
			for (int n = 0; n < nf; ++n){
				fout << "        ret_[" << k << "] = " << func_result(n,i) << ";" << endl;
				k++;
			}
		}
		
		/* end of derivatives for the user-defined functions. */
		fout << "    return 0;\n";
		fout << "    }" << endl;
	}
	
	
	/* =========================================
	 * Static variable definitions for initial
	 * conditions, default parameters, variable
	 * names etc.
	 * =========================================
	 */
	fout <<  endl;
	fout <<"/* Initial conditions, names and ode_model structure */" <<endl;
	
	/* default initial conditions */
	fout << "const static double " << Name() << "_init_v[" << nv << "]={" <<endl;
	for (int i = 0; i < nv; ++i){
		fout << "                 RCONST(" << vardefic_list[i] << ")";
		if (i != nv-1)
			fout << ",\n" ;
	}
	fout << "\n            };\n" ;
	fout <<  endl;
	
	/* default parameter values */
	fout << "const static double " << Name() << "_init_p[" << np << "]={" <<endl;
	for (int i = 0; i < np; ++i){             
		fout << "                 RCONST(" << pardefval_list[i] << ")";
		if (i != np-1)
			fout << ",\n" ;
	}
	fout << "\n            };\n" ;
	
	/* variable names */
	string tmp_names = Name()+"_varnames";
	MakeConstCArrayOfStrings(fout, tmp_names, varname_list);

	/* parameter names */
	if (np > 0){
		tmp_names = Name()+"_parnames";
		MakeConstCArrayOfStrings(fout, tmp_names, parname_list);
	}
	else
		fout << "const static char *"<< Name() << "_parnames[] = {\"\"};\n";
	
	/* function names */
	if (nf > 0){
		tmp_names = Name()+"_funcnames";
		MakeConstCArrayOfStrings(fout, tmp_names, funcname_list);
	}
	else
		fout << "const static char *" << Name() << "_funcnames[] = {\"\"};\n";
	fout <<  endl;
	
	/* ode_model data structure */
	fout << "ode_model  " << Name() << "_odeModel = {"<< nv <<", "<< np <<", "<< nf <<",\n";
	fout << "                     "<< Name() <<"_init_v, "<< Name() <<"_init_p,\n";
	fout << "                     &"<< Name() <<"_vf, &"<< Name() <<"_jac, &" << Name() <<"_jacp,\n";
	if(options["sens"] == "yes")
		fout << "                     &"<< Name() <<"_sens,\n";
	else 
		fout << "                     0,\n";
	if(options["func"] == "yes" & nf > 0 )
		fout << "                     &"<< Name() <<"_func,\n";
	else
		fout << "                     0,\n";
	if(options["func"] == "yes" & nf >0 & options["sens"] == "yes")
		fout << "                     &"<< Name() <<"_func_sens,\n";
	else
		fout << "                     0,\n";
	fout << "                     "<< Name()<<"_varnames, "<< Name()<<"_parnames,\n";
	fout << "                     "<< Name()<<"_funcnames, 0\n";
	fout << "                     };\n";
		
	
    fout.close();
    pout.close();

	
    if (options["demo"] == "yes"){
		
        /* =========================================================
         * Create a self-contained ODE solver for this vector field
         * that allows the user to give the initial conditions,
         * parameter values, and some solver control parameters
         * on the command line.
		 * =========================================================
         */

        string tfilename = Name()+"_cvsdemo.c";
        ofstream tout;
        tout.open(tfilename.c_str());
        tout << csrc << left;
        tout << "/*\n";
        tout << " *  " << tfilename << endl;
        tout << " *\n" ;
        tout << " *" << endl;
        tout << " *  CVODES ODE solver for the vector field named: " << Name() << endl;
        tout << " *" << endl;
        PrintVFGENComment(tout," *  ");
        tout << " *\n";
        tout << " */" << endl;
        tout << endl;
        tout << "#include <stdlib.h>" << endl;
        tout << "#include <string.h>" << endl;
        tout << "#include <math.h>" << endl;
        // Default is CVODE v2.6.0
		tout << "/* Include headers for SUNDIALS v2.4.0, CVODE v2.6.0 */" << endl;
		tout << "#include <sundials/sundials_types.h>" << endl;
		tout << "#include <nvector/nvector_serial.h>" << endl;
		tout << "#include <cvodes/cvodes.h>" << endl;
		tout << "#include <cvodes/cvodes_dense.h>" << endl;

        tout << endl;
        tout << "#include \"" << pfilename << "\"\n";
        tout << endl << endl;
		tout << "    #define ZERO  RCONST(0.0)" << endl;
		tout << endl << endl;
			
        tout << "int use(int argc, char *argv[], int nv, char *vname[], double y_[], int np, char *pname[], const double p_[])\n" ;
        tout << "    {\n" ;
        tout << "    int i;\n" ;
        tout << "    printf(\"use: %s [options]\\n\",argv[0]);\n" ;
        tout << "    printf(\"options:\\n\");\n" ;
        tout << "    printf(\"    -h    Print this help message.\\n\");\n" ;
        tout << "    for (i = 0; i < nv; ++i)\n" ;
        tout << "        printf(\"    %s=<initial_condition>   Default value is %e\\n\",vname[i],y_[i]);\n";
        tout << "    for (i = 0; i < np; ++i)\n" ;
        tout << "        printf(\"    %s=<parameter_value>   Default value is %e\\n\",pname[i],p_[i]);\n";
        tout << "    printf(\"    abserr=<absolute_error_tolerance>\\n\");\n" ;
        tout << "    printf(\"    relerr=<relative_error_tolerance>\\n\");\n" ;
        tout << "    printf(\"    stoptime=<stop_time>\\n\");\n" ;
        tout << "    return 0;\n";
        tout << "    }\n" ;
        tout << endl << endl;
        tout << "int assign(char *str[], int ns, double v[], char *a)\n" ;
        tout << "    {\n" ;
        tout << "    int i;\n" ;
        tout << "    char name[256];\n" ;
        tout << "    char *e;\n" ;
        tout << "\n" ;
        tout << "    e = strchr(a,'=');\n" ;
        tout << "    if (e == NULL)\n" ;
        tout << "        return(-1);\n" ;
        tout << "    *e = '\\0';\n" ;
        tout << "    strcpy(name,a);\n" ;
        tout << "    *e = '=';\n" ;
        tout << "    ++e;\n" ;
        tout << "    for (i = 0; i < ns; ++i)\n" ;
        tout << "        if (strcmp(str[i],name)==0)\n" ;
        tout << "            break;\n" ;
        tout << "    if (i == ns)\n" ;
        tout << "        return(-1);\n" ;
        tout << "    v[i] = atof(e);\n" ;
        tout << "    return(i);\n" ;
        tout << "    }\n" ;
        tout << endl << endl;
        tout << "int main (int argc, char *argv[])\n" ;
        tout << "    {\n" ;
        tout << "    int i,j;\n";
        tout << "    int flag;\n";
        tout << "    const int N_ = " << nv << ";\n" ;
        tout << "    const int P_ = " << np << ";\n" ;
        for (int i = 0; i < nc; ++i)
            {
            tout << "    const realtype " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
            }
        tout << "    const realtype def_p_[" << np << "] = {" ;
        for (int i = 0; i < np; ++i)
            {
            tout << "RCONST(" << pardefval_list[i] << ")" ;
            if (i != np-1)
                tout << ", " ;
            }
        tout << "};\n" ;
        // CDeclare(tout,"realtype",parname_list);
        GetFromVector(tout,"    const realtype ",parname_list,"def_p_","[]",0,";");
        tout << "    realtype def_y_[" << nv << "] = {";
        for (int i = 0; i < nv; ++i)
            {
            // tout << def_var_value.at(i) ;
            // tout << "RCONST(0.0)" ;
            tout << "RCONST(" << vardefic_list[i] << ")";
            if (i != nv-1)
                tout << ", " ;
            }
        tout << "};\n" ;
        tout << "    realtype y_[" << nv << "];\n" ;

        tout << "    realtype p_[" << np << "];\n" ;
        tout << "    realtype solver_param_[3] = {RCONST(1.0e-6),RCONST(0.0),RCONST(10.0)};\n" ;
        MakeCArrayOfStrings(tout,"varnames_",varname_list);
        if (np > 0)
            MakeCArrayOfStrings(tout,"parnames_",parname_list);
        else
            tout << "    char *parnames_[] = {\"\"};\n";
        tout << "    char *solver_param_names_[3] = {\"abserr\",\"relerr\",\"stoptime\"};\n" ;
        tout << endl;
        tout << "    for (i = 0; i < N_; ++i)\n" ;
        tout << "        y_[i] = def_y_[i];\n" ;
        tout << "    for (i = 0; i < P_; ++i)\n" ;
        tout << "        p_[i] = def_p_[i];\n" ;
        tout << "    for (i = 1; i < argc; ++i)\n" ;
        tout << "        {\n" ;
        tout << "        int j;\n" ;
        tout << "        if (strcmp(argv[i],\"-h\") == 0)\n" ;
        tout << "            {\n" ;
        tout << "            use(argc,argv,N_,varnames_,def_y_,P_,parnames_,def_p_);\n" ;
        tout << "            exit(0);\n" ;
        tout << "            }\n" ;
        tout << "        j = assign(varnames_,N_,y_,argv[i]);\n" ;
        tout << "        if (j == -1)\n" ;
        tout << "            {\n" ;
        tout << "            j = assign(parnames_,P_,p_,argv[i]);\n" ;
        tout << "            if (j == -1)\n" ;
        tout << "                {\n" ;
        tout << "                j = assign(solver_param_names_,3,solver_param_,argv[i]);\n" ;
        tout << "                if (j == -1)\n" ;
        tout << "                    {\n" ;
        tout << "                    fprintf(stderr,\"unknown argument: %s\\n\",argv[i]);\n" ;
        tout << "                    use(argc,argv,N_,varnames_,def_y_,P_,parnames_,def_p_); \n";
        tout << "                    exit(-1);\n" ;
        tout << "                    }\n" ;
        tout << "                }\n" ;
        tout << "            }\n" ;
        tout << "        }\n" ;
        tout << endl;
        tout << "    N_Vector y0_;\n";
        tout << "    y0_ = N_VNew_Serial(N_);\n";
        tout << "    for (i = 0; i < N_; ++i)\n";
        tout << "        NV_Ith_S(y0_,i) = y_[i];\n";
        tout << endl;

        tout << "    /* For non-stiff problems:   */\n";
        tout << "    /* void *cvode_mem = CVodeCreate(CV_ADAMS,CV_FUNCTIONAL);*/\n";
        tout << "    /* For stiff problems:       */\n";
        tout << "    void *cvode_mem = CVodeCreate(CV_BDF,CV_NEWTON);\n";
        tout << endl;
		tout << "    realtype t = RCONST(0.0);\n";
		tout << "    flag = CVodeInit(cvode_mem,"<< Name() << "_vf, t, y0_ );\n";	
		tout << "    flag = CVodeSStolerances(cvode_mem, solver_param_[1], solver_param_[0]);\n";	
		tout << "    flag = CVodeSetUserData(cvode_mem, &(p_[0]));\n";	
		tout << "    flag = CVDense(cvode_mem, N_);\n";
		tout << "    flag = CVDlsSetDenseJacFn(cvode_mem," << Name() << "_jac);\n";
		
		if (options["sens"] == "yes"){
			tout << "     realtype pbar[P_];\n";
			tout << "     int is;\n";
			tout << "     N_Vector *yS;\n";
			tout << "     for (i = 0; i < P_; ++i)\n";
			tout << "         pbar[i] = p_[i];\n";
			tout << "     yS = N_VCloneVectorArray_Serial(P_, y0_);\n";
			tout << "     for (is = 0; is < P_; is++) N_VConst(ZERO, yS[is]);\n";
			tout << "     flag = CVodeSensInit1(cvode_mem, P_, CV_STAGGERED1," << Name() << "_sens, yS);\n";
			tout << "     flag = CVodeSensEEtolerances(cvode_mem);\n";
			tout << "     flag = CVodeSetSensErrCon(cvode_mem, TRUE);\n";
			tout << "     flag = CVodeSetSensParams(cvode_mem, NULL, pbar, NULL);\n";
		}
		
    
        tout << endl;
        tout << "    realtype t1 = solver_param_[2];\n" ;
        tout << "    /* Print the initial condition  */\n";
        tout << "    printf(\"%.8e\",t);\n" ;
        tout << "    for (j = 0; j < N_; ++j)\n" ;
        tout << "        printf(\" %.8e\",NV_Ith_S(y0_,j));\n" ;
        if (options["func"] == "yes" & nf > 0)
            {
            tout << "    realtype funcval[" << nf << "];\n";
            tout << "    " << Name() << "_func(t,y0_,funcval,(void *) p_);\n";
            for (int i = 0; i < nf; ++i)
                 tout << "    printf(\" %.8e\",funcval[" << i << "]);\n";
            }
        tout << "    printf(\"\\n\");\n";
        tout << "    flag = CVodeSetStopTime(cvode_mem,t1);\n";
        tout << "    while (t < t1)\n";
        tout << "        {\n" ;
        tout << "        /* Advance the solution */\n";
        tout << "        flag = CVode(cvode_mem,t1,y0_,&t,CV_ONE_STEP);\n";
        tout << "        if (flag != CV_SUCCESS && flag != CV_TSTOP_RETURN)\n";
        tout << "            {\n" ;
        tout << "            fprintf(stderr,\"flag=%d\\n\",flag);\n" ;
        tout << "            break;\n" ;
        tout << "            }\n" ;
        tout << "        /* Print the solution at the current time */\n";
        tout << "        printf(\"%.8e\",t);\n" ;
        tout << "        for (j = 0; j < N_; ++j)\n" ;
        tout << "            printf(\" %.8e\",NV_Ith_S(y0_,j));\n" ;

        if (options["func"] == "yes" & nf > 0)
            {
            tout << "        " << Name() << "_func(t,y0_,funcval,(void *) p_);\n";
            for (int i = 0; i < nf; ++i)
                 tout << "        printf(\" %.8e\",funcval[" << i << "]);\n";
            }
        tout << "        printf(\"\\n\");\n";
        tout << "        }\n" ;
        tout << endl;
        tout << "    N_VDestroy_Serial(y0_);\n";
		if (options["sens"] == "yes"){
			tout << "    N_VDestroyVectorArray_Serial(yS, P_);\n";
		}
			
		tout << "    CVodeFree(&cvode_mem);\n";
        tout << "    }\n" ;
        tout.close();
		
        /* ============================================
		 * Create a Makefile for the CVODE demo program
		 * ============================================
         */
        string mfilename = "Makefile-"+Name()+"_cvsdemo";
        ofstream mout;
        mout.open(mfilename.c_str());
        mout << "#\n";
        mout << "# " << mfilename << endl;
        mout << "#\n";
        mout << "# This is the Makefile for the " << Name() << "_cvsdemo program.\n";
			// Default is CVODES v2.6.0
		mout << "# This file is configured for SUNDIALS v2.4.0, CVODES v2.6.0.\n";
        mout << "#\n";
        PrintVFGENComment(mout,"# ");
        mout << "#\n";
        mout << "# This Makefile is not guaranteed to work in all operating systems.\n";
        mout << "# You may have to edit this file to meet the conventions of your operating system.\n";
        mout << "#\n\n";
        mout << endl;
        mout << "SUNDIALS_DIR=/usr/local\n";
        mout << "SUNDIALS_LIB_DIR=$(SUNDIALS_DIR)/lib\n";
        mout << "SUNDIALS_INC_DIR=$(SUNDIALS_DIR)/include\n";
        // CVODE v2.4.0 and v2.6.0
		mout << "SUNDIALS_INCS=-I$(SUNDIALS_INC_DIR)\n";
        mout << "LIBS=-lm\n";
        mout << endl;
        mout << "all: " << Name() << "_cvsdemo" << endl;
        mout << endl;
        mout << Name() << "_cvsdemo: " << Name() << "_cvsdemo.o " << Name() << "_cvs.o" << endl;
        mout << "\t$(CC) $(LDFLAGS) -o " << Name() << "_cvsdemo ";
        mout <<       Name() << "_cvsdemo.o " << Name() << "_cvs.o -L$(SUNDIALS_LIB_DIR) $(SUNDIALS_LIBS) $(LIBS)" << endl;
        mout << endl;
        mout << Name() << "_cvsdemo.o: " << Name() << "_cvsdemo.c " << Name() << "_cvs.h" << endl;
        mout << "\t$(CC) $(CPPFLAGS) $(SUNDIALS_INCS) -c " << Name() << "_cvsdemo.c" << endl;
        mout << endl;
        mout << Name() << "_cvs.o: " << Name() << "_cvs.c " << Name() << "_cvs.h" << endl;
        mout << "\t$(CC) $(CPPFLAGS) $(SUNDIALS_INCS) -c " << Name() << "_cvs.c" << endl;
        mout << endl;
        mout << "clean:\n";
        mout << "\trm -f " << Name() << "_cvsdemo " << Name() << "_cvsdemo.o " << Name() << "_cvs.o" << endl;
        mout << endl;
        mout.close();
	}
}

