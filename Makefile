CC=gcc
CPP=c++

CPPFLAGS = -Wall -O2
CFLAGS = -Wall -O2 -fgnu89-inline -march=native -std=c11 
LDFLAGS = 
IFLAGS = 
CBLAS = `pkg-config --cflags cblas` `pkg-config --libs cblas` `pkg-config --libs lapack`

#RMHMC_SRC=./src/mcmc/RMHMC.c ./src/mcmc/mcmc_kernel.c ./src/ode/ode_model.c ./src/app/read_cnf.c ./src/app/model_parameters_rmhmc.c 
MALA_SRC=./src/mcmc/smmala.c ./src/mcmc/mcmc_kernel.c ./src/ode/ode_model.c ./src/app/read_cnf.c ./src/mcmc/mv_norm.c ./src/app/model_parameters_smmala.c 

ODE_SOURCE = ./src/app/odeSolver_main.c ./src/ode/ode_model.c
VFGEN_SOURCE = `ls ./src/vfgen/*.cpp`

BIN=bin
.PHONY: all

all: vfgen bin/ode_smmala ODE*cvs.c ODE*.so


bin/ode_smmala: ./src/app/ode_smmala.c $(MALA_SRC)
	$(CC) -D_SMMALA $(CFLAGS) $(MALA_SRC) -o $(BIN)/ode_smmala $(CBLAS) -lm -lgsl -lsundials_cvodes -lsundials_nvecserial -ldl src/app/ode_smmala.c

# bin/ode_smmala: ./src/app/ode_smmala.c $(MALA_SRC)
# 	$(CC) -D_SMMALA $(CFLAGS) $(MALA_SRC) -o $(BIN)/ode_smmala -lm -lgsl -lgslcblas -lcblas -lsundials_cvodes -lsundials_nvecserial -ldl src/app/ode_smmala.c

#bin/ode_rmhmc: ./src/app/ode_rmhmc.c $(RMHMC_SRC)
#	$(CC) -D_RMHMC $(CFLAGS) $(RMHMC_SRC) -o $(BIN)/ode_rmhmc -lm -lgsl -lcblas -lsundials_cvodes -lsundials_nvecserial -ldl src/app/ode_rmhmc.c

vfgen: 
	$(CPP) $(CPPFLAGS) $(VFGEN_SOURCE) -o vfgen $(IFLAGS) -DVERSION=\"2.4.1\" $(LDFLAGS) -lcln -lginac -lmxml -lpthread

ODE%_cvs.c: ODE%.xml
	./vfgen cvodes:sens=yes,func=yes ODE%.xml

ODE%.so: ODE%_cvs.c
	$(CC) -shared -fPIC $(CFLAGS) -o $@  ODE$*_cvs.c


