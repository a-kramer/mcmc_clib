CC=gcc
CPP=c++

CPPFLAGS = -Wall -O2
CFLAGS = -Wall -O2 -fgnu89-inline -march=native
LDFLAGS = 
IFLAGS = 

RMHMC_SRC=./src/mcmc/RMHMC.c ./src/mcmc/mcmc_kernel.c ./src/ode/ode_model.c ./src/app/read_cnf.c ./src/app/model_parameters_rmhmc.c 
MALA_SRC=./src/mcmc/smmala.c ./src/mcmc/mcmc_kernel.c ./src/ode/ode_model.c ./src/app/read_cnf.c ./src/mcmc/mv_norm.c ./src/app/model_parameters_smmala.c 

ODE_SOURCE = ./src/app/odeSolver_main.c ./src/ode/ode_model.c
VFGEN_SOURCE = `ls ./src/vfgen/*.cpp`

BIN=bin

all: vfgen bin/ode_smmala bin/ode_rmhmc ODEmodel11S26P4U.so

bin/ode_smmala: ./src/app/ode_smmala.c $(MALA_SRC)
	$(CC) -D_SMMALA $(CFLAGS) $(MALA_SRC) -o $(BIN)/ode_smmala -lm -lgsl -lcblas -lsundials_cvodes -lsundials_nvecserial -ldl src/app/ode_smmala.c

bin/ode_rmhmc: ./src/app/ode_rmhmc.c $(RMHMC_SRC)
	$(CC) -D_RMHMC $(CFLAGS) $(RMHMC_SRC) -o $(BIN)/ode_rmhmc -lm -lgsl -lcblas -lsundials_cvodes -lsundials_nvecserial -ldl src/app/ode_rmhmc.c

vfgen: 
	$(CPP) $(CPPFLAGS) $(VFGEN_SOURCE) -o vfgen $(IFLAGS) -DVERSION=\"2.4.1\" $(LDFLAGS) -lcln -lginac -lmxml

ODEmodel11S26P4U_cvs.c: ODEmodel11S26P4U.xml
	./vfgen cvodes:sens=yes,func=yes ODEmodel11S26P4U.xml

ODEmodel11S26P4U.so: ODEmodel11S26P4U_cvs.c
	$(CC) -shared -fPIC $(CFLAGS) -o ODEmodel11S26P4U.so  ODEmodel11S26P4U_cvs.c


