all: Log_Likelihood_Fitter

Log_Likelihood_Fitter: Log_Likelihood_Fitter.o
	g++ `root-config --glibs` -lMinuit Log_Likelihood_Fitter.o -o Log_Likelihood_Fitter

Log_Likelihood_Fitter.o: Log_Likelihood_Fitter.C
	g++ -c -g  -Wall `root-config --cflags` Log_Likelihood_Fitter.C -o Log_Likelihood_Fitter.o

clean:
	rm -rf *o Log_Likelihood_Fitter
