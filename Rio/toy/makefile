all: ToyMCGenerator

ToyMCGenerator: ToyMCGenerator.o
	clang++ `root-config --libs` -lMinuit ToyMCGenerator.o -o ToyMCGenerator

ToyMCGenerator.o: ToyMCGenerator.C
	clang++ -c -g -Wall `root-config --cflags` ToyMCGenerator.C -o ToyMCGenerator.o

clean:
	rm -rf *o ToyMCGenerator
