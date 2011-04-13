optimization = -O3 -fopenmp
objects = CmdLine.o MCMC.o MCMC_params.o mvnorm.o CountData.o gzstream.o State.o BinaryTree.o Node.o Tree.o iterator.o
libs = -lgsl -lgslcblas -lm -lz
inc = -I/opt/local/include

phylopop: PhyloPop.o $(objects)
	g++-mp-4.3 $(libs) -o phylopop PhyloPop.o $(objects) $(optimization)

test: test.o $(objects)
	g++-mp-4.3 $(libs) -o test test.o $(objects) $(optimization)

test2: test2.o $(objects)
	g++-mp-4.3 $(libs) -o test test.o $(objects) $(optimization)

%.o: %.cpp
	g++-mp-4.3 -c $< -o $@ $(inc) $(optimization)

clean:
	rm *.o
