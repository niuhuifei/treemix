optimization = -O3 -fopenmp
objects = PhyloPop_params.o GraphState2.o GraphState3.o CountData2.o PopGraph.o CmdLine.o mvnorm.o CountData.o gzstream.o
libs = -L/opt/local/lib -lgsl -lgslcblas -lm -lz
inc = -I/opt/local/include

phylopop_g: PhyloPop_graph.o $(objects)
	g++-mp-4.3 $(libs) -o phylopop_g PhyloPop_graph.o $(objects) $(optimization)

treemix: TreeMix.o $(objects)
	g++ $(libs) -o treemix TreeMix.o $(objects) $(optimization)

fst: mean_fst.o $(objects)
	g++-mp-4.3 $(libs) -o fst mean_fst.o gzstream.o CmdLine.o CountData.o $(optimization)

test: test.o $(objects)
	g++-mp-4.3 $(libs) -o test test.o $(objects) $(optimization)

test2: test2.o $(objects)
	g++-mp-4.3 $(libs) -o test2 test2.o $(objects) $(optimization)

%.o: %.cpp
	g++ -c $< -o $@ $(inc) $(optimization)

clean:
	rm *.o
