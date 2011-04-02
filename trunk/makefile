optimization = -O3 -fopenmp
objects = mvnorm.o CountData.o gzstream.o State.o BinaryTree.o Node.o Tree.o iterator.o
libs = -lgsl -lgslcblas -lm -lz
inc = -I/opt/local/include

test: test.o $(objects)
	g++-mp-4.3 $(libs) -o test test.o $(objects) $(optimization)

%.o: %.cpp
	g++-mp-4.3 -c $< -o $@ $(inc) $(optimization)

clean:
	rm *.o
