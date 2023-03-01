all: compile clean run 

compile:
	gcc -c -Wall -g libs/control_panel.c
	gcc -c -Wall -g libs/preprocessor.c
	gcc -c -Wall -g libs/hash_table.c -lm
	gcc -c -Wall -g libs/DC.c -lgsl -lgslcblas
	gcc -c -Wall -g libs/sol_methods.c -lgsl -lgslcblas -lm
	gcc -c -Wall -g libs/csparse.c
	gcc -c -Wall -g libs/transient.c -lgsl -lgslcblas
	gcc -c -Wall -g libs/AC.c -lm -libs/CXSparse/Sparse/CXSparse/Lib/libcxsparse.a 
	gcc -c -Wall -g main.c
	gcc -o main main.o control_panel.o preprocessor.o hash_table.o DC.o AC.o sol_methods.o csparse.o transient.o -lm -lgsl -lgslcblas libs/CXSparse/Sparse/CXSparse/Lib/libcxsparse.a 

clean:
	rm *.o
	
run:
	./main tests/test_ac.txt
