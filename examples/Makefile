all: test_complexe2 test_complexe test_dot test_dot_dyn test_copy test_iamax test_iamin

INC_DIR=../include

LIB_DIR=../lib

LIBST=-lmnblas
LIBDYN=-lmnblasdyn

OPTIONS_COMPIL  =-Wall -O2 -fPIC -I$(INC_DIR)

OPTIONS_LINK_STATIC  =   -fopenmp -L$(LIB_DIR) $(LIBST)
OPTIONS_LINK_DYN  = -fopenmp -L$(LIB_DIR) $(LIBDYN)

test_dot: test_dot.o flop.o $(LIB_DIR)/libmnblas.a
	gcc -o test_dot test_dot.o flop.o $(OPTIONS_LINK_STATIC)

flop.o: flop.c flop.h
	gcc $(OPTIONS_COMPIL) -c flop.c

test_dot_dyn: test_dot.o flop.o
	gcc -o test_dot_dyn flop.o test_dot.o $(OPTIONS_LINK_DYN)

test_dot.o: test_dot.c $(INC_DIR)/mnblas.h
	gcc $(OPTIONS_COMPIL) -c test_dot.c

test_complexe.o: test_complexe.c  $(INC_DIR)/complexe.h
	gcc $(OPTIONS_COMPIL) -c test_complexe.c

test_complexe: test_complexe.o flop.o
	gcc -o test_complexe test_complexe.o flop.o $(OPTIONS_LINK_STATIC)

test_complexe2.o: test_complexe2.c  $(INC_DIR)/complexe2.h
	gcc $(OPTIONS_COMPIL) -c test_complexe2.c

test_complexe2: test_complexe2.o flop.o
	gcc -o test_complexe2 test_complexe2.o flop.o $(OPTIONS_LINK_STATIC)

test_copy.o: test_copy.c  $(INC_DIR)/complexe.h
	gcc $(OPTIONS_COMPIL) -c test_copy.c

test_copy: test_copy.o flop.o
	gcc -o test_copy test_copy.o flop.o $(OPTIONS_LINK_STATIC)

test_iamax.o: test_iamax.c  $(INC_DIR)/mnblas.h $(INC_DIR)/complexe.h
	gcc $(OPTIONS_COMPIL) -c test_iamax.c

test_iamax: test_iamax.o flop.o $(LIB_DIR)/libmnblas.a
	gcc -o test_iamax test_iamax.o flop.o $(OPTIONS_LINK_STATIC)

test_iamin.o: test_iamin.c  $(INC_DIR)/mnblas.h $(INC_DIR)/complexe.h
	gcc $(OPTIONS_COMPIL) -c test_iamin.c

test_iamin: test_iamin.o flop.o $(LIB_DIR)/libmnblas.a
	gcc -o test_iamin test_iamin.o flop.o $(OPTIONS_LINK_STATIC)


clean:
	rm -f *.o test_dot test_dot_dyn test_complexe test_iamax test_iamin test_complexe2 *~
