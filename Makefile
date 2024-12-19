all:
	gcc -Wall --pedantic src/backsubst.c src/gauss_B.c src/main.c src/mat_io.c -o bin/gauss

test: all
	bin/gauss dane/A dane/b
test1: all 
	bin/gauss dane/A1 dane/b1
test2: all 
	bin/gauss dane/A2 dane/b2
test3: all 
	bin/gauss dane/A3 dane/b3
