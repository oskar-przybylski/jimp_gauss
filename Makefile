all:
	gcc -Wall --pedantic src/backsubst.c src/gauss_B.c src/main.c src/mat_io.c -o bin/gauss

2by2_a: all
	bin/gauss dane/A dane/b
2by2_b: all 
	bin/gauss dane/A1 dane/b1
2by2_c: all 
	bin/gauss dane/A2 dane/b2
3by3: all 
	bin/gauss dane/A3 dane/b3
div_by_zero: all 
	bin/gauss dane/A4 dane/b4
non_square: all 
	bin/gauss dane/A5 dane/b5
negative_dim: all 
	bin/gauss dane/A6 dane/b6
