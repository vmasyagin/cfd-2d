#include <stdlib.h>
#include <stdio.h>
#include <solver.h>
#include "global.h"
#include "mpi.h"
#include <float.h>
#include <iostream>
#include <algorithm>
#include <fstream>

#include "SolverZeidel.h"

int main(int argc, char** argv)
{

	/*FILE* fin = fopen("rae2822.dat", "r");
	FILE* fout = fopen("rae2822.out", "w");

	int n = 5;
	double x, y;
	for (int i = 0; i < 128; i++)
	{
		fscanf(fin, "%lf %lf", &x, &y);
		fprintf(fout, "Point( %3d) = {  %10.8f,  %10.8f, 0, lc3};\n", n, x, y);
		n++;
	}
	fclose(fin);
	fclose(fout);*/
//    setbuf(stdout, 0);

	//CSRMatrix mtr(4);
	//mtr.set(0, 0, 1);
	//mtr.set(0, 1, 1);
	//mtr.add(0, 0, 1);
	//mtr.set(1, 1, 1);
	//mtr.add(3, 3, 4);
	//mtr.add(2, 3, 9);
	//mtr.set(3, 1, 10);
	//mtr.assemble();

#ifdef _DEBUG
	_controlfp(~(_MCW_EM & (~_EM_INEXACT) & (~_EM_UNDERFLOW)), _MCW_EM);
#endif

	Parallel::init(&argc, &argv);

	hLog = fopen("task.log", "w"); // ��������� ���� ��� ������ ����; ������ printf(...) ���������� ������������ log(...)
	
	Method * method = Solver::initMethod( "task.xml" ); 
	
	Solver :: runMethod		( method );
	Solver :: destroyMethod	( method );

	fclose(hLog);

	Parallel::done();

	return 0;
}