#include <mpi.h>
#include <random>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <chrono>
#include <cmath>
using namespace std;
using namespace std::chrono;
double* thera_grid;
double* phi_grid;
double alpha = 0.1;
double beta = 0.2;

int Nt = 10000; // need to x3
int Np = 10000; // need to x3
double Rmin = 1.0;
double Rmax = 2.0;
double pi = acos(-1.0);
void potTerms(double R, double t, double p, double V)
{
	auto R2 = pow(R,2);
	V = 1.0/R + cos(t)*exp(-R*alpha)/R2 + sin(t)*exp(-R*beta)*cos(p)/R2/R;
}
void linSpacing(double xmin, double xmax, double* array, int Len)
{
	double st = (xmax-xmin)/(Len-1.0);
	for(int i=0;i<Len;i++)
	{
		array[i] = xmin + st*double(i);
	};
} 
void prepareArrays()
{
	thera_grid = new double[Nt];
	phi_grid = new double[Np];
	linSpacing(0,2*pi,phi_grid,Np);
	linSpacing(0,pi,thera_grid,Nt);		
}
double computeSomeValue(int i, int j, double R )
{
	double value;
	potTerms(R,thera_grid[j],phi_grid[j], value);
	return value;
}
int main(){
   MPI_Init(NULL,NULL);
   int world_rank = 0; // just some definition for world_rank
   int world_size = 1;// just some definition for world_size
   int error1 = MPI_Comm_size(MPI_COMM_WORLD, &world_size);
   int error2 = MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
   double R  = Rmax;
   auto t1 = MPI_Wtime();

	if(world_size>1)
		 R = Rmin + (Rmax-Rmin)/(world_size-1.0)*world_rank;
	double** matrix = new double*[Np];
	for(int i=0; i<Np;i++)
		matrix[i] = new double[Nt];

	for(int i=0; i<Np; i++)
		for(int j=0;j<Nt;j+=2)
		{
			matrix[i][j] = computeSomeValue(i,j,R);
			matrix[i][j+1] = computeSomeValue(i,j+1,R); // compiler is smart enough to unroll the loop on its own
		}
    auto t2 = MPI_Wtime(); 
	cout<<"ID "<<world_rank << " TIME IS "<< t2 -t1<<endl;
	if(world_rank==0)
		cout<<"Buffer size is "<< double(Np*Nt)/1024/1024/1024*8.0 << "GB " <<endl;
	
   MPI_File fileToWrite;
   MPI_File singleFile;
   MPI_Status status;
   MPI_Offset filesize;
       auto t3 = MPI_Wtime(); 
   auto error = MPI_File_open(MPI_COMM_WORLD, "testfile", MPI_MODE_CREATE | MPI_MODE_WRONLY,
	MPI_INFO_NULL, &fileToWrite); 
//   MPI_File_set_view(fileToWrite, world_rank * Np*Nt * sizeof(double),
//	MPI_DOUBLE, MPI_DOUBLE, "native",MPI_INFO_NULL);
	MPI_File_write_at_all(fileToWrite, world_rank * Np*Nt * sizeof(double), matrix, Np*Nt, MPI_DOUBLE, &status);
//	MPI_File_write(fileToWrite, matrix, Np*Nt, MPI_DOUBLE,
//	MPI_STATUS_IGNORE);
    auto t4 = MPI_Wtime(); 
	if(world_rank==0)
		cout<<"ID "<<world_rank << " TIME TO WRITE COLLECTIVE IS "<< t4 -t2<<endl;
	MPI_File_close(&fileToWrite); 
	
	
	//testing single file
		error = MPI_File_open(MPI_COMM_WORLD, "testfileSingle", MPI_MODE_CREATE | MPI_MODE_WRONLY,
	MPI_INFO_NULL, &singleFile);
	       auto t5 = MPI_Wtime();
	if(world_rank==0)
		MPI_File_write(singleFile, matrix, Np*Nt, MPI_DOUBLE,
	MPI_STATUS_IGNORE);
		MPI_File_close(&singleFile); 
	    auto t6 = MPI_Wtime(); 
			if(world_rank==0)
			cout<<"ID "<<world_rank << " TIME TO WRITE SINGLE IS "<< t6 -t5<<endl;
		
   MPI_Barrier( MPI_COMM_WORLD);	
   MPI_Finalize();
   return 0;
}