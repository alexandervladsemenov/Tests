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

class brownianTrajectory
{
  private:
        double step, position, time, timeStep;
        const double Rmax = 3e4;
        int myid;
  public:
  brownianTrajectory(int myid)
  {
      this->myid = myid;
      mt19937 gen;
      uniform_real_distribution<double> dis(-1,1);
      gen.seed(myid);
      step = abs(dis(gen))+0.1; // avoid accidnetall very small step
      timeStep = 10.0;
      position = 0;
      while(abs(position)<Rmax)
      {
        double rand = dis(gen);
        position+=step*rand/abs(rand);
        time+=timeStep;
      }
   }
   double getTime()
   {
     return time;
   } 
  ~brownianTrajectory()
  {
    cout<<"Time is "<<time<<endl;
    cout<<"The position is "<<position<<endl;
    cout<<"Myid is "<< myid <<"\n"<<endl;
  }
  
};
template <typename T>
void makeItZero(T* y, int len)
{
	for(int i=0;i <len; i++)
		y[i] = 0;
}
int main(){
    MPI_Init(NULL,NULL);
    MPI_Win win; // window RMA
   int NSizeOfTrajectoryBank = 10;
   int* trajects = new int[NSizeOfTrajectoryBank]; // this will be shared
   double* times = new double[NSizeOfTrajectoryBank]; // this will be share
   double* localtimes = new double[NSizeOfTrajectoryBank]; // this will be share
   double* idTime = new double[1];
    int* localbuff = new int[NSizeOfTrajectoryBank];// local buffer
    // zeroing everything
    makeItZero(trajects,NSizeOfTrajectoryBank);
    makeItZero(localbuff,NSizeOfTrajectoryBank);
    makeItZero(times,NSizeOfTrajectoryBank);
    makeItZero(localtimes,NSizeOfTrajectoryBank);
    makeItZero(idTime,1);   
    // done
   int world_rank = 0; // just some definition for world_rank
   int world_size = 1;// just some definition for world_size
   int error1 = MPI_Comm_size(MPI_COMM_WORLD, &world_size);
   int error2 = MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
   int len; // processsor name
   char* procNames = new char(world_size);
   MPI_Get_processor_name(procNames, &len);
   MPI_Win_create(idTime, 1, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
   MPI_Win_fence(0, win);
   brownianTrajectory myTrajectory(world_rank);
   idTime[0] = myTrajectory.getTime();
   localtimes[world_rank] = idTime[0];
   cout<< " DONE IS " << world_rank << "\n" << endl;
   for(int i =0; i<world_size; i++)
   {
     if((world_rank!=i) && (world_rank==3))
     MPI_Get(&localtimes[i], 1, MPI_DOUBLE, i, 0, 1, MPI_DOUBLE, win);
        }   

     for(int i =0; i<world_size; i++)
     {  
       if(world_rank==3)
       {
       cout<<"Time in "<<i<< " CPU IS "<<localtimes[i] <<endl;
       cout<< "THIS WAS ID "<<world_rank<<"\n"<<endl;
       }
     }
   MPI_Win_fence(0, win);
   MPI_Barrier( MPI_COMM_WORLD);	
   MPI_Finalize();

   return 0;
}
