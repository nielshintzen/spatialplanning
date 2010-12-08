#include <fstream>

#define X_MAX  144          // Max X dimension of Lorna map/grid 144//
#define Y_MAX  120          // Max X dimension of Lorna map/grid 144///

#define PLUSGROUP  15       // age of plusgroup

using namespace std ;
typedef float    (*FTYPE)[X_MAX][Y_MAX];

//----------------------------------------------------------------------------//
// Reading in- and output files // 
//----------------------------------------------------------------------------//

FTYPE theFood    = (FTYPE) malloc((size_t)sizeof(*theFood)  * 52);
FTYPE theTemp    = (FTYPE) malloc((size_t)sizeof(*theTemp)  * 52);
FTYPE theLMort   = (FTYPE) malloc((size_t)sizeof(*theLMort) * 52);



#ifdef __linux__
  //OS = Linux
  fstream GridFood     ("/media/n/Projecten/SpatialPlanning/svnjjp/data/food7d.dat", ios::in);
  fstream GridTemp     ("/media/n/Projecten/SpatialPlanning/svnjjp/data/temp7d.dat", ios::in);
  fstream GridLMort    ("/media/n/Projecten/SpatialPlanning/svnjjp/data/larvalmortality7d.dat", ios::in);
  ofstream mypopulation("~/mypopulation.csv", ios::out);
  ofstream myfile;
  string filename      ("~/testoutputspat_mut");  
#else 
  //OS = Windows
  fstream GridFood     ("N:\\Projecten\\SpatialPlanning\\svnjjp\\data\\food7d.dat", ios::in);
  fstream GridTemp     ("N:\\Projecten\\SpatialPlanning\\svnjjp\\data\\temp7d.dat", ios::in);
  fstream GridLMort    ("N:\\Projecten\\SpatialPlanning\\svnjjp\\data\\larvalmortality7d.dat", ios::in);
  ofstream mypopulation("D:\\mypopulation.csv", ios::out);
  ofstream myfile;
  string filename      ("D:\\testoutputspat_mut");
#endif  

//----------------------------------------------------------------------------//
// Functions to deal with output //
//----------------------------------------------------------------------------//

void          output           (struct ind x[], int t,     int number);
void          writePopStruct   (ofstream &mypopulation, struct ind x[], int Indvs, int time);
int           writeOutput      (int t,  int write2file[]);
void          readgrid         (fstream * aFile, int anXmax, int anYmax, int anTmax, FTYPE agrid);

int writeOutput(int time, int file[]){
        int res = 0;
        for(int Yr = 0; Yr < T_STEP; Yr++){ 
          if(time == file[Yr]){ res = 2;}
          if((time - file[Yr]) < A_MAX & (time >= file[Yr])){ res = max(1,res);
          } else {res = max(0,res);}
        }
        return((int) res);
    }

void writePopStruct (ofstream &mypop, struct ind x[], int Indvs, int time){
    
    double cohort[PLUSGROUP] = {0};
    double weight[PLUSGROUP] = {0};
    double sex[PLUSGROUP]    = {0};
    double stage[PLUSGROUP]  = {0};
    for(int n = 0 ; n < Indvs ; n++){
      int age = (int)(x[n].age / 52);
      if(age > PLUSGROUP){ age = PLUSGROUP;}
      cohort[age]++;
      weight[age] += x[n].weight;
      if(x[n].sex == 1){ sex[age]++;}
      if(x[n].stage == 1){stage[age]++;}
    }
    
    for(int nn = 0; nn < PLUSGROUP; nn++){
      mypop << time << "," << cohort[nn] << "," << weight[nn] << "," << sex[nn] << "," << stage[nn] << endl;
    }
}

void output(struct ind x[],int t, int number){
  for(int n = 0 ; n < number ; n++) {
     cout <<t<<","<<(int) x[n].sex <<","<<(int) x[n].stage<<","<<(int) x[n].age<<","<<x[n].weight<<","<<x[n].u_f<<","<<x[n].u_m <<","<< x[n].X <<","<<x[n].Y << endl;
  }
}

void readgrid (fstream * aFile, int anXmax, int anYmax, int anTmax, FTYPE agrid)
{
   for (int tt = 0; tt < anTmax; tt++){
  	 for (int yy = 0; yy < anYmax; yy++){
	     for (int xx = 0; xx < anXmax; xx++){
 				  *aFile >> (agrid[tt][xx][yy]); 
       }
      }
		}		   
    cout << ("Read grid finished") << endl << flush;
}
