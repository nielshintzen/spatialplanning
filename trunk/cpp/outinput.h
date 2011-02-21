#define X_MAX     144   // Max X dimension of Lorna map//
#define Y_MAX     120   // Max X dimension of Lorna map//
#define PLUSGROUP  15

using namespace std ;
typedef float    (*FTYPE)[X_MAX][Y_MAX];

//----------------------------------------------------------------------------//
// Reading in- and output files // 
//----------------------------------------------------------------------------//
FTYPE theFood    = (FTYPE) malloc((size_t)sizeof(*theFood)  * 52);
FTYPE theTemp    = (FTYPE) malloc((size_t)sizeof(*theTemp)  * 52);
FTYPE theLMort   = (FTYPE) malloc((size_t)sizeof(*theLMort) * 52);
double theGrowthGam[52];

#ifdef __linux__
  //OS = Linux
  fstream GridFood     ("/media/n/Projecten/SpatialPlanning/svnjjp/data/food7d.dat", ios::in);
  fstream GridTemp     ("/media/n/Projecten/SpatialPlanning/svnjjp/data/temp7d.dat", ios::in);
  fstream GridLMort    ("/media/n/Projecten/SpatialPlanning/svnjjp/data/larvalmortality7d.dat", ios::in);
  fstream WeekPropFood ("/media/n/Projecten/SpatialPlanning/svnjjp/data/growthgam7d.dat", ios::in);                                                          
  ofstream mypopulation;
  string popname       ("~/mypopulation_mut"); 
  ofstream myfile;
  string filename      ("~/testoutputspat_mut");  
#else 
  //OS = Windows
  fstream GridFood     ("N:\\Projecten\\SpatialPlanning\\svnjjp\\data\\food7d.dat", ios::in);
  fstream GridTemp     ("N:\\Projecten\\SpatialPlanning\\svnjjp\\data\\temp7d.dat", ios::in);
  fstream GridLMort    ("N:\\Projecten\\SpatialPlanning\\svnjjp\\data\\larvalmortality7d.dat", ios::in);
  fstream WeekPropFood ("N:\\Projecten\\SpatialPlanning\\svnjjp\\data\\growthgam7d.dat", ios::in);
  ofstream mypopulation;
  string popname       ("D:\\mypopulation");
  ofstream myfile;
  string filename      ("D:\\testoutputspat_mut");
#endif  

//----------------------------------------------------------------------------//
// Functions to deal with output //
//----------------------------------------------------------------------------//
void  writePopStruct (ofstream &mypopulation, struct ind x[], int Indvs, int time);
void  output         (struct ind x[],int t, int number);
void  readgrid       (fstream * aFile, int anXmax, int anYmax, int anTmax, FTYPE agrid);
void  readgrowthgam  (fstream * aFile, int weekMax,double agrowthgam[]);

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

void readgrid (fstream * aFile, int anXmax, int anYmax, int anTmax, FTYPE agrid){
   for (int tt = 0; tt < anTmax; tt++){
  	 for (int yy = 0; yy < anYmax; yy++){
	     for (int xx = 0; xx < anXmax; xx++){
 				  *aFile >> (agrid[tt][xx][yy]); 
       }
     }
	}		   
}

void readgrowthgam (fstream * aFile, int weekMax,double agrowthgam[]){
     for(int tt = 0; tt < weekMax; tt++){
       *aFile >> (agrowthgam[tt]);
     }
}
             


