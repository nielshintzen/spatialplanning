#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <cstdlib> 

using namespace std ;

#define POPMAX 5000         // Numbers of individuals to start simulation with maximum on Geertcomputer = 15000000 but the flag -mcmodel=large// 
#define X_MAX  144          // Max X dimension of Lorna map/grid 144//
#define Y_MAX  120          // Max X dimension of Lorna map/grid 144///

#define T_MAX  3000         // Maximum number of years that sim runs // 
#define T_STEP 10           // Number of times output is written to disk
#define A_MAX  780          // Number of timesteps output will be written to disk
#define P_WRITE 5000        // Maximum number of individuals written to disk

#define L_CHR1       157    // length of 1st chromosome (with juvenile strategy)     // 
#define L_CHR2        52    // length of 2nd chromosome (with adult seasonal strategy) //

#define k           0.85    // Lorna DEB parm //
#define Pam        390.0    // Lorna DEB parm //
#define pm          19.4    // Lorna DEB parm //
#define Em        2500.0    // Lorna DEB parm //
#define Eg        5600.0    // Lorna DEB parm //
#define TA        7000.0    // Lorna DEB parm //
#define Tref       283.0    // Lorna DEB parm //
#define TAlow    50000.0    // Lorna DEB parm //
#define TL         277.0    // Lorna DEB parm //
#define TAhigh  100000.0    // Lorna DEB parm //
#define TH         297.0    // Lorna DEB parm //
#define m          0.219    // Lorna DEB parm //
#define foodh   0.000069    // Lorna DEB parm //

#define BETA       3.000

#define O_f       -1.340    // PMRN slope females // probalistic maturation reaction norm
#define O_m       -0.500    // PMRN slope males   //
#define U_M        14.60    // intercept PMRN by sex //
#define U_F        24.70    // intercept PMRN by sex //
#define OMEGA      2.275598 // PMRN vertical width //

#define R1ple      9.0E2    // stock recruitment parameters, This means that max 900 individuals are born (being 900 million individuals in field), matches up with weights//
#define R2ple      3.0E3    // stock recruitment parameters //
#define R1sol      9.0E2    // stock recruitment parameters //
#define R2sol      3.0E3    // stock recruitment parameters //
#define MUT_RATE   0.050    // mutation rate //


#define M_B        1.5E-4   // baseline mortality (old) divided by 365 to go to days //
#define ETA       -0.280    // exponent of size dependent natural mortality //
#define F_MAX      0.0      // fishing mortality //
#define PI         0.5938   // fishing selectivity constant //
#define LAMBDAple  2.200    // mesh size selection factor //
#define LAMBDAsol  2.200    // mesh size selection factor //
#define S_mesh     7.000    // mesh size //

#define BORNWGHT   30       // ton for 1 000 000 individuals, because for one individual it is in gram // //Individuals get born at age 1//
#define PLUSGROUP  15       // age of plusgroup


typedef float    (*FTYPE)[X_MAX][Y_MAX];

class ind
{
      public:     
       double length() ;
       double Lp50()   ;
       double u()      ;
       double growth(double fd, double tmp) ;
       double swim()  ;
       
       unsigned long int id ;       // age in days        
       int    age ;       // age in days
       char   sex ;       // 1= male, 2=female //
       char   stage ;     // stage1= immature, stage 2= mature, stage 3&4 : dead, stage3 = fishing, stage4 = natural//
       double weight ;
       double u_m ;       // intercept PMRN by sex  (males)  //
       double u_f ;       // intercept PMRN by sex  (females //  
      
       int    X ;         // current Xpos //
       int    Y ;         // current Xpos //
       char   juvXdir[L_CHR1] ; // juvenile strategy, which direction does fish go, chosen every 7 days (to reduce lenght of vector) length= 365*3/7=157)  
       char   juvYdir[L_CHR1] ; // juvenile strategy, which direction does fish go, chosen every 7 days (to reduce lenght of vector) length= 365*3/7=157)    
       char   adultXdir[L_CHR2] ; // adult strategy, which direction does fish go, chosen every 7 days (to reduce lenght of vector) length= 365*1/7=52)  
       char   adultYdir[L_CHR2] ; // adult strategy, which direction does fish go, chosen every 7 days (to reduce lenght of vector) length= 365*1/7=52)    

}; 

double ind::u()
{  
   if (sex == 1) {return u_m ; }
   else          {return u_f ; }
}

double ind::length()
{ 
   return (pow(weight,0.3333333)/m) ;  /* This now comes from Lornas R file, TO BE FIXED, INSTEAD OF WEIGHT LORNA USES V (volume) */   
}
double ind::Lp50()
{ 
   double O_p  ;  
   if (sex == 1) {O_p = O_m ; } 
   else          {O_p = O_f ; }
   return u() + O_p * (age/52) ; // not sure if division by 52 (weeks) is the most correct/convenient way, but I changed age to be weeks instead years)  
}

double ind::growth(double fd, double tmp )
{ 
  return( 7 * ((k*((46.0*fd)/(foodh+(46.0*fd)))*(Pam*
                 (exp((TA/Tref)-(TA/(273+(tmp+0.0))))*((1+exp((TAlow/Tref)-(TAlow/TL))+exp((TAhigh/(TH-(0.1* length() )))-(TAhigh/Tref)))
                 /(1+exp((TAlow/(273.0+(tmp+0.0)))-(TAlow/TL))+exp((TAhigh/(TH-(0.1* length() )))-(TAhigh/(273.0+(tmp+0)))))))))
                 *pow((pow((m* length() ),3)),0.666666) -((pm*(exp((TA/Tref)-(TA/(273+(tmp+0.0))))))*(pow((m* length() ),3))))/
                 ((k*((46.0*fd)/(foodh+(46.0*fd)))*Em)+Eg));
}

double ind::swim()
{ 
  return(length()/20);
}


struct ind ple[POPMAX]; 
struct ind sol[POPMAX]; 

void   move             (struct ind x[], int Indvs, int tofy, FTYPE temp) ; 
void   growth           (struct ind x[], int Indvs, double B, int tofy, FTYPE food, FTYPE temp) ;
void   age              (struct ind x[], int Indvs)            ;
void   maturation       (struct ind x[], int Indvs)           ;
void   mortality        (struct ind x[], double lambda, int Indvs, double B );
void   output           (struct ind x[], int t,     int number);
void   popStruct        (ofstream &mypopulation, struct ind x[], int Indvs, int t);
int    writeOutput      (int t,  int write2file[]);
int    alive2front      (struct ind x[])                       ;
int    reproduction     (struct ind x[], double R1, double R2, int Indvs, double Bspawn, FTYPE temp);
void   larvalmortality  (struct ind x[], int Indvs, FTYPE larvmort);
void   readgrid         (fstream * aFile, int anXmax, int anYmax, int anTmax, FTYPE agrid);
double rnorm            (double mu, double sd)                 ;      // function random 
double max              (double first, double second)          ;      // function max 
double min              (double first, double second)          ;      // function min 
double rand_sex         ();

unsigned long int id=0;

int main (int argc, char* argv[]) {

  srand(atoi(argv[1]));
  double Btotalple, Bnurseple, Bspawnple, Btotalsol, Bnursesol, Bspawnsol ;    /* biomass on nursery, total biomass */

  FTYPE theFood    = (FTYPE) malloc((size_t)sizeof(*theFood)  * 52);
  FTYPE theTemp    = (FTYPE) malloc((size_t)sizeof(*theTemp)  * 52);
  FTYPE theLMort   = (FTYPE) malloc((size_t)sizeof(*theLMort) * 52);

  fstream GridFood  ("N:\\Projecten\\SpatialPlanning\\svnjjp\\data\\food7d.dat", ios::in);
  fstream GridTemp  ("N:\\Projecten\\SpatialPlanning\\svnjjp\\data\\temp7d.dat", ios::in);
  fstream GridLMort ("N:\\Projecten\\SpatialPlanning\\svnjjp\\data\\larvalmortality7d.dat", ios::in);
  
  ofstream mypopulation;
  mypopulation.open("d:\\mypopulation.csv");
  //fstream GridFood  ("d:\\data\\food7d.dat", ios::in);
  //fstream GridTemp  ("d:\\data\\temp7d.dat", ios::in);
  //fstream GridLMort ("d:\\data\\larvalmortality7d.dat", ios::in);
  //fstream GridFood ("/media/n/Projecten/SpatialPlanning/svnjjp/data/food7d.dat", ios::in);
  //fstream GridTemp ("/media/n/Projecten/SpatialPlanning/svnjjp/data/temp7d.dat", ios::in);
  //fstream GridLMort ("/media/n/Projecten/SpatialPlanning/svnjjp/data/larvalmortality7d.dat", ios::in);

  readgrid(&GridFood , X_MAX, Y_MAX, 52, theFood);
  cout << "thefood pos d16,x50,y60 " << theFood[15][49][59]<<" should be xxx "  << endl;
 	
  readgrid(&GridTemp , X_MAX, Y_MAX, 52, theTemp);	
  cout << "theTemp pos d16,x50,y60 " << theTemp[15][49][59]<<" should be xxx "  << endl;
 
  readgrid(&GridLMort , X_MAX, Y_MAX, 52, theLMort);	
  cout << "theLmort pos d16,x50,y60 " << theLMort[15][49][59]<<" should be xxx "  << endl;
  
  //Define sequence of timesteps to write to file
  int write2file[T_STEP] = {6};
  for(int Yr = 1; Yr < T_STEP; Yr++){ write2file[Yr] = ((int) (T_MAX / (T_STEP - 1))) + write2file[Yr -1];}

  /* INITIALISE INDIVIDUALS AT START, FIRST PLAICE, THEN SOLE */

  for(int i=0; i < POPMAX; i++) {
    ple[i].sex    = (i%2)+1;
    ple[i].weight = BORNWGHT;
    ple[i].id     = id ; 
    ple[i].stage  = 1 ; /* everybody should be mature */
    ple[i].age    = 0 ;
    ple[i].u_m    = U_M ;
    ple[i].u_f    = U_F;
    ple[i].X      = 75 ;
    ple[i].Y      = 53 ;
    int X         = ple[i].X;
    int Y         = ple[i].Y;
    int resX, resY;
    for(int dd=0; dd < L_CHR1;  dd++){ //check juvenile strategy    
        do{ple[i].juvXdir[dd] = (char)( (rand()% 11) -5); // Movement of maximum 5 left or right //
           ple[i].juvYdir[dd] = (char)( (rand()% 11) -5); // Movement of maximum 5 up or down    //
           resX = (int) (ple[i].juvXdir[dd] * ple[i].swim()) ;
           resY = (int) (ple[i].juvYdir[dd] * ple[i].swim()) ;            
        } while ((theTemp[1][X + resX][Y + resY ] < -15) ||(( X + resX) <0) || (( X + resX) > X_MAX) ||(( Y + resY) < 0) || ((Y + resY) > Y_MAX));
        X += resX;
        Y += resY;
        ple[i].weight  = ple[i].weight  + ple[i].growth(theFood[(dd+6)%52][X][Y], theTemp[(dd+6)%52][X][Y]);      
    }    
    for(int dd=0; dd <L_CHR2;  dd++){ //check juvenile strategy    
      ple[i].adultXdir[dd] = (char)( (rand()% 11) -5); // Movement of maximum 5 left or right //
      ple[i].adultYdir[dd] = (char)( (rand()% 11) -5); // Movement of maximum 5 up or down    //
    }
    ple[i].weight = BORNWGHT;
    id++;
  }
  
  for(int i=0; i < POPMAX; i++){
    sol[i].sex    = (i%2)+1;
    sol[i].weight = BORNWGHT ;
    sol[i].id     = id ; 
    sol[i].stage  = 1 ; /* everybody should be mature */
    sol[i].age    = 0 ;
    sol[i].u_m    = U_M ;
    sol[i].u_f    = U_F;
    sol[i].X      = 75 ;
    sol[i].Y      = 53 ;
    int X         = sol[i].X;
    int Y         = sol[i].Y;
    int resX, resY;
    for(int dd=0; dd < L_CHR1;  dd++){ //check juvenile strategy    
        do{sol[i].juvXdir[dd] = (char)( (rand()% 11) -5); // Movement of maximum 5 left or right //
           sol[i].juvYdir[dd] = (char)( (rand()% 11) -5); // Movement of maximum 5 up or down    //
           resX = (int) (sol[i].juvXdir[dd] * sol[i].swim()) ;
           resY = (int) (sol[i].juvYdir[dd] * sol[i].swim()) ;            
        } while ((theTemp[1][X + resX][Y + resY ] < -15) ||(( X + resX) <0) ||((X + resX) > X_MAX) ||((Y + resY )< 0) ||((Y + resY) > Y_MAX));
        X += resX;
        Y += resY;
        sol[i].weight  = sol[i].weight  + sol[i].growth(theFood[(dd+6)%52][X][Y], theTemp[(dd+6)%52][X][Y]);      
    }   
    for(int dd=0; dd < L_CHR2 ;  dd++){ //check juvenile strategy    
      sol[i].adultXdir[dd] = (char)((rand()% 11) -5); // Movement of maximum 5 left or right //
      sol[i].adultYdir[dd] = (char)((rand()% 11) -5); // Movement of maximum 5 up or down    //
    }   
    sol[i].weight = BORNWGHT; 
    id++;
  }                                                    // end for loop over individuals /
  cout << "init plepop and solpop done" << endl;

  int aliveple = POPMAX, alivesol = POPMAX; 

  ofstream myfile;
  string filename("D:\\testoutputspat_mut");
  string ext(".csv");
  filename += ( argv[1] + ext) ;
  myfile.open (filename.c_str() );
  

  /* START SIM */
  for(int t = 6; t < T_MAX; t++){    
    cout << "Start of simulation" << endl;
   /* CALCULATE TOTAL BIOMASS AND BIOMASS ON NURSERY FOR TWO SPECIES */
    Bnurseple = Btotalple = Bspawnple = Bnursesol = Btotalsol = Bspawnsol = 0;
        
    for(int n = 0 ; n < aliveple ; n++) {
        if (ple[n].stage < 3 ) {
            Btotalple  += ple[n].weight ; 
            if (ple[n].stage < 2 ) Bnurseple += ple[n].weight;
        }        
    } 
     
    for(int n = 0 ; n < alivesol ; n++) {
        if (sol[n].stage < 3 ) {
            Btotalsol  += sol[n].weight ; 
            if (sol[n].stage < 2 ) Bnursesol += sol[n].weight;
        }        
    } 
    
    Bspawnple =  Btotalple - Bnurseple;
    Bspawnsol =  Btotalsol - Bnursesol;

    move(ple, aliveple, t%52, theTemp);                                                       // Move individuals every tenth timestep //
    move(sol, alivesol, t%52, theTemp);                                                       // Move individuals every tenth timestep //  
 
    age        (ple, aliveple)    ;                                                           // Function of ageing //
    age        (sol, alivesol)    ;                                                           // Function of ageing //    

    mortality(ple, LAMBDAple, aliveple ,Bnurseple ) ;                                         // Function mortality //
    mortality(sol, LAMBDAsol, alivesol ,Bnursesol ) ;                                         // Function mortality */
    
    growth     (ple, aliveple, Bnurseple, t % 52, theFood, theTemp) ;                        // Function of growth //   
    growth     (sol, alivesol, Bnursesol, t % 52, theFood, theTemp) ;                        // Function of growth //    
   
    maturation (ple, aliveple)   ;                                                            // Function of maturation //
    maturation (sol, alivesol)   ;                                                            // Function of maturation //
    
    cout<<"ssb ple " << Bspawnple<<" num ple "<<aliveple<< endl; output(ple,t, 3);            // Write biomass and number to screen, followed by data for 10 indivuals // 
    cout<<"ssb sol " << Bspawnsol<<" num sol "<<alivesol<< endl; output(sol,t, 3);            // Write biomass and number to screen, followed by data for 10 indivuals //
  
    aliveple = alive2front (ple)  ;                                                           // shuffle so that alives are in front*/
    alivesol = alive2front (sol)  ;                                                           // shuffle so that alives are in front*/

    if(t%52 == 5) aliveple = reproduction(ple, R1ple, R2ple, aliveple, Bspawnple, theTemp); // Function of reproduction  in week 5*/
    if(t%52 == 5) alivesol = reproduction(sol, R1sol, R2sol, alivesol, Bspawnsol, theTemp); // Function of reproduction  in week 5*/

    if(t%52 == 5){ larvalmortality (ple, aliveple, theLMort); aliveple = alive2front (ple);} // larvalmortality depends on field, now uniform field where everybody survives //
    if(t%52 == 5){ larvalmortality (sol, alivesol, theLMort); alivesol = alive2front (sol);} // larvalmortality depends on field, now uniform field where everybody survives // 

    if ((t>428 && t <481)||( t > (T_MAX-53))){
      for ( int nn = 0 ; nn <aliveple; nn++){ 
        if(ple[nn].id > POPMAX){
          myfile <<t <<"," << nn << ","<< ple[nn].id <<"," << (int) ple[nn].sex <<"," <<ple[nn].age<< ","<<(int) ple[nn].stage << "," << ple[nn].X<<","<<ple[nn].Y <<"," << ple[nn].weight <<   endl;
        }
      }
    }

    //Write output, settings at top of file
    int idx;
    idx = writeOutput(t, write2file);
    int maxInd=1; vector<int> printInd(maxInd);
    if(idx > 0){                     //If you need to print output
      if(idx == 2){                  //If first point in time to print output
        int maxInd, mm = 0;
        for(int nn = 0; nn < aliveple; nn++){
          if(ple[n].age < 52){ maxInd++;}
        }
        vector<int> printInd(maxInd);
        for(int nn = 0; nn <aliveple; nn++){
          if(ple[nn].age < 52){ printInd[mm] = ple[nn].id;}
          mm++;
        }                               //Choose randomly a number of individuals
      }
        
      for(int nn = 0; nn < maxInd; nn++){
        int mm = printInd[nn];
        myfile <<t << "," <<       ple[mm].id          <<"," << (int) ple[mm].sex          <<"," <<       ple[mm].age          << "," << (int) ple[mm].stage 
                   << "," <<       ple[mm].X           <<"," <<       ple[mm].Y            <<"," <<       ple[mm].weight       << "," << (int) ple[mm].juvXdir[t] 
                   << "," << (int) ple[mm].juvYdir[t]  <<"," << (int) ple[mm].adultXdir[t] <<"," << (int) ple[mm].adultYdir[t] << endl;
      }   
    }
    
    //Write output every 15 years (cycle of complete new population        
    if(t % (52*PLUSGROUP) == 0){ popStruct(mypopulation, ple,aliveple,t);}
 
  } //end of timeloop
     
  myfile.close() ; mypopulation.close();
  return 0 ;  
} 

int writeOutput(int time, int file[]){
        int res = 0;
        for(int Yr = 0; Yr < T_STEP; Yr++){ 
          if(time == file[Yr]){ res = 2;}
          if((time - file[Yr]) < A_MAX & (time >= file[Yr])){ res = max(1,res);
          } else {res = max(0,res);}
        }
        return((int) res);
    }

void popStruct (ofstream &mypop, struct ind x[], int Indvs, int t){
    
    double cohort[PLUSGROUP] = {0};
    double weight[PLUSGROUP] = {0};
    double sex[PLUSGROUP]    = {0};
    double stage[PLUSGROUP]  = {0};
    for(int n = 0 ; n < Indvs ; n++){
      int age = ((int) x[n].age / 52);
      if(age > PLUSGROUP){ age = PLUSGROUP;}
      cohort[age]++;
      weight[age] += x[n].weight;
      if(x[n].sex == 1){ sex[age]++;}
      if(x[n].stage == 1){stage[age]++;}
    }
    
    for(int nn = 0; nn <= PLUSGROUP; nn++){
      mypop << t << "," << cohort[nn] << "," << weight[nn] << "," << sex[nn] << "," << stage[nn] << endl;
    }
}


void age (struct ind x[], int Indvs){
  for(int n = 0 ; n < Indvs ; n++) {x[n].age += 1 ; }       
}     

void move (struct ind x[], int Indvs, int tofy, FTYPE temp){ 
  int Xtemp, Ytemp;
  for(int nn = 0 ; nn < Indvs ; nn++){
    if (x[nn].stage < 2){             
      Xtemp = (int) (x[nn].juvXdir[(int) (x[nn].age)] * x[nn].swim()) ;
      Ytemp = (int) (x[nn].juvYdir[(int) (x[nn].age)] * x[nn].swim()) ;
    } else {
      Xtemp = (int) x[nn].adultXdir[(int) (tofy)];
      Ytemp = (int) x[nn].adultYdir[(int) (tofy)];
    }
    x[nn].X += ((temp[tofy][x[nn].X + Xtemp][x[nn].Y + Ytemp ] < -15 )||(x[nn].X + Xtemp < 0 )||(x[nn].X + Xtemp > X_MAX )) ? 0 : Xtemp ; 
    x[nn].Y += ((temp[tofy][x[nn].X + Xtemp][x[nn].Y + Ytemp ] < -15 )||(x[nn].Y + Ytemp < 0 )||(x[nn].Y + Ytemp > Y_MAX )) ? 0 : Ytemp ;                   
  }                                                           // end for loop over individuals //
}  

void growth (struct ind x[], int Indvs, double B, int tofy, FTYPE food, FTYPE temp){
 
  for(int n = 0 ; n < Indvs ; n++){	  
     x[n].weight  = x[n].weight + x[n].growth(food[tofy][x[n].X][x[n].Y], temp[tofy][x[n].X][x[n].Y]);
   }                
}      
    
void maturation (struct ind x[], int Indvs){  
  double pmat, Lp50 ;                                                         // probability of maturing
                                                                              // length when p=0.5 (50% individuals maturing
  for(int n = 0 ; n < Indvs ; n++) {
    if (x[n].stage < 2) { 
      pmat= 1/(1 + exp(-(x[n].length() - x[n].Lp50()) / OMEGA)) ;  
      if(((double)rand()/((double)RAND_MAX+1)) < pmat) x[n].stage  = 2;                                            
    }
  }
}

void mortality (struct ind x[], double LAMBDA, int Indvs, double B )  {
  double m_p, F, M_tot, F_tot, psurv ;
  
  for(int n = 0 ; n < Indvs ; n++)  {
    if (x[n].stage <3)    {
      if ( x[n].age > A_MAX){x[n].stage = 4 ;                       // Animals older than 15 yrs die (necessary becaus they dont have movement strategy anymore 
      } else{
             
        m_p     = (pow(x[n].weight , ETA)/52) ;                   // predation mortality due to foraging   
        M_tot   = M_B +  m_p  ;                                    // TOTAL natural mortality = background + predation mortality    
                                                                   // Fishing mortality NOTE SHOULD DEPEND ON LOCATION 
        F_tot       = F_MAX / ( 1 + exp(- PI * (x[n].length() - LAMBDA * S_mesh))) ;  
                                                                   // Probability total mortality (natural + fishing)      
        psurv   = exp(-(M_tot + F_tot)) ;                                       
        if(((double)rand()/((double)RAND_MAX+1)) > psurv ) 
        {
          if(((double)rand()/((double)RAND_MAX+1)) >  M_tot / (M_tot + F_tot)){x[n].stage = 3 ; } // dead by fishery
          else                                                                {x[n].stage = 4 ; } // dead by natural causes  
        }  
      }  
    }
  }
}

int reproduction (struct ind x[], double R1, double R2, int Indvs, double SSB, FTYPE temp){

  int up_nu, rndAdult1, rndAdult2, inher;
 
  int N_r = (int) ( R1* SSB / (R2 + SSB )) ; /*NOTE SSB IS SCALED */
  if (N_r > 0 ){ 
    up_nu = (int) min((Indvs + N_r),POPMAX);
    for(int nu = Indvs; nu < up_nu  ; nu++){
      //create new individuals 
      if(nu % 2 == 1){ //Males//
        x[nu].sex  = 1;
        do{ rndAdult1 = (int) rand() % Indvs;
            rndAdult2 = (int) rand() % Indvs;  
        } while ( (int) x[rndAdult1].sex !=1 || (int) x[rndAdult2].sex != 1);
        if(x[rndAdult1].weight >= x[rndAdult2].weight){inher = rndAdult1; } else {inher = rndAdult2;} 
      } else{ //Females//
        x[nu].sex  = 2;
        do{ rndAdult1 = (int) rand() % Indvs;
            rndAdult2 = (int) rand() % Indvs;  
        } while ((int) x[rndAdult1].sex !=2 || (int) x[rndAdult2].sex != 2);
        if(x[rndAdult1].weight >= x[rndAdult2].weight){inher = rndAdult1; } else {inher = rndAdult2;} 
      }  // end males or females
      x[nu].id   = id;
      x[nu].u_m  = x[inher].u_m;
      x[nu].u_f  = x[inher].u_f;
      x[nu].X    = x[inher].X;
      x[nu].Y    = x[inher].Y;
      for(int dd=0; dd < L_CHR1;  dd++){ //check juvenile strategy    
        x[nu].juvXdir[dd] = x[inher].juvXdir[dd];
        x[nu].juvYdir[dd] = x[inher].juvYdir[dd];
      }
      for(int dd=0; dd < L_CHR2;  dd++){ //check juvenile strategy    
        x[nu].adultXdir[dd] = x[inher].adultXdir[dd];
        x[nu].adultYdir[dd] = x[inher].adultYdir[dd];
      } 
      //mutate movement strategies with mutation rate of 10%
      if((double)rand()/((double)RAND_MAX+1) < MUT_RATE){
         int mpos = rand()% ((L_CHR1 + L_CHR2) + 1);
         if  (mpos < L_CHR1){
         x[nu].juvXdir[mpos] = (char)((rand()% 11) -5); 
         x[nu].juvYdir[mpos] = (char)((rand()% 11) -5); 
        } else { 
         x[nu].adultXdir[mpos-L_CHR1] = (char)((rand()% 11) -5); 
         x[nu].adultYdir[mpos-L_CHR1] = (char)((rand()% 11) -5);
        }
      }
      x[nu].age    = 0;
      x[nu].stage  = 1;
      x[nu].weight = BORNWGHT;
      
      id++;
    }                     // end creating individual   
  }    
  return (up_nu) ;
}

void larvalmortality (struct ind x[], int Indvs, FTYPE larvmort  )  {
  double psurv ;
  
  for(int n = 0 ; n < Indvs ; n++)  {
    if (x[n].age ==0)    {                                                 // only individuals of 0 days old suffer from larval mortality 
      psurv   = exp(-larvmort[5][x[n].X][x[n].Y]) ;                                       
      if(((double)rand()/((double)RAND_MAX+1)) > psurv){x[n].stage = 4 ; } // dead by natural causes  
     
      x[n].X    = 75;                                                     // now that we know which larvae died because of wrong position, move to standard position // 
      x[n].Y    = 53;  
    }
  }
}

int alive2front(struct ind x[]) 
{
  int alv = 0 ;                                                              /* bring alive ind to front so recruits can go to 
                                                                                back to have status 1 and 2 at beginning */
  for(int n = 0 ; n < POPMAX ; n++) {
    if(x[n].stage <3) {
      x[alv].id     = x[n].id    ;
      x[alv].sex    = x[n].sex    ;
      x[alv].stage  = x[n].stage  ;
      x[alv].age    = x[n].age    ;
      x[alv].weight = x[n].weight ;
      x[alv].u_f    = x[n].u_f    ;
      x[alv].u_m    = x[n].u_m    ;     	
      x[alv].X      = x[n].X      ;
      x[alv].Y      = x[n].Y      ;
      for(int dd=0; dd < L_CHR1 ; dd++){
        x[alv].juvXdir[dd] = x[n].juvXdir[dd] ;
        x[alv].juvYdir[dd] = x[n].juvYdir[dd] ;
      } 
      for(int dd=0; dd < L_CHR2; dd++){
        x[alv].adultXdir[dd] = x[n].adultXdir[dd] ;
        x[alv].adultYdir[dd] = x[n].adultYdir[dd] ;
      } 
      alv++ ;
    }
  }
  for(int n = alv ; n < POPMAX ; n++) x[n].stage  = 4 ;                    // the rest is unimportant so give stage 4 
  return (alv);
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

void output(struct ind x[],int t, int number){
  for(int n = 0 ; n < number ; n++) {
     cout <<t<<","<<(int) x[n].sex <<","<<(int) x[n].stage<<","<<(int) x[n].age<<","<<x[n].weight<<","<<x[n].u_f<<","<<x[n].u_m <<","<< x[n].X <<","<<x[n].Y << endl;
  }
}
double rnorm(double mu, double sd)                                          // function random for evolution 
{
  double Z1= (double) rand()/RAND_MAX ;
  double Z2= (double) rand()/RAND_MAX ;
  return(max(min( (sin(2.0*3.141592654*Z1)*sqrt(-2.0*log(Z2))*sd + mu), mu+4*sd), mu-4*sd)) ; // ugly fix to restrain max output to mu +sd*4 
}        
 
double max(double first, double second)                                     // fonction max 
{
  if (first < second) {return (second) ; }
  else                {return (first)  ; }
}

double min(double first, double second)                                     // fonction min 
{
  if (first > second) {return (second) ; }
  else                {return (first)  ; }
}

double rand_sex()     {return rand()%RAND_MAX + 1 ; }                       // function rand ] 0,+ INF ] 
