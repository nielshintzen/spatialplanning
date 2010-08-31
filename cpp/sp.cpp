#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib> 

using namespace std ;

#define POPMAX 15000         // Numbers of individuals to start simulation with 600000 // 
#define T_MAX  900          // Maximum number of years that sim runs // 
#define X_MAX  144          // Max X dimension of Lorna map/grid 144//
#define Y_MAX  120          // Max X dimension of Lorna map/grid 144///

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

#define M_B        1.5E-4   // baseline mortality (old) divided by 365 to go to days //
#define ETA       -0.280    // exponent of size dependent natural mortality //
#define F_MAX      0.0      // fishing mortality //
#define PI         0.5938   // fishing selectivity constant //
#define LAMBDAple  2.200    // mesh size selection factor //
#define LAMBDAsol  2.200    // mesh size selection factor //
#define S_mesh     7.000    // mesh size //

#define EGGWGHT    4.2E-3   // ton for 1 000 000 individuals, because for one individual it is in gram // 


typedef float    (*FTYPE)[X_MAX][Y_MAX];

class ind
{
      public:     
       double length() ;
       double Lp50()   ;
       double u()      ;
        
       int    age ;       // age in days
       char   sex ;       // 1= male, 2=female //
       char   stage ;     // stage1= immature, stage 2= mature, stage 3&4 : dead, stage3 = fishing, stage4 = natural//
       double weight ;
       double u_m ;       // intercept PMRN by sex  (males)  //
       double u_f ;       // intercept PMRN by sex  (females //  
      
       int    X ;         // current Xpos //
       int    Y ;         // current Xpos //
       char   Xdir[550] ; // strategy, which direction does fish go, chosen every 10 days (to reduce lenght of vector) length= 365*15/10=550)  
       char   Ydir[550] ; // strategy, which direction does fish go, chosen every 10 days (to reduce lenght of vector) length= 365*15/10=550)    
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
   return u() + O_p * (age/365) ; // not sure if division by 365 is the most correct/convenient way, but I changed age to be days instead years)  
}


struct ind ple[POPMAX]; 
struct ind sol[POPMAX]; 

class cell
{
  public:
  double kg_sol()   ;
  double effort;
  double depth; 
  double kg_sol1;
  double kg_sol2;
  double kg_sol3;
};

double cell::kg_sol(){ 
  return kg_sol1 + kg_sol2 + kg_sol3; 
}

void   move                 (struct ind x[], int Indvs) ; 
void   growth               (struct ind x[], int Indvs, double B, int tofy, FTYPE food, FTYPE temp) ;
void   age                  (struct ind x[], int Indvs)            ;
void   maturation           (struct ind x[], int Indvs)           ;
void   mortality            (struct ind x[], double lambda, int Indvs, double B );
void   output               (struct ind x[], int t, int number);
int    alive2front          (struct ind x[])                       ;
int    reproduction         (struct ind x[], double R1, double R2, int Indvs, double Bspawn, FTYPE temp);
void   aggregate            (struct cell g[][30], struct ind x[], int Indvs); 
void   readgrid             (fstream * aFile, int anXmax, int anYmax, int anTmax, FTYPE agrid);
double rnorm                (double mu, double sd)                 ;      // function random 
double max                  (double first, double second)          ;      // function max 
double min                  (double first, double second)          ;      // function min 
double rand_sex             ();


int main (){   
 double Btotalple, Bnurseple, Bspawnple, Btotalsol, Bnursesol, Bspawnsol ;    /* biomass on nursery, total biomass */

 FTYPE theFood    = (FTYPE) malloc((size_t)sizeof(*theFood) * 365);
 FTYPE theTemp    = (FTYPE) malloc((size_t)sizeof(*theTemp) * 365);

 fstream GridFood ("N:\\Projecten\\SpatialPlanning\\svnjjp\\data\\food.dat", ios::in);
 fstream GridTemp ("N:\\Projecten\\SpatialPlanning\\svnjjp\\data\\temp.dat", ios::in);
 //fstream GridFood ("media/n/Projecten/SpatialPlanning/svnjjp/data/food.dat", ios::in);
 //fstream GridTemp ("media/n/Projecten/SpatialPlanning/svnjjp/data/temp.dat", ios::in);

 cell grid[X_MAX][Y_MAX];

 readgrid(&GridFood , X_MAX, Y_MAX, 365, theFood);
 cout << "thefood pos d16,x50,y60 " << theFood[15][49][59]<<" should be 0.888712 "  << endl;
 	
 readgrid(&GridTemp , X_MAX, Y_MAX, 365, theTemp);	
 cout << "theTemp pos d16,x50,y60 " << theTemp[15][49][59]<<" should be 7.362849 "  << endl;
 
 /* INITIALISE INDIVIDUALS AT START, FIRST PLAICE, THEN SOLE */

  for(int i=0; i < POPMAX; i++) {
    if( i > (POPMAX/2)){
      ple[i].sex    = 2;
      ple[i].weight = 80 ;
    } else {
      ple[i].sex    = 1;
      ple[i].weight = 60 ;
    }  
    ple[i].stage  = 1 ; /* everybody should be mature */
    ple[i].age    = 1 ;
    ple[i].u_m    = U_M ;
    ple[i].u_f    = U_F;
    ple[i].X      = 50 ;
    ple[i].Y      = 60 ;
    int X =ple[i].X ;
    int Y =ple[i].Y ; 
    for(int dd=0; dd < 550; dd++){    
        do{ple[i].Xdir[dd] = (char)( (rand()% 11) -5); // Movement of maximum 5 left or right //
           ple[i].Ydir[dd] = (char)( (rand()% 11) -5); // Movement of maximum 5 up or down    //
        } while (theTemp[1][X + (int) ple[i].Xdir[dd]][Y + (int) ple[i].Ydir[dd] ] < -15 ||( X + (int) ple[i].Xdir[dd]) <0 || X + (int) ple[i].Xdir[dd] > X_MAX || Y + (int) ple[i].Ydir[dd] < 0 || Y + (int) ple[i].Ydir[dd] > Y_MAX);
        X = X + (int) ple[i].Xdir[dd];
        Y = Y + (int) ple[i].Ydir[dd];
    }
  }
  cout << "init plepop done" << endl;

  for(int i=0; i < POPMAX; i++){
    if( i > POPMAX/2)    {
      sol[i].sex    = 2;
      sol[i].weight = 80 ;
    } else{
      sol[i].sex    = 1;
      sol[i].weight = 60 ;
    }  
    sol[i].stage  = 1 ; /* everybody should be mature */
    sol[i].age    = 1 ;
    sol[i].u_m    = U_M ;
    sol[i].u_f    = U_F;
    sol[i].X      = 50 ;
    sol[i].Y      = 60 ;
    int X =sol[i].X ;
    int Y =sol[i].Y ; 
    for(int dd=0; dd < 550; dd++){
        do{sol[i].Xdir[dd] = (char)( (rand()% 11) -5);
           sol[i].Ydir[dd] = (char)( (rand()% 11) -5);
        } while( theTemp[1][X + (int) sol[i].Xdir[dd]][Y + (int) sol[i].Ydir[dd] ] < -15  ||( X + (int) sol[i].Xdir[dd])<0 || (X + (int) sol[i].Xdir[dd]) > X_MAX || (Y + (int) sol[i].Ydir[dd]) < 0 ||( Y + (int) sol[i].Ydir[dd]) > Y_MAX);
        X = X + (int) sol[i].Xdir[dd];
        Y = Y + (int) sol[i].Ydir[dd];
    }
  }                                                           // end for loop over individuals /
  cout << "init solpop done" << endl;

  int aliveple = POPMAX, alivesol = POPMAX; 

  /* START SIM */
  for(int t = 0; t < T_MAX; t++){    

  
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

    if (t%10 == 0 ) move(ple, aliveple) ;                                            // Move individuals every tenth timestep //
    if (t%10 == 0 ) move(sol, alivesol);                                             // Move individuals every tenth timestep //  
 
    age        (ple, aliveple)    ;                                                  // Function of ageing //
    age        (sol, alivesol)    ;                                                  // Function of ageing //    
 
    mortality(ple, LAMBDAple, aliveple ,Bnurseple ) ;                                // Function mortality //
    mortality(sol, LAMBDAsol, alivesol ,Bnursesol ) ;                                // Function mortality */
    
    growth     (ple, aliveple, Bnurseple, t % 365, theFood, theTemp) ;               // Function of growth //   
    growth     (sol, alivesol, Bnursesol, t % 365, theFood, theTemp) ;               // Function of growth //    
 
    maturation (ple, aliveple)   ;                                                   // Function of maturation //
    maturation (sol, alivesol)   ;                                                   // Function of maturation //
    
    cout<<"ssb ple " << Bspawnple<<" num ple "<<aliveple<< endl; output(ple,t, 10);  // Write biomass and number to screen, followed by data for 10 indivuals // 
    cout<<"ssb sol " << Bspawnsol<<" num sol "<<alivesol<< endl; output(sol,t, 10);  // Write biomass and number to screen, followed by data for 10 indivuals //
  
    aliveple = alive2front (ple)       ;                                             // shuffle so that alives are in front*/
    alivesol = alive2front (sol)       ;                                             // shuffle so that alives are in front*/

    if(t%365 == 30) aliveple = reproduction(ple, R1ple, R2ple, aliveple, Bspawnple, theTemp); // Function of reproduction */
    if(t%365 == 30) alivesol = reproduction(sol, R1sol, R2sol, alivesol, Bspawnsol, theTemp); // Function of reproduction */
   
  }  
  
  
  ofstream myfile;
  myfile.open ("d:\\testoutputspat.csv" );
  
  double tw=0;
  for(int nn = 0 ; nn < aliveple ; nn++) {
     int X=50;
     int Y=60;
     if (ple[nn].weight > tw){
       tw=ple[nn].weight; 
       for(int dd=0; dd < 550; dd++){
         X = X + ple[nn].Xdir[dd];
         Y = Y + ple[nn].Ydir[dd];
         myfile << nn <<"," << (int) ple[nn].sex <<","<< dd <<"," << X<<","<<Y <<"," << ple[nn].weight <<   endl;
       }
     }  
  }
  myfile.close() ;
  return 0 ;  
} 

void age (struct ind x[], int Indvs){
  for(int n = 0 ; n < Indvs ; n++) {x[n].age += 1 ; }       
}     

void move (struct ind x[], int Indvs){ 
  for(int nn = 0 ; nn < Indvs ; nn++){	  
    x[nn].X += (int) x[nn].Xdir[(int) (x[nn].age/10)] ; 
    x[nn].Y += (int) x[nn].Ydir[(int) (x[nn].age/10)] ;         
  }                                                           // end for loop over individuals //
}  

void growth (struct ind x[], int Indvs, double B, int tofy, FTYPE food, FTYPE temp){
 
  for(int n = 0 ; n < Indvs ; n++){	  
     x[n].weight  = x[n].weight + ((k*((46.0*food[tofy][x[n].X][x[n].Y])/(foodh+(46.0*food[tofy][x[n].X][x[n].Y])))*(Pam*
                 (exp((TA/Tref)-(TA/(273+(temp[tofy][x[n].X][x[n].Y]+0.0))))*((1+exp((TAlow/Tref)-(TAlow/TL))+exp((TAhigh/(TH-(0.1* x[n].length() )))-(TAhigh/Tref)))
                 /(1+exp((TAlow/(273.0+(temp[tofy][x[n].X][x[n].Y]+0.0)))-(TAlow/TL))+exp((TAhigh/(TH-(0.1* x[n].length() )))-(TAhigh/(273.0+(temp[tofy][x[n].X][x[n].Y]+0)))))))))
                 *pow((pow((m* x[n].length() ),3)),0.666666) -((pm*(exp((TA/Tref)-(TA/(273+(temp[tofy][x[n].X][x[n].Y]+0.0))))))*(pow((m* x[n].length() ),3))))/
                 ((k*((46.0*food[tofy][x[n].X][x[n].Y])/(foodh+(46.0*food[tofy][x[n].X][x[n].Y])))*Em)+Eg);
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
                                                                 // Natural mortality 
      m_p     = (pow(x[n].weight , ETA)/365) ;                   // predation mortality due to foraging   
      M_tot   = M_B +  m_p  ;                                    // TOTAL natural mortality = background + predation mortality    
                                                                 // Fishing mortality NOTE SHOULD DEPEND ON LOCATION 
      F_tot       = F_MAX / ( 1 + exp(- PI * (x[n].length() - LAMBDA * S_mesh))) ;  
                                                                 // Probability total mortality (natural + fishing)      
      psurv   = exp(-(M_tot + F_tot)) ;                                       
      if(((double)rand()/((double)RAND_MAX+1)) > psurv)
      {
        if(((double)rand()/((double)RAND_MAX+1)) >  M_tot / (M_tot + F_tot)){x[n].stage = 3 ; } // dead by fishery
        else                                                                {x[n].stage = 4 ; } // dead by natural causes  
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
      x[nu].u_m  = x[inher].u_m;
      x[nu].u_f  = x[inher].u_f;
      for(int dd=0; dd < 550; dd++){
        x[nu].Xdir[dd] = x[inher].Xdir[dd];
        x[nu].Ydir[dd] = x[inher].Ydir[dd];
      }                                                                    
      x[nu].age    = 0;
      x[nu].stage  = 1;
      x[nu].weight = EGGWGHT;
    }                                                                        // end creating individual   
  }    
  return (up_nu) ;
}

void   aggregate (struct cell g[][30], struct ind x[], int Indvs){
  for(int n = 0 ; n < Indvs ; n++)  {
    if (x[n].stage <3)    {
      g[x[n].X][x[n].Y].kg_sol1++;   /*incomplete, just a placeholder, this should depend on sizeclass of animal and have weight times number  */
    }
  }
}

int alive2front(struct ind x[]) 
{
  int alv = 0 ;                                                              /* bring alive ind to front so recruits can go to 
                                                                                back to have status 1 and 2 at beginning */
  for(int n = 0 ; n < POPMAX ; n++) {
    if(x[n].stage <3) {
      x[alv].sex    = x[n].sex    ;
      x[alv].stage  = x[n].stage  ;
      x[alv].age    = x[n].age    ;
      x[alv].weight = x[n].weight ;
      x[alv].u_f    = x[n].u_f    ;
      x[alv].u_m    = x[n].u_m    ;     	
      x[alv].X      = x[n].X      ;
      x[alv].Y      = x[n].Y      ;
      for(int dd=0; dd < 550; dd++){
        x[alv].Xdir[dd] = x[n].Xdir[dd] ;
        x[alv].Ydir[dd] = x[n].Ydir[dd] ;
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
