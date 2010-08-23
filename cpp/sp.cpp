#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib> 

using namespace std ;

#define POPMAX 600000
#define T_MAX  5000

#define ALPHA      0.750
#define BETA       3.000
#define k          0.010    // Fulton factor //
#define O_f       -1.340    // PMRN slope //
#define O_m       -0.500    
#define R1ple      9.0E2    // stock recruitment parameters //
#define R2ple      3.0E3     // stock recruitment parameters //
#define R1sol      9.0E2    // stock recruitment parameters //
#define R2sol      3.0E3     // stock recruitment parameters //
#define OMEGA      2.275598 // PMRN vertical width //
#define M_B        0.055    // baseline mortality (old) //
#define b          0.600    
#define a          5.500    // potential genotypic energy acquisition // 
#define ETA       -0.250 
#define F_MAX      1.5      // fishing mortality //
#define PI         0.5938   // fishing selectivity constant //
#define LAMBDAple  2.200    // mesh size selection factor //
#define LAMBDAsol  2.200    // mesh size selection factor //

#define S_mesh     7.000    // mesh size //
#define EGGWGHT    4.2E-3   // ton for 1 000 000 individuals, because for one individual it is in gram // 

#define U_M        14.60    // intercept PMR by sex //
#define U_F        24.70    // intercept PMRN by sex //

class ind
{
      public:     
       double length() ;
       double Lp50()   ;
       double u()      ;
        
       int age ;
       char sex ;        // 1= male, 2=female //
       char stage ;      // stage1= immature, stage 2= mature, stage 3&4 : dead, stage3 = fishing, stage4 = natural//
       double weight ;
       double u_m ;     // intercept PMRN by sex //
       double u_f ;
      
       int X ;     // reproductive investment rate by sex //
       int Y ;     
};

double ind::u()
{  
   if (sex == 1) {return u_m ; }
   else          {return u_f ; }
}

double ind::length()
{ 
   return  pow(weight / k, (1/BETA)) ;   
}
double ind::Lp50()
{ 
   double O_p  ;  
   if (sex == 1) {O_p = O_m ; } 
   else          {O_p = O_f ; }
   return u() + O_p * age ;   
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

void   move          (struct ind x[], int Indvs) ; 
void   growth        (struct ind x[], int Indvs, double B) ;
void   age           (struct ind x[], int Indvs)            ;
void   maturation    (struct ind x[], int Indvs)           ;
void   mortality     (struct ind x[], double lambda, int Indvs, double B );
void   output        (struct ind x[], int t, int number);
int    alive2front   (struct ind x[])                       ;
int    reproduction  (struct ind x[], double R1, double R2, int Indvs, double Bspawn);
void   aggregate     (struct cell g[][30], struct ind x[], int Indvs); 
double rnorm         (double mu, double sd)                 ;      // function random 
double max           (double first, double second)          ;      // function max 
double min           (double first, double second)          ;      // function min 
double rand_sex      ();


int main (){   
 double Btotalple, Bnurseple, Bspawnple, Btotalsol, Bnursesol, Bspawnsol ;    /* biomass on nursery, total biomass */
 int aliveple = POPMAX, alivesol = POPMAX; 

 cell grid[30][30];

 /* INITIALISE INDIVIDUALS AT START, FIRST PLAICE, THEN SOLE */

  for(int i=0; i < POPMAX; i++) {
    if( i > POPMAX/2){
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
      ple[i].X      = 15 ;
      ple[i].Y      = 15 ;
  }
  
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
      sol[i].X      = 15 ;
      sol[i].Y      = 15 ;
  }
  

  /* START SIM */
  for(int t = 0; t < T_MAX; t++){    

    /* INITIALISE GRID IN EACH TIMESTEP */
 
    for(int xx=0; xx < 30; xx++) {
      for(int yy=0; yy < 30; yy++) {
        grid[xx][yy].effort= 0;
        grid[xx][yy].kg_sol1 = 0 ;
        grid[xx][yy].kg_sol2 = 0 ;
        grid[xx][yy].kg_sol3 = 0 ;
      }
    }

   /* CALCULATE TOTAL BIOMASS AND BIOMASS ON NURSERY FOR TWO SPECIES */
    Bnurseple = Btotalple = Bspawnple = Bnursesol = Btotalsol = Bspawnsol = 0;
        
    for(int n = 0 ; n < POPMAX ; n++) {
        if (ple[n].stage < 3 ) {
            Btotalple  += ple[n].weight ; 
            if (ple[n].stage < 2 ) Bnurseple += ple[n].weight;
        }        
    } 
    Bspawnple =  Btotalple - Bnurseple;
    
    for(int n = 0 ; n < POPMAX ; n++) {
        if (sol[n].stage < 3 ) {
            Btotalsol  += sol[n].weight ; 
            if (sol[n].stage < 2 ) Bnursesol += sol[n].weight;
        }        
    } 
    Bspawnsol =  Btotalsol - Bnursesol;
    
    cout<<  Bspawnple<< " "<<Bspawnsol << endl; 

    move       (ple, aliveple) ; 
    growth     (ple, aliveple, Bnurseple) ;                                    // Function of growth //   
    age        (ple, aliveple)    ;                                            // Function of ageing //
    maturation (ple, aliveple)   ;                                             // Function of maturation //
    mortality  (ple, LAMBDAple, aliveple ,Bnurseple ) ;                        // Function mortality //
    
    move       (sol, alivesol);
    growth     (sol, alivesol, Bnursesol) ;                                    // Function of growth //    
    age        (sol, alivesol)    ;                                            // Function of ageing //    
    maturation (sol, alivesol)   ;                                             // Function of maturation //
    mortality  (sol, LAMBDAsol, alivesol ,Bnursesol ) ;                        // Function mortality */
    
    cout<< Bspawnple<<endl;
    output     (ple,t, 40);
    cout<< Bspawnsol<<endl;
    output     (sol,t, 40);
  
    aliveple = alive2front (ple)       ;                                       // shuffle so that alives are in front*/
    aliveple = reproduction(ple, R1ple, R2ple, aliveple, Bspawnple);           // Function of reproduction */
    alivesol = alive2front (sol)       ;                                       // shuffle so that alives are in front*/
    alivesol = reproduction(sol, R1sol, R2sol, alivesol, Bspawnsol);           // Function of reproduction */
  }  
  return 0 ;  
} 

void age (struct ind x[], int Indvs){
  for(int n = 0 ; n < Indvs ; n++) {x[n].age += 1 ; }       
}     

void move (struct ind x[], int Indvs){ 
    
  for(int n = 0 ; n < Indvs ; n++){	  
    x[n].X += rnorm(0,2)    ;
    x[n].Y += rnorm(0,2)    ;        
  }                                                           // end for loop over individuals //
}  

void growth (struct ind x[], int Indvs, double B){
 
  for(int n = 0 ; n < Indvs ; n++)	  
    x[n].weight = pow( a/b  - (a/b - pow(x[n].weight, 1 - ALPHA))*exp(-b *( 1- ALPHA)), 1/(1-ALPHA))  ; 
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
      m_p     = pow(x[n].weight , ETA) ;                         // predation mortality due to foraging   
      M_tot   = M_B +  m_p  ;                                    // TOTAL natural mortality = background + predation mortality    
                                                                 // Fishing mortality NOTE SHOULD DEPEND ON LOCATION 
      F       = F_MAX / ( 1 + exp(- PI * (x[n].length() - LAMBDA * S_mesh))) ;  
      F_tot   = F;                                               // TOTAL fishing mortality                                 
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

int reproduction (struct ind x[], double R1, double R2, int Indvs, double SSB){

  int N_r = (int) ( R1* SSB / (R2 + SSB )) ; /*NOTE SSB IS SCALED */
  if (N_r > 0 ){
    for(int nu = Indvs; nu < (Indvs + N_r); nu++){
      x[nu].sex    = (int) ceil(2*((double)rand_sex()/((double)RAND_MAX)));
      x[nu].age    = 0;
      x[nu].stage  = 1;
      x[nu].weight = EGGWGHT;
      x[nu].u_m    = U_M ;			
      x[nu].u_f    = U_F ; 
      x[nu].X      = 15 ;
      x[nu].Y      = 15 ;
    }                                                                        // end creating individual   
  }                                                                          // end if more than 0 recruits
 return (Indvs + N_r) ;
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
                                                                                back to have staus 1 and 2 at beginning */
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
      alv++ ;
    }
  }
  for(int n = alv ; n < POPMAX ; n++) x[n].stage  = 4 ;                    // the rest is unimportant so give stage 4 
  
  return (alv);
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
