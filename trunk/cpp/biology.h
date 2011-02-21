#define M_B        1.5E-4   // baseline mortality (old) divided by 365 to go to days //
#define ETA       -0.280    // exponent of size dependent natural mortality //
#define F_MAX      0.0      // fishing mortality //
#define PI         0.5938   // fishing selectivity constant //
#define LAMBDAple  2.200    // mesh size selection factor //
#define LAMBDAsol  2.200    // mesh size selection factor //
#define S_mesh     7.000    // mesh size //

#define BORNWGHT   30       // ton for 1 000 000 individuals, because for one individual it is in gram // //Individuals get born at age 1//
#define PLUSGROUP  15       // age of plusgroup

using namespace std ;
typedef float    (*FTYPE)[X_MAX][Y_MAX];

//----------------------------------------------------------------------------//
// Define biological functions
//----------------------------------------------------------------------------//
void   move             (struct ind x[], int Indvs, int tofy, FTYPE temp) ; 
void   growth           (struct ind x[], int Indvs, double B, int tofy, FTYPE food, FTYPE temp, double f1[]) ;
void   age              (struct ind x[], int Indvs)            ;
void   maturation       (struct ind x[], int Indvs)           ;
void   mortality        (struct ind x[], double lambda, int Indvs, double B );
int    alive2front      (struct ind x[])                       ;
int    reproduction     (struct ind x[], double R1, double R2, int Indvs, double Bspawn, FTYPE temp);
void   larvalmortality  (struct ind x[], int Indvs, FTYPE larvmort);
double rand_sex         ();

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

void growth (struct ind x[], int Indvs, double B, int tofy, FTYPE food, FTYPE temp, double f1[]){
  for(int n = 0 ; n < Indvs ; n++){	  
     x[n].weight  = x[n].weight + x[n].growth(food[tofy][x[n].X][x[n].Y], temp[tofy][x[n].X][x[n].Y],f1[tofy]);
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
      x[n].age  = 52;                                                     // they are born as 1 year olds, so set age at 52  
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

double rand_sex()     {return rand()%RAND_MAX + 1 ; }                       // function rand ] 0,+ INF ] 
