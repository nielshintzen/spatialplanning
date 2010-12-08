//----------------------------------------------------------------------------//
// Define class of individual & function for class ind // 
//----------------------------------------------------------------------------//


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

#define O_f       -1.340    // PMRN slope females // probalistic maturation reaction norm
#define O_m       -0.500    // PMRN slope males   //
#define U_M        14.60    // intercept PMRN by sex //
#define U_F        24.70    // intercept PMRN by sex //
#define OMEGA      2.275598 // PMRN vertical width //

unsigned long int id=0;

//----------------------------------------------------------------------------//
// Define class of individual //
//----------------------------------------------------------------------------//
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

//----------------------------------------------------------------------------//
// Define function of class individual
//----------------------------------------------------------------------------//

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
