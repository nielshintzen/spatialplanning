#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <cstdlib>

#define POPMAX 5000         // Numbers of individuals to start simulation with maximum on Geertcomputer = 15000000 but the flag -mcmodel=large// 
#define T_MAX  1000         // Maximum number of years that sim runs // 
#define T_STEP 10           // Number of times output is written to disk
#define A_MAX  780          // Number of timesteps output will be written to disk
#define P_WRITE 6000        // Maximum number of individuals written to disk
#define SPAREA  2           // Spawning area: 1 = English Channel, 2 = Duitse bocht

#define R1ple      4.5E5    // stock recruitment parameters, This means that max 900 individuals are born (being 900 million individuals in field), matches up with weights//
#define R2ple      3.0E3    // stock recruitment parameters //
#define R1sol      9.0E2    // stock recruitment parameters //
#define R2sol      3.0E3    // stock recruitment parameters //
#define MUT_RATE   0.050    // mutation rate //

unsigned long int id=0;
unsigned long int minid, maxid;

#include "class.h"          //Defines the class of the individual
#include "outinput.h"       //Functions to deal with input and output
#include "biology.h"        //Functions to deal with the biology as reproduction and maturation
#include "utils.h"          //Functions to get min, max etc

struct ind ple[POPMAX];
//struct ind sol[POPMAX];

int main (int argc, char* argv[]) {

  srand(atoi(argv[1])); //Take argument to write special extension to file
  double Btotalple, Bnurseple, Bspawnple, Btotalsol, Bnursesol, Bspawnsol ;    /* biomass on nursery, total biomass */

  //Read in the data
  readgrid(&GridFood , X_MAX, Y_MAX, 52, theFood);
  cout << "Read Food completed" << endl;

  readgrid(&GridTemp , X_MAX, Y_MAX, 52, theTemp);
  cout << "Read Temp completed" << endl;

  readgrid(&GridLMort , X_MAX, Y_MAX, 52, theLMort);
  cout << "Read Larval Mortality completed" << endl;

  readgrowthgam(&WeekPropFood,52,theGrowthGam);
  cout << "Read growth gam completed" << endl;

  /* INITIALISE INDIVIDUALS AT START, FIRST PLAICE, THEN SOLE */
  for(int i=0; i < POPMAX; i++) {
    ple[i].sex    = (i%2)+1;
    ple[i].weight = BORNWGHT;
    ple[i].id     = id ;
    ple[i].stage  = 1 ; /* everybody should be mature */
    ple[i].age    = 52 ;
    ple[i].u_m    = U_M ;
    ple[i].u_f    = U_F;
    if(SPAREA == 1){ 
      ple[i].X      = 75 ;
      ple[i].Y      = 53 ;
    } else if(SPAREA == 2){
        ple[i].X      = 91 ;
        ple[i].Y      = 67 ;
    }             
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
        ple[i].weight  = ple[i].weight  + ple[i].growth(theFood[(dd+6)%52][X][Y], theTemp[(dd+6)%52][X][Y],theGrowthGam[(dd+6)%52]);      
    }    
    for(int dd=0; dd <L_CHR2;  dd++){ //check juvenile strategy
      ple[i].adultXdir[dd] = (char)( (rand()% 11) -5); // Movement of maximum 5 left or right //
      ple[i].adultYdir[dd] = (char)( (rand()% 11) -5); // Movement of maximum 5 up or down    //
    }
    ple[i].weight = BORNWGHT;
    id++;
  }

//  for(int i=0; i < POPMAX; i++){
//    sol[i].sex    = (i%2)+1;
//    sol[i].weight = BORNWGHT ;
//    sol[i].id     = id ;
//    sol[i].stage  = 1 ;
//    sol[i].age    = 52 ;
//    sol[i].u_m    = U_M ;
//    sol[i].u_f    = U_F;
//    if(SPAREA == 1){ 
//      sol[i].X      = 75 ;
//      sol[i].Y      = 53 ;
//    } else if(SPAREA == 2){
//        sol[i].X      = 91 ;
//        sol[i].Y      = 67 ;
//    } 
//    int X         = sol[i].X;
//    int Y         = sol[i].Y;
//    int resX, resY;
//    for(int dd=0; dd < L_CHR1;  dd++){                         //check juvenile strategy
//        do{sol[i].juvXdir[dd] = (char)( (rand()% 11) -5);      // Movement of maximum 5 left or right //
//           sol[i].juvYdir[dd] = (char)( (rand()% 11) -5);      // Movement of maximum 5 up or down    //
//           resX = (int) (sol[i].juvXdir[dd] * sol[i].swim()) ;
//           resY = (int) (sol[i].juvYdir[dd] * sol[i].swim()) ;
//        } while ((theTemp[1][X + resX][Y + resY ] < -15) ||(( X + resX) <0) ||((X + resX) > X_MAX) ||((Y + resY )< 0) ||((Y + resY) > Y_MAX));
//        X += resX;
//        Y += resY;
//        sol[i].weight  = sol[i].weight  + sol[i].growth(theFood[(dd+6)%52][X][Y], theTemp[(dd+6)%52][X][Y],theGrowthGam[(dd+6)%52]);      
//    }
//    for(int dd=0; dd < L_CHR2 ;  dd++){                        //check juvenile strategy
//      sol[i].adultXdir[dd] = (char)((rand()% 11) -5);          // Movement of maximum 5 left or right //
//      sol[i].adultYdir[dd] = (char)((rand()% 11) -5);          // Movement of maximum 5 up or down    //
//    }
//    sol[i].weight = BORNWGHT;
//    id++;
//  }                                                            // end for loop over individuals /
  cout << "Initialisation of Plaice and Sole done" << endl;

  int aliveple = POPMAX, alivesol = POPMAX;

  string ext(".csv");                                          //Open file to write output to disk
  string SPname("_SPAREA");
  char buffer [4];
  string SP(itoa(SPAREA,buffer,10));
  filename += ( argv[1] + SPname + SP + ext);
  cout << filename << endl;
  popname += (argv[1] + SPname + SP + ext);
  myfile.open (filename.c_str() );
  mypopulation.open(popname.c_str());

  /* START SIM */
  for(int t = 6; t < T_MAX; t++){
   /* CALCULATE TOTAL BIOMASS AND BIOMASS ON NURSERY FOR TWO SPECIES */
    Bnurseple = Btotalple = Bspawnple = Bnursesol = Btotalsol = Bspawnsol = 0;

    for(int n = 0 ; n < aliveple ; n++) {
        if (ple[n].stage < 3 ) {
            Btotalple  += ple[n].weight ;
            if (ple[n].stage < 2 ) Bnurseple += ple[n].weight;
        }
    }

//    for(int n = 0 ; n < alivesol ; n++) {
//        if (sol[n].stage < 3 ) {
//            Btotalsol  += sol[n].weight ; 
//            if (sol[n].stage < 2 ) Bnursesol += sol[n].weight;
//        }
//    }

    Bspawnple =  Btotalple - Bnurseple;
//    Bspawnsol =  Btotalsol - Bnursesol;

    move(ple, aliveple, t%52, theTemp);                                                       // Move individuals every tenth timestep //
//    move(sol, alivesol, t%52, theTemp);                                                       // Move individuals every tenth timestep //  

    age        (ple, aliveple)    ;                                                           // Function of ageing //
//    age        (sol, alivesol)    ;                                                           // Function of ageing //    

    mortality(ple, LAMBDAple, aliveple ,Bnurseple ) ;                                         // Function mortality //
//    mortality(sol, LAMBDAsol, alivesol ,Bnursesol ) ;                                         // Function mortality */

    growth     (ple, aliveple, Bnurseple, t % 52, theFood, theTemp, theGrowthGam) ;                        // Function of growth //   
//    growth     (sol, alivesol, Bnursesol, t % 52, theFood, theTemp, theGrowthGam) ;                        // Function of growth //

    if(t%52 == 10 ) maturation (ple, aliveple)   ; //Checked with Cindy, gonads start to develop in March // Function of maturation //
//    if(t%52 == 10 ) maturation (sol, alivesol)   ;                                           // Function of maturation //

    if(t%52 == 5) cout<<"i " << argv[1]  <<" t " << t << " ssb ple " << Bspawnple<<" num ple "<<aliveple<< endl; //output(ple,t, 3);            // Write biomass and number to screen, followed by data for 10$
//    if(t%52 == 5) cout<<"i " << argv[1]  <<  " t " << t << " ssb sol " << Bspawnsol<<" num sol "<<alivesol<< endl; //output(sol,t, 3);            // Write biomass and number to screen, followed by data fo$

    aliveple = alive2front (ple)  ;                                                           // shuffle so that alives are in front*/
//    alivesol = alive2front (sol)  ;                                                           // shuffle so that alives are in front*/

    if(t%52 == 5) aliveple = reproduction(ple, R1ple, R2ple, aliveple, Bspawnple, theTemp); // Function of reproduction  in week 5*/
//    if(t%52 == 5) alivesol = reproduction(sol, R1sol, R2sol, alivesol, Bspawnsol, theTemp); // Function of reproduction  in week 5*/

    if(t%52 == 5){ larvalmortality (ple, aliveple, theLMort); aliveple = alive2front (ple);} // larvalmortality depends on field, now uniform field where everybody survives //
//    if(t%52 == 5){ larvalmortality (sol, alivesol, theLMort); alivesol = alive2front (sol);} // larvalmortality depends on field, now uniform field where everybody survives // 

    //Write output
    if ((t==6) ||( (t+A_MAX) % (int)(T_MAX/(T_STEP-1)) < 52 && t % 52 == 6)){
        int nn  = aliveple;
        int age = ple[nn].age;
        do{ age = ple[nn].age;
          nn--;
        } while ((nn > (aliveple - P_WRITE)) && (age <= 53));
        minid = ple[nn + 1].id;
        maxid = ple[aliveple - 1].id;
    } else if ((t < 6 + A_MAX) ||( (t + A_MAX)% (int)(T_MAX/(T_STEP-1)) < A_MAX +52)){
      for(int nn = 0; nn < aliveple; nn++){
        if(ple[nn].stage < 3 && (ple[nn].id > minid & ple[nn].id < maxid)){
          myfile <<t << "," <<       ple[nn].id          << "," << (int) ple[nn].sex          << "," <<       ple[nn].age                           << "," << (int) ple[nn].stage 
                     << "," <<       ple[nn].X           << "," <<       ple[nn].Y            << "," <<       ple[nn].weight       
                     << "," << (int) ple[nn].juvXdir[(int) (ple[nn].age)]                     << "," << (int) ple[nn].juvYdir[(int) (ple[nn].age)]  
                     << "," << (int) ple[nn].adultXdir[(t+1)%52]                              << "," << (int) ple[nn].adultYdir[(t+1)%52]           << endl;
        }
      }
    }

    //Write output every 15 years (cycle of complete new population)        
    if(t % (A_MAX) == 5){ writePopStruct(mypopulation, ple,aliveple,t);}

  } //end of timeloop

  myfile.close() ; mypopulation.close();
  return 0 ;
}




