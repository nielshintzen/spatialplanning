#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <cstdlib>

#define OS 1                // OS =1: windows, OS = 2: linux
#define POPMAX 10000        // Numbers of individuals to start simulation with maximum on Geertcomputer = 15000000 but the flag -mcmodel=large// 
#define T_MAX  4000         // Maximum number of years that sim runs // 
#define T_STEP 5            // Number of times output is written to disk
#define A_MAX  780          // Number of timesteps output will be written to disk
#define P_WRITE 5000        // Maximum number of individuals written to disk

#define R1ple      9.0E3    // stock recruitment parameters, This means that max 900 individuals are born (being 900 million individuals in field), matches up with weights//
#define R2ple      3.0E3    // stock recruitment parameters //
#define R1sol      9.0E2    // stock recruitment parameters //
#define R2sol      3.0E3    // stock recruitment parameters //
#define MUT_RATE   0.050    // mutation rate //

#include "class.h"          //Defines the class of the individual  
#include "outinput.h"       //Functions to deal with input and output
#include "biology.h"        //Functions to deal with the biology as reproduction and maturation
#include "utils.h"          //Functions to get min, max etc


struct ind ple[POPMAX]; 
struct ind sol[POPMAX]; 

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
  
  //Define sequence of timesteps to write to file
  int write2file[T_STEP] = {6};
  for(int Yr = 1; Yr < T_STEP; Yr++){ write2file[Yr] = ((int) ((T_MAX - 6) / (T_STEP - 1))) + write2file[Yr -1];}
  write2file[T_STEP]     = write2file[T_STEP] - A_MAX;  //Starting writing to file at end of step makes no sense therefore substract A_MAX
  int maxInd             = 0; 
  vector<int> printInd;

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
  cout << "Initialisation of Plaice and Sole done" << endl;

  int aliveple = POPMAX, alivesol = POPMAX; 

  //Open file to write output to disk
  string ext(".csv");
  filename += ( argv[1] + ext) ;
  myfile.open (filename.c_str() );
  //mypopulation.open();
  
  
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
    
    cout<<"t " << t << " ssb ple " << Bspawnple<<" num ple "<<aliveple<< endl; //output(ple,t, 3);            // Write biomass and number to screen, followed by data for 10 indivuals // 
    cout<<"t " << t << " ssb sol " << Bspawnsol<<" num sol "<<alivesol<< endl; //output(sol,t, 3);            // Write biomass and number to screen, followed by data for 10 indivuals //
  
    aliveple = alive2front (ple)  ;                                                           // shuffle so that alives are in front*/
    alivesol = alive2front (sol)  ;                                                           // shuffle so that alives are in front*/

    if(t%52 == 5) aliveple = reproduction(ple, R1ple, R2ple, aliveple, Bspawnple, theTemp); // Function of reproduction  in week 5*/
    if(t%52 == 5) alivesol = reproduction(sol, R1sol, R2sol, alivesol, Bspawnsol, theTemp); // Function of reproduction  in week 5*/

    if(t%52 == 5){ larvalmortality (ple, aliveple, theLMort); aliveple = alive2front (ple);} // larvalmortality depends on field, now uniform field where everybody survives //
    if(t%52 == 5){ larvalmortality (sol, alivesol, theLMort); alivesol = alive2front (sol);} // larvalmortality depends on field, now uniform field where everybody survives // 

    //Write output
    int idx =0;
    idx = writeOutput(t, write2file);
    if(idx > 0){                     //If you need to print output
      if(idx == 2){                  //If first point in time to print output
        for(int nn = 0; nn < maxInd; nn++) printInd.pop_back(); //Reset length of printInd to 0           
        maxInd = 0;                  //Reset maxInd to 0;
        for(int nn = 0; nn < aliveple; nn++){
          if(ple[nn].age < 52){ maxInd++;}
        }
        maxInd = min(P_WRITE,maxInd);//Only take max P_WRITE individuals
        int mm = 0;
        for(int nn = 0; nn < aliveple; nn++){
          if((int) ple[nn].age < 52 & (int) mm < (int) P_WRITE){
            printInd.push_back(mm);   //Extend length of vector
            printInd[mm] = ple[nn].id;
            mm++;
          }
        }                         
      }
      for(int nn = 0; nn < maxInd; nn++){
        int mm = (aliveple + 1);
        for(int p = 0; p < aliveple; p++){ if((int) ple[p].id == (int) printInd[nn]){ mm = p;}}
        if(ple[mm].stage < 3){
          myfile <<t << "," <<       ple[mm].id          << "," << (int) ple[mm].sex          << "," <<       ple[mm].age                           << "," << (int) ple[mm].stage 
                     << "," <<       ple[mm].X           << "," <<       ple[mm].Y            << "," <<       ple[mm].weight       
                     << "," << (int) ple[mm].juvXdir[(int) (ple[mm].age)]                     << "," << (int) ple[mm].juvYdir[(int) (ple[mm].age)]  
                     << "," << (int) ple[mm].adultXdir[t%52]                                  << "," << (int) ple[mm].adultYdir[t%52]               << endl;
        }
      }   
    }
    
    //Write output every 15 years (cycle of complete new population)        
    if(t % (A_MAX) == 5){ popStruct(mypopulation, ple,aliveple,t);}
 
  } //end of timeloop
     
  myfile.close() ; mypopulation.close();
  return 0 ;  
} 




