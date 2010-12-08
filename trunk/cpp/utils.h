//----------------------------------------------------------------------------//
// Generic functions
//----------------------------------------------------------------------------//

double rnorm            (double mu, double sd)                 ;      // function random 
double max              (double first, double second)          ;      // function max 
double min              (double first, double second)          ;      // function min 


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
