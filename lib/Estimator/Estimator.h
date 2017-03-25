#ifndef Estimator_h
#define Estimator_h

#include <complex>
#include <iostream>
//#include <Arduino.h>
using namespace std;

class Estimator
{
private:

  float * value;
  complex _previousPhasor;
  complex phasor;
  int N;

public:
	// Constructor
	 Estimator(int);
   Algorithm algorithm(int);
   virtual  complex computePhasor();
};




#endif
