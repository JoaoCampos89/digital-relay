#ifndef Estimator_h
#define Estimator_h

#include <complex>
#include <iostream>
//#include <Arduino.h>
using namespace std;

class Fourier
private:

  float * value;
  complex _previousPhasor;
  complex phasor;
  int N;

public:
	// Constructor
	 Estimator(int);
   virtual  complex computePhasor();




};




#endif
