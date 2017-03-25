#ifndef Phasor_h
#define Phasorh

//#include <Arduino.h>


class Phasor
{
private:

  float _real;
  float _imag;


public:
	// Constructor
  Phasor(float, float);
  float getReal();
  float getImag();

};




#endif
