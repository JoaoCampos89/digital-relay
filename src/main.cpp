/**
 *  Digital Relay for using in teaching digital relay protection using low cost
 *  devices
 *
 *
 * Author: João Campos
 */

#include <Arduino.h>
#include <math.h>
//#include <iostream>     // std::cout
//#include <complex>      // std::complex
// computes phasors
void computePhasors(int);

void computeRMS(int);
// computes Fourier algorithm
//complex<double> Fourier(double *, int);
//complex<double> Fourier(double *, int, int);
// computes
double Fourier_Real(double *, int);
double Fourier_Imag(double *, int);
void Sampler();
// fundamental frequency
int f = 60;
// Samples per period
unsigned int N = 4;
// circularBuffer
unsigned int cb = 0;
bool sampler = false;
double * V[3];
double * I[3];
double * V_PHASOR[3];
double * I_PHASOR[3];

using namespace std;


void setup(){

  // initialize timer1
  // http://www.robotshop.com/letsmakerobots/arduino-101-timers-and-interrupts
   noInterrupts();           // disable all interrupts
   TCCR1A = 0;
   TCCR1B = 0;
   TCNT1  = 0;

   OCR1A = (16000000/256)/(N*f);            // compare match register (16MHz/256)/Desired Frequency
   TCCR1B |= (1 << WGM12);   // CTC mode
   TCCR1B |= (1 << CS12);    // 256 prescaler
   TIMSK1 |= (1 << OCIE1A);  // enable timer compare interrupt
   interrupts();             // enable all interrupts

}

ISR(TIMER1_COMPA_vect)          // timer compare interrupt service routine
{
  sampler = true;
}




void loop(){
  if(sampler){

    Sampler();

    computePhasors(N);


    sampler = false;
  }

  //Fourier();

//Protection();

};
// Sampling all analog values
void Sampler(){
  for (int i = 0; i < 3; i++) {
    V[i][N-1] = analogRead(i);
  }
  for (int i = 0; i < 3; i++) {
    I[i][N-1] = analogRead(i+3);
  }
  // translating elements
  for (int i= 0; i<3; i++){
    for (unsigned int j = 0; j < N-1; j++) {
          V[i][j] =  V[i][j+1];
          I[i][j] =  I[i][j+1];
    }

  }
};

// Sampling all analog values
// N must be multiple of 2
// for cb equal
void circularSampler(unsigned int &cb, unsigned int N){
  for (int i = 0; i < 3; i++) {
    V[i][cb] = analogRead(i);
  }
  for (int i = 0; i < 3; i++) {
    I[i][cb] = analogRead(i+3);
  }
  cb++;
  cb = cb & N;
};




// computer real Fourier coeficientes
double Fourier_Real(double * s, int N){
  double teta = 2*PI/N;
  double twoOverN = 2/N;

  double Xre = 0;


  for (int i = 0; i < N-1; i++) {
    Xre = Xre + twoOverN*s[i]*cos(teta*i);

  }
  return Xre;
};
// computer imag Fourier coeficientes
double Fourier_Imag(double * s, int N){
  double teta = 2*PI/N;
  double twoOverN = 2/N;
  double Xim = 0;

  for (int i = 0; i < N-1; i++) {
    Xim = Xim - twoOverN*s[i]*sin(teta*i);
  }
  return Xim;
};

void computePhasors(int N){
  for (int i = 0; i < 3; i++) {
    V_PHASOR[i][0] = Fourier_Real(V[i], N);
    V_PHASOR[i][1] = Fourier_Imag(V[i], N);
    I_PHASOR[i][0] = Fourier_Real(I[i], N);
    I_PHASOR[i][1] = Fourier_Imag(I[i], N);
  }

}


double computeRMS(double * s, unsigned int totalSamples){
  double rms = 0;
  for (unsigned int i = 0; i < totalSamples; i++) {
      rms = rms + s[i]*s[i];
  }
  rms = sqrt(rms/totalSamples);
  return rms;

}
/**
 * [computePowerRMS computes RMS power]
 * @param  c            [current array]
 * @param  v            [voltage array]
 * @param  totalSamples [description]
 * @return              [rms]
 */
double computePowerRMS(double * c, double * v, unsigned int totalSamples){
  double rms = 0;
  for (unsigned int i = 0; i < totalSamples; i++) {
      rms = rms + c[i]*v[i]*c[i]*v[i];
  }
  rms = sqrt(rms/totalSamples);
  return rms;

}








/**
 * without circular buffer
 */
/*complex<double> Fourier(double * s, int N){
  double teta = 2*pi()/N;

  double Xre = 0;
  double Xim = 0;

  for (int i = 0; i < N; i++) {
    Xre = Xre + (2/N)*s(m)*cos(teta*m);
    Xim = Xim - (2/N)*s(m)*sin(teta*m);
  }

return complex<double> phasor(Xre,Xim);
}

complex<double> Fourier(double * s, int cb, int N){
  double teta = 2*pi()/N;

  double Xre = 0;
  double Xim = 0;

  for (int i = 0; i < N; i++) {
    Xre = Xre + (2/N)*s(m)*cos(teta*m);
    Xim = Xim - (2/N)*s(m)*sin(teta*m);
  }

return complex<double> phasor(Xre,Xim);
}*/
