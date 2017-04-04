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
void computePhasor(double *[], unsigned int);
double computeRMS(double * , unsigned int );
void setupTimer1Interrupt(unsigned int, int);
//void computeRMS(int);
// computes Fourier algorithm
//complex<double> Fourier(double *, int);
//complex<double> Fourier(double *, int, int);
// computes
double Fourier_Real(double *, unsigned int);
//double Fourier_Imag(double *, unsigned int);
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
double V_RMS[3];
double I_RMS[3];
double * V_PHASOR[3];
double * I_PHASOR[3];
long timerLast = 0;
char msg[100];
using namespace std;


void setup(){
  for (int i = 0; i < 6; i++) {
    pinMode(i, INPUT);
    /* code */
  }
  setupTimer1Interrupt(N, f);
  Serial.begin(115200);

}

ISR(TIMER1_COMPA_vect)          // timer compare interrupt service routine
{
  sampler = true;
}




void loop(){
  long timer = millis();
  if(sampler){

    Sampler();
  //  V_PHASOR[0][0] =  computeRMS(V[0], N);
  //  V_PHASOR[0][0] = Fourier_Real(V[0], N);
    for (unsigned int i = 0; i < 1; i++) {
    V_RMS[i] = computeRMS(V[i], N);
    I_RMS[i] = computeRMS(I[i], N);
  }
  // computePhasor(V_PHASOR,N);


    sampler = false;
  //    Serial.println("RMS Value: " + String(V[0][0]));
  }

  // updates each 3 seconds
 if(timer-timerLast>1000){
    timerLast = timer;
    //sprintf(msg, "RMS Value: %f", 5.00);
    Serial.println("RMS Value: " + String(V_RMS[0]));
 }








  //Fourier();

//Protection();

}
// Sampling all analog values
void Sampler(){
  for (int i = 0; i < 1; i++) {
    V[i][N-1] = analogRead(i);
  }
/*  for (int i = 0; i < 1; i++) {
    I[i][N-1] = analogRead(i+3);
  }*/
  // translating elements
  for (int i= 0; i<1; i++){
    for (unsigned int j = 0; j < N-1; j++) {
          V[i][j] =  V[i][j+1];
          //I[i][j] =  I[i][j+1];
    }

  }
}

// Sampling all analog values
// N must be multiple of 2
// for cb equal
/*void circularSampler(unsigned int &cb, unsigned int N){
  for (int i = 0; i < 3; i++) {
    V[i][cb] = analogRead(i);
  }
  for (int i = 0; i < 3; i++) {
    I[i][cb] = analogRead(i+3);
  }
  cb++;
  cb = cb & N;
}*/

void setupTimer1Interrupt(unsigned int N, int f){
  // initialize timer1
  // http://www.robotshop.com/letsmakerobots/arduino-101-timers-and-interrupts
   noInterrupts();           // disable all interrupts
   TCCR1A = 0;
   TCCR1B = 0;
   TCNT1  = 0;

   OCR1A = (16000000/256)/(4*60);            // compare match register (16MHz/256)/Desired Frequency
   TCCR1B |= (1 << WGM12);   // CTC mode
   TCCR1B |= (1 << CS12);    // 256 prescaler
   TIMSK1 |= (1 << OCIE1A);  // enable timer compare interrupt
   interrupts();             // enable all interrupts



}


// computer real Fourier coeficientes
double Fourier_Real(double * s, unsigned int N){
  double teta = 2*PI/N;
  double twoOverN = 2/N;

  double Xre = 0;


  for (unsigned int i = 0; i < N-1; i++) {
    Xre = Xre + twoOverN*s[i]*cos(teta*i);

  }
  return Xre;
}
// computer imag Fourier coeficientes
double Fourier_Imag(double * s , unsigned int N){
  double teta = 2*PI/N;
  double twoOverN = 2/N;
  double Xim = 0;

  for (unsigned int i = 0; i < N-1; i++) {
    Xim = Xim - twoOverN*s[i]*sin(teta*i);
  }
  return Xim;
}





void computePhasor(double * s , double * v, unsigned int N){

    s[0] = Fourier_Real(v, N);
    s[1] = Fourier_Imag(v, N);

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
/*double computePowerRMS(double * c, double * v, unsigned int totalSamples){
  double rms = 0;
  for (unsigned int i = 0; i < totalSamples; i++) {
      rms = rms + c[i]*v[i]*c[i]*v[i];
  }
  rms = sqrt(rms/totalSamples);
  return rms;

}*/








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
