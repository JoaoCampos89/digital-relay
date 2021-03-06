/**
 *  Digital Relay for using in teaching digital relay protection using low cost
 *  devices
 *
 *
 * Author: João Campos
 */

#include <Arduino.h>
//#include <math.h>
//#include <iostream>     // std::cout
//#include <complex>      // std::complex
// computes phasors
void computePhasor(double *, double *, unsigned int);
double computeRMS(double * , unsigned int );
double computePowerRMS(double * , double * , unsigned int );
void setupTimer1Interrupt(unsigned int, int);
//void computeRMS(int);
// computes Fourier algorithm
//complex<double> Fourier(double *, int);
//complex<double> Fourier(double *, int, int);
// computes
double Fourier_Real(double *, unsigned int, unsigned int);
double Fourier_Imag(double *, unsigned int, unsigned int);
void Sampler();
// fundamental frequency
int f = 60;
// Samples per period
unsigned int N = 16;
// circularBuffer
unsigned int cb = 0;
bool sampler = false;
double V[3][16];
double I[3][16];
double P_RMS[3];
double V_RMS[3];
double I_RMS[3];
double V_PHASOR[3][2];
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
    P_RMS[i] = computePowerRMS(V[i], I[i], N);
  }
   computePhasor(V_PHASOR[0], V[0], N);


    sampler = false;
  //    Serial.println("RMS Value: " + String(V[0][0]));
  }

  // updates each 3 seconds
 if(timer-timerLast>1000){
    timerLast = timer;
    //double test = Fourier_Real(V[0],  N, 1);
    //sprintf(msg, "RMS Value: %f", 5.00);
    Serial.println("V RMS Value: " + String(V_RMS[0]));
    Serial.println("I RMS Value: " + String(I_RMS[0]));
    Serial.println("Z RMS Value: " + String(V_RMS[0]/I_RMS[0]));
    Serial.println("P RMS Value: " + String(P_RMS[0]));
    Serial.println("S RMS Value: " + String(V_RMS[0]*I_RMS[0]));
    Serial.println("Va Fourier Real Value: " + String(Fourier_Real(V[0],  N, 0)));
      Serial.println("Va Fourier Imag Value: " + String(Fourier_Imag(V[0],  N, 0)));
 }








  //Fourier();

//Protection();

}
// Sampling all analog values
void Sampler(){
  for (int i = 0; i < 1; i++) {
    V[i][N-1] = (analogRead(i)*5.0)/1023.0;
  }
  for (int i = 0; i < 1; i++) {
    I[i][N-1] = (analogRead(i+3)*5.0)/1023.0;
  }
  // translating elements
  for (int i= 0; i<1; i++){
    for (unsigned int j = 0; j < N-1; j++) {
          V[i][j] =  V[i][j+1];
          I[i][j] =  I[i][j+1];
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

   OCR1A = (16000000/256)/(N*f);            // compare match register (16MHz/256)/Desired Frequency
   TCCR1B |= (1 << WGM12);   // CTC mode
   TCCR1B |= (1 << CS12);    // 256 prescaler
   TIMSK1 |= (1 << OCIE1A);  // enable timer compare interrupt
   interrupts();             // enable all interrupts



}


// computer real Fourier coeficientes
double Fourier_Real(double * s, unsigned int N, unsigned int harmonic){
  float teta = 2*3.14/N;
  float twoOverN = 2/N;

  double Xre = 0;
  double coeficents[N];
  for (unsigned int i = 0; i < N-1; i++){
    coeficents[i] = twoOverN*cos( teta*harmonic*i);
  }


  for (unsigned int i = 0; i < N-1; i++) {
    Xre = Xre + s[i]*coeficents[i];

  }
  return Xre;
}
// computer imag Fourier coeficientes
double Fourier_Imag(double * s , unsigned int N, unsigned int harmonic){
  double teta = 2*3.14/N;
  double twoOverN = 2/N;
  double Xim = 0;
  double coeficents[N];
  for (unsigned int i = 0; i < N-1; i++){
    coeficents[i] = twoOverN*sin( teta*harmonic*i);
  }

  for (unsigned int i = 0; i < N-1; i++) {
    Xim = Xim - s[i]*coeficents[i];
  }
  return Xim;
}





void computePhasor(double * s , double * v, unsigned int N ){

    s[0] = Fourier_Real(v, N, 0);
    s[1] = Fourier_Imag(v, N, 0);

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
