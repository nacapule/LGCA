#include "WELL1024a.h"
#include <iostream>
#include <time.h>

#define BINS 100
#define MIT 1000000000
int main(){
  
unsigned int init[32];
unsigned int b[BINS];
double val;

for(int i=0;i<BINS;i++){
 b[i]=0;
}

time_t t0=time(NULL);

InitWELLRNG1024a(&init[0]);
for(int i=0;i<MIT;i++){
// std::cout << WELLRNG1024a() << "\n";
 val=WELLRNG1024a();
 /*for(int j=0;j<BINS;j++){
  if(val<(j+1)*1.0/BINS) {
   b[j]++; 
   break;
  }
 }*/
} 
time_t t1=time(NULL);

 for(int i=0;i<BINS;i++){
 std::cout << i << " "<< b[i]/(double)MIT << "\n";
 }
 std::cout << t1-t0 << "\n"; 
 
 



}