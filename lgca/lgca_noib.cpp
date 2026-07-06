#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <list>
#include <cstdlib>
#include <cstring>  /*libraries*/
#include <vector>
#include <bitset>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdlib.h>

#ifdef CIMG
#include "../CImg/CImg.h"  /*library for the graphics output*/
#endif

// #define HISTOG

#define sqr(a) ((a)*(a)) /*macro for the square (more efficient)*/

#include "../rng/WELL1024a.h"    /*random number generator*/

#define NODES 6             /*number of velocity channels, for now only hexagonal*/

#define xdim 120            /*lattice size*/
#define ydim 120
// #define xdim 1
// #define ydim 1

#define D 1                 /* diffusion coefficient for scaling*/
#define v 16                 /* particle velocity for scaling*/

// LGCA object
template <class T, int DX, int DY> class lgca{
    
    public:
        
        std::vector< std::vector< T > > lattice;        /*pre interaction configurations*/
        std::vector< std::vector< T > > lattice_temp;   /*post interaction configurations*/
        std::vector< double> probi;                     /*channel probabilities (boson)*/
	std::vector< std::vector <double> > ci;        /*veclocity channel vectors*/
        
        lgca(int dimx, int dimy){
            
            probi.resize(NODES);
            
            ci.resize(6);
            for(int i=0; i<6; i++){
                ci[i].resize(2);
            }

            ci[0][0]=1;
            ci[0][1]=0;

            ci[1][0]=0.5;
            ci[1][1]=sqrt(3.0)/2.0;

            ci[2][0]=-0.5;
            ci[2][1]=sqrt(3.0)/2.0;
                                            /*defition of velocity channels*/
            ci[3][0]=-1;
            ci[3][1]=0;

            ci[4][0]=-0.5;
            ci[4][1]=sqrt(3.0)/(-2.0);

            ci[5][0]=0.5;
            ci[5][1]=sqrt(3.0)/(-2.0);

            
            lattice.resize(dimx);
            lattice_temp.resize(dimx);
            for(int i=0; i<dimx; i++){
                lattice[i].resize(dimy);
                lattice_temp[i].resize(dimy);
            }
        }
        T get(int x, int y){
            return lattice[x][y];
        }                                   /*helping functions*/
        
        T get_temp(int x, int y){
            return lattice_temp[x][y];
        }
        
//         TODO void place_parts(
        
        
//creation of a new random post-interaction configuration FOR NUMBER CONSERVING INTERACTIONS
        void shuffle(std::bitset<NODES>& newconf,const std::bitset<NODES>& oldconf){
            int chosen;
            std::bitset<NODES>counted;
            newconf.reset();
            counted.reset();
            for(int i=0; i<NODES; i++){
                chosen=int(round(WELLRNG1024a()*(NODES-1)));
//         std::cout << "yeschosen " << chosen << std::flush;
                while(counted[chosen]){
                    chosen=int(round(WELLRNG1024a()*(NODES-1)));
//             std::cout << "nochosen " << chosen << std::flush;
                }
                newconf.set(chosen,oldconf[i]);
                counted.set(chosen);
            }
        }

//creation of a new random post-interaction configuration FOR NUMBER CONSERVING INTERACTIONS for BOSON models (unused)
        void shuffle(std::vector<int>& newconf, const std::vector<int>& oldconf){
            int chosen;
            int nuparts=0;
            int parts;
            newconf.assign(NODES,0);
            for(int i=0; i<NODES; i++){
                nuparts+=oldconf[i];
            }
            while(nuparts>0){
                chosen=int(round(WELLRNG1024a()*(NODES-1)));
                parts=int(round(WELLRNG1024a()*nuparts));
                newconf[chosen]+=parts;
                nuparts-=parts;
            }
        }
        
        
        //calculation of the probability of a given post interaction configuration (FERMIONS)
        double interaction(const T& conf, double sens, int X, int Y){
            
            double prob=0, parts=0;
            
            for(int i=0; i<NODES; i++){
                for(int j=0; j<NODES; j++){
//                     prob+=conf[i]*lattice[X][Y][j]*sqr((ci[i][0]*ci[j][0]+ci[i][1]*ci[j][1]))*sens;
//                     prob+=conf[i]*lattice[(X+1)%DX][Y][j]*sqr((ci[i][0]*ci[j][0]+ci[i][1]*ci[j][1]))*sens;
//                     prob+=conf[i]*lattice[X][(Y+1)%DX][j]*sqr((ci[i][0]*ci[j][0]+ci[i][1]*ci[j][1]))*sens;
//                     prob+=conf[i]*lattice[(X-1+DX)%DX][Y][j]*sqr((ci[i][0]*ci[j][0]+ci[i][1]*ci[j][1]))*sens;
//                     prob+=conf[i]*lattice[X][(Y-1+DX)%DX][j]*sqr((ci[i][0]*ci[j][0]+ci[i][1]*ci[j][1]))*sens;
//                     prob+=conf[i]*lattice[(X+1)%DX][(Y-1+DX)%DX][j]*sqr((ci[i][0]*ci[j][0]+ci[i][1]*ci[j][1]))*sens;
//                     prob+=conf[i]*lattice[(X-1+DX)%DX][(Y+1)%DX][j]*sqr((ci[i][0]*ci[j][0]+ci[i][1]*ci[j][1]))*sens;
                    
//                     prob+=conf[i]*lattice[(X+1)%DX][Y][j]*(0.5*(ci[i][0]*ci[j][0]+ci[i][1]*ci[j][1])-0.5*sqr(ci[i][0]*ci[0][0]+ci[i][1]*ci[0][1])+sqr(cos(double(i+j)*M_PI/6.0)*ci[0][0]+sin(double(i+j)*M_PI/6.0)*ci[0][1])+0.5*sqr(ci[0][0]))*(fmod(-1.0*double(i+j)*M_PI/6.0+2*M_PI,2*M_PI)<M_PI?sens:-1.0*sens);
//                     prob+=conf[i]*lattice[(X+1)%DX][(Y-1+DX)%DX][j]*(0.5*(ci[i][0]*ci[j][0]+ci[i][1]*ci[j][1])-0.5*sqr(ci[i][0]*ci[1][0]+ci[i][1]*ci[1][1])+sqr(cos(double(i+j)*M_PI/6.0)*ci[1][0]+sin(double(i+j)*M_PI/6.0)*ci[1][1])+0.5*sqr(ci[1][0]))*(fmod(M_PI/3.0-double(i+j)*M_PI/6.0+2*M_PI,2*M_PI)<M_PI?sens:-1.0*sens);
//                     prob+=conf[i]*lattice[X][(Y-1+DX)%DX][j]*(0.5*(ci[i][0]*ci[j][0]+ci[i][1]*ci[j][1])-0.5*sqr(ci[i][0]*ci[2][0]+ci[i][1]*ci[2][1])+sqr(cos(double(i+j)*M_PI/6.0)*ci[2][0]+sin(double(i+j)*M_PI/6.0)*ci[2][1])+0.5*sqr(ci[2][0]))*(fmod(2.0*M_PI/3.0-double(i+j)*M_PI/6.0+2*M_PI,2*M_PI)<M_PI?sens:-1.0*sens);
//                     prob+=conf[i]*lattice[(X-1+DX)%DX][Y][j]*(0.5*(ci[i][0]*ci[j][0]+ci[i][1]*ci[j][1])-0.5*sqr(ci[i][0]*ci[3][0]+ci[i][1]*ci[3][1])+sqr(cos(double(i+j)*M_PI/6.0)*ci[3][0]+sin(double(i+j)*M_PI/6.0)*ci[3][1])+0.5*sqr(ci[3][0]))*(fmod(M_PI-double(i+j)*M_PI/6.0+2*M_PI,2*M_PI)<M_PI?sens:-1.0*sens);
//                     prob+=conf[i]*lattice[(X-1+DX)%DX][(Y+1)%DX][j]*(0.5*(ci[i][0]*ci[j][0]+ci[i][1]*ci[j][1])-0.5*sqr(ci[i][0]*ci[4][0]+ci[i][1]*ci[4][1])+sqr(cos(double(i+j)*M_PI/6.0)*ci[4][0]+sin(double(i+j)*M_PI/6.0)*ci[4][1])+0.5*sqr(ci[4][0]))*(fmod(4.0*M_PI/3.0-double(i+j)*M_PI/6.0+2*M_PI,2*M_PI)<M_PI?sens:-1.0*sens);
//                     prob+=conf[i]*lattice[X][(Y+1)%DX][j]*(0.5*(ci[i][0]*ci[j][0]+ci[i][1]*ci[j][1])-0.5*sqr(ci[i][0]*ci[5][0]+ci[i][1]*ci[5][1])+sqr(cos(double(i+j)*M_PI/6.0)*ci[5][0]+sin(double(i+j)*M_PI/6.0)*ci[5][1])+0.5*sqr(ci[5][0]))*(fmod(5.0*M_PI/3.0-double(i+j)*M_PI/6.0+2*M_PI,2*M_PI)<M_PI?sens:-1.0*sens);
                    
//                     prob+=conf[i]*lattice[X][Y][j]*(ci[i][0]*ci[j][0]+ci[i][1]*ci[j][1])*sens;
                    prob+=conf[i]*lattice[(X+1)%DX][Y][j]*(ci[i][0]*ci[j][0]+ci[i][1]*ci[j][1])*sens;
                    prob+=conf[i]*lattice[X][(Y+1)%DX][j]*(ci[i][0]*ci[j][0]+ci[i][1]*ci[j][1])*sens;
                    prob+=conf[i]*lattice[(X-1+DX)%DX][Y][j]*(ci[i][0]*ci[j][0]+ci[i][1]*ci[j][1])*sens;
                    prob+=conf[i]*lattice[X][(Y-1+DX)%DX][j]*(ci[i][0]*ci[j][0]+ci[i][1]*ci[j][1])*sens;
                    prob+=conf[i]*lattice[(X+1)%DX][(Y-1+DX)%DX][j]*(ci[i][0]*ci[j][0]+ci[i][1]*ci[j][1])*sens;
                    prob+=conf[i]*lattice[(X-1+DX)%DX][(Y+1)%DX][j]*(ci[i][0]*ci[j][0]+ci[i][1]*ci[j][1])*sens;
		    
		    parts+=lattice[(X+1)%DX][Y][j]+lattice[X][(Y+1)%DX][j]+lattice[(X-1+DX)%DX][Y][j];
		    parts+=lattice[X][(Y-1+DX)%DX][j]+lattice[(X+1)%DX][(Y-1+DX)%DX][j]+lattice[(X-1+DX)%DX][(Y+1)%DX][j];
                }
            }
//             if(parts!=0) prob/=parts;
            return prob;
            
        }
        
        
        //calculation of the channel probabilities for BOSON models
        void pdfchannels(int X, int Y, double sens){
            
            int numparts=0,parts,partneigh;
            
            double s1,s2,mx,z=0;
            
            probi.assign(NODES,0);
            
//             for(int i=0;i<NODES;i++){
//                 
//                 parts=lattice[X][Y][i]/*+lattice[(X+1)%DX][Y][i]+lattice[X][(Y+1)%DX][i]+lattice[(X-1+DX)%DX][Y][i]+lattice[X][(Y-1+DX)%DX][i]+lattice[(X+1)%DX][(Y-1+DX)%DX][i]+lattice[(X-1+DX)%DX][(Y+1)%DX][i]*/;
// 		partneigh=lattice[(X+1)%DX][Y][i]+lattice[X][(Y+1)%DX][i]+lattice[(X-1+DX)%DX][Y][i]+lattice[X][(Y-1+DX)%DX][i]+lattice[(X+1)%DX][(Y-1+DX)%DX][i]+lattice[(X-1+DX)%DX][(Y+1)%DX][i];
//                 numparts+=parts;
// //                 for(int j=0;j<NODES;j++){
// // //                     probi[j]-=2*parts*(ci[j][0]*ci[i][0]+ci[j][1]*ci[i][1])*sens*sqr(v/4*D);
// // // 		    probi[j]+=partneigh*(ci[j][0]*ci[i][0]+ci[j][1]*ci[i][1])*sens*sqr(v/4*D);
// //                     probi[j]-=2*parts*sqr(ci[j][0]*ci[i][0]+ci[j][1]*ci[i][1])*sens*sqr(v/4*D);
// // 		    probi[j]+=partneigh*sqr(ci[j][0]*ci[i][0]+ci[j][1]*ci[i][1])*sens*sqr(v/4*D);
// //                 }
//                 
//             }
// 	    std::cout << "preint" << "\n";
	    for(int i=0;i<NODES;i++){
                
                parts=/*lattice[X][Y][i]+*/lattice[(X+1)%DX][Y][i]+lattice[X][(Y+1)%DX][i]+lattice[(X-1+DX)%DX][Y][i]+lattice[X][(Y-1+DX)%DX][i]+lattice[(X+1)%DX][(Y-1+DX)%DX][i]+lattice[(X-1+DX)%DX][(Y+1)%DX][i];
// 		std::cout << "c[" << i <<"]="<<parts ; 
                numparts+=parts;
                for(int j=0;j<NODES;j++){
                    probi[j]+=parts*(ci[j][0]*ci[i][0]+ci[j][1]*ci[i][1])*sens;
//                     probi[j]+=parts*sqr(ci[j][0]*ci[i][0]+ci[j][1]*ci[i][1])*sens;
                    
//                     probi[j]+=lattice[(X+1)%DX][Y][i]*(/*0.5*(ci[j][0]*ci[i][0]+ci[j][1]*ci[i][1])*/-0.5*sqr((ci[j][0]*ci[0][0]+ci[j][1]*ci[0][1]))/*+sqr(cos(double(i+j)*M_PI/6.0)*ci[0][0]+sin(double(i+j)*M_PI/6.0)*ci[0][1])+0.5*sqr(ci[0][0])*/)*sens;
//                     probi[j]+=lattice[(X+1)%DX][(Y-1+DX)%DX][i]*(/*0.5*(ci[j][0]*ci[i][0]+ci[j][1]*ci[i][1])*/-0.5*sqr((ci[j][0]*ci[1][0]+ci[j][1]*ci[1][1]))/*+sqr(cos(double(i+j)*M_PI/6.0)*ci[1][0]+sin(double(i+j)*M_PI/6.0)*ci[1][1])+0.5*sqr(ci[1][0])*/)*sens;
//                     probi[j]+=lattice[X][(Y-1+DX)%DX][i]*(/*0.5*(ci[j][0]*ci[i][0]+ci[j][1]*ci[i][1])*/-0.5*sqr((ci[j][0]*ci[2][0]+ci[j][1]*ci[2][1]))/*+sqr(cos(double(i+j)*M_PI/6.0)*ci[2][0]+sin(double(i+j)*M_PI/6.0)*ci[2][1])+0.5*sqr(ci[2][0])*/)*sens;
//                     probi[j]+=lattice[(X-1+DX)%DX][Y][i]*(/*0.5*(ci[j][0]*ci[i][0]+ci[j][1]*ci[i][1])*/-0.5*sqr((ci[j][0]*ci[3][0]+ci[j][1]*ci[3][1]))/*+sqr(cos(double(i+j)*M_PI/6.0)*ci[3][0]+sin(double(i+j)*M_PI/6.0)*ci[3][1])+0.5*sqr(ci[3][0])*/)*sens;
//                     probi[j]+=lattice[(X-1+DX)%DX][(Y+1)%DX][i]*(/*0.5*(ci[j][0]*ci[i][0]+ci[j][1]*ci[i][1])*/-0.5*sqr((ci[j][0]*ci[4][0]+ci[j][1]*ci[4][1]))/*+sqr(cos(double(i+j)*M_PI/6.0)*ci[4][0]+sin(double(i+j)*M_PI/6.0)*ci[4][1])+0.5*sqr(ci[4][0])*/)*sens;
//                     probi[j]+=lattice[X][(Y+1)%DX][i]*(/*0.5*(ci[j][0]*ci[i][0]+ci[j][1]*ci[i][1])*/-0.5*sqr((ci[j][0]*ci[5][0]+ci[j][1]*ci[5][1]))/*+sqr(cos(double(i+j)*M_PI/6.0)*ci[5][0]+sin(double(i+j)*M_PI/6.0)*ci[5][1])+0.5*sqr(ci[5][0])*/)*sens;
                }
                
            }
// 	    std::cout << "\n";
            
            //subtract maximum element so that the exponential does not explode, does not affect probabilties
	    mx=*std::max_element( probi.begin(), probi.end() );
	    for(int i=0; i<NODES; i++){
	      probi[i]-=mx;
	    }
            
            //calculation of the partition function
            for(int i=0; i<NODES;i++){
                if(numparts!=0) probi[i]/=numparts;
                probi[i]=exp(probi[i]);
                z+=probi[i];
            }
            //normalization
            for(int i=0; i<NODES; i++){
                probi[i]/=z;
// 		std::cout << "prob(c[" << i<< "])=" << probi[i]; 
            }
            //calculation of the CDF
            for(int i=1; i<NODES; i++){
                probi[i]+=probi[i-1];
            }
            
            
        }
        
        // random channel generator accordgin to the CDF calculated in previous function
        int randchangen(){
            
            
            double r=WELLRNG1024a();
            
            for(int i=0; i<NODES; i++){
                
                if(r <= probi[i]) return i;
                
            }
            
        }
        
        // particle channel reassignment (BOSON)
        void newconfgen(std::vector<int>& temp,int X,int Y,double sens){
            
//             double P,Pnew,alpha;
//             std::vector<int>temp2;
//     
//             temp=lattice[X][Y];
//                     
//             P=exp(interaction(temp,sens,X,Y));
//                     
// //                     std::cout << "temp " << "c1= " << temp[0] << " " << "c2= " << temp[1] << " " << "c3= " << temp[2] << " " << "c4= " << temp[3] << " " << "c5= " << temp[4] << " " << "c6= " << temp[5] << " " << "P= " << P << "\n";
//                     
//             for(int k=0; k<12; k++){
//                         
// //                         std::cout << "t: " << t <<"i: " << i << "j: " << j << "\n" << std::flush;
//                         
//                 shuffle(temp2,temp);
// //                         std::cout << "shuffle" << "\n" << std::flush;
//                 Pnew=exp(interaction(temp2,sens,X,Y));
//                         
// //                         std::cout << "temp " << "c1= " << temp2[0] << " " << "c2= " << temp2[1] << " " << "c3= " << temp2[2] << " " << "c4= " << temp2[3] << " " << "c5= " << temp2[4] << " " << "c6= " << temp2[5] << " " << "P= " << Pnew << "\n";
//                         
//                 alpha=Pnew/P;
//                 if(WELLRNG1024a()<alpha){
//                     temp=temp2;
//                     P=Pnew;
//                 }
//                         
//             }
            
            temp.assign(NODES,0);
            int partno=0;
            
            for(int i=0;i<NODES;i++){
                partno+=lattice[X][Y][i];
            }
            
            pdfchannels(X,Y,sens);
            
            for(int i=0; i<partno; i++){
                temp[randchangen()]+=1;
            }
            
//             std::cout<< "postint"<<"\n";
	    for(int i=0; i<NODES; i++){
// 	      std::cout << "c[" << i <<"]="<<temp[i] ; 
	    }
//             std::cout << "\n" << std::flush;
        }

        //stochastic post interaction configuration selection (FERMION)
        void newconfgen(std::bitset<NODES>& temp,int X,int Y,double sens){
            
            double P,Pnew,alpha;
    
            std::bitset<NODES>temp2;
            
            temp=lattice[X][Y];
    
                    
//             P=exp(interaction(temp,sens,X,Y));
	    P=interaction(temp,sens,X,Y);
                    
//                     std::cout << "temp " << "c1= " << temp[0] << " " << "c2= " << temp[1] << " " << "c3= " << temp[2] << " " << "c4= " << temp[3] << " " << "c5= " << temp[4] << " " << "c6= " << temp[5] << " " << "P= " << P << "\n";
                    
            //Metropolis algorithm
            for(int k=0; k<12; k++){
                        
//                         std::cout << "t: " << t <<"i: " << i << "j: " << j << "\n" << std::flush;
                        
                shuffle(temp2,temp);
//                         std::cout << "shuffle" << "\n" << std::flush;
//                 Pnew=exp(interaction(temp2,sens,X,Y));
                Pnew=interaction(temp2,sens,X,Y);
//                         std::cout << "temp " << "c1= " << temp2[0] << " " << "c2= " << temp2[1] << " " << "c3= " << temp2[2] << " " << "c4= " << temp2[3] << " " << "c5= " << temp2[4] << " " << "c6= " << temp2[5] << " " << "P= " << Pnew << "\n";
                        
                alpha=Pnew-P;
                if(log(WELLRNG1024a())<alpha){
                    temp=temp2;
                    P=Pnew;
                }
                        
            }
    
        }
        
};


//test if nodes are empty or not
bool empty(const std::vector<int>& site){
    for(int i=0; i<site.size(); i++){
        if(site[i]>0) return false;
    }
    return true;
}

bool empty(const std::bitset<NODES>& site){
    return site.none();
}

//initialization of the lattice (FERMION)
void initialize(std::bitset<NODES>& conf, int val=0){
    conf.reset();
}

//initialization of the lattice (BOSON)
void initialize(std::vector<int>& conf,int val=0){
    conf.assign(NODES,val);
//     conf.assign(NODES,0);
//     conf[0]=val;
//     conf[3]=val;
}

//fucntions for adding single particles to channels
void add_part(std::vector<int>& conf, int chann){
    conf[chann]+=1;
}

void add_part(std::bitset<NODES>& conf, int chann){
    conf.set(chann);
}

//main program
int main(int argc, char** args){
// std::cout << "inicia" << "\n" << std::flush;
// int xdim=atoi(args[1]);
// int ydim=atoi(args[2]);
int tsteps=atoi(args[1]);       //duration to be input in command line at execution
int iters=atoi(args[2]);        //realizations to be input in command line at execution
double sens=atof(args[3]);      //sensitivity to be input in command line at execution
double dens=atof(args[4]);      //initial density to be input in command line at execution
// double syncp=atof(args[5]);  //asynchonicity to be input in command line at execution

// double heating=2.0/(1.0*tsteps); //for gradual heating
// double senseff;


unsigned int init[32];
for(int i=0;i<32;i++){init[i]=time(NULL);}  //rng initialization
InitWELLRNG1024a(&init[0]);

// using states = std::bitset<NODES>;   //uncomment for FERMIONS
using states = std::vector<int>;        //uncomment for BOSONS

states temp;

//new lattice
lgca<states,xdim,ydim> *s=new lgca<states,xdim,ydim>(xdim,ydim);

// std::cout << "declaracion" << "\n" << std::flush;

//initialization of graphic output
#ifdef CIMG   
// cimg_library::CImg<unsigned char> p(4*(xdim+1)+(ydim-1)*2,4*ydim+4,1,3,0);
cimg_library::CImg<unsigned char> p(4*xdim+4,4*ydim+4,1,3,0);

cimg_library::CImgDisplay dis(p,"result");
   
  p.fill(0);
   
  int col=0;
 
  unsigned char white[] = {255,255,255};
  
#endif
  
  //initialization of variables to store observable outputs
long int nuparts;

std::vector <double> polaro,sqpolaro,nemato,sqnemato,entro,sqentro,band,sqband;
polaro.assign(tsteps,0);
sqpolaro.assign(tsteps,0);
nemato.assign(tsteps,0);
sqnemato.assign(tsteps,0);
entro.assign(tsteps,0);
sqentro.assign(tsteps,0);
band.assign(tsteps,0);
sqband.assign(tsteps,0);

std::vector <double> polmax,nemmax,entromax,bandmax;
polmax.assign(tsteps,0);
nemmax.assign(tsteps,0);
entromax.assign(tsteps,0);
bandmax.assign(tsteps,0);

//loop for every realization
for(int r=0; r<iters; r++){
  
  nuparts=int(dens)*NODES*xdim*ydim;    //number of particles in the lattice (for density>1)
//     nuparts=int(dens)*xdim*ydim;
//   senseff=0;
  
    //lattice initialization
    for(int i=0; i<xdim; i++){
        for(int j=0; j<ydim; j++){
            if(dens>=1){
                initialize(s->lattice[i][j],int(dens));
                initialize(s->lattice_temp[i][j]);
            }else{
                initialize(s->lattice[i][j]);
                initialize(s->lattice_temp[i][j]);
            }
            
            
        }
    
    }
    
    //assignment of random individual particles according to density value
    for(int i=0; i<xdim; i++){
        for(int j=0; j<ydim; j++){
            for(int k=0; k<NODES; k++){
                if(WELLRNG1024a()<( dens-int(dens) ) ) {
		  add_part(s->lattice[i][j],k);
//                     add_part(s->lattice[i][j],1);
		  nuparts+=1;
		}
            }
            
//             std::cout << "i= " << i << " " << "j= " << j << " " << "c1= " << s->lattice[i][j][0] << " " << "c2= " << s->lattice[i][j][1] << " " << "c3= " << s->lattice[i][j][2] << " " << "c4= " << s->lattice[i][j][3] << " " << "c5= " << s->lattice[i][j][4] << " " << "c6= " << s->lattice[i][j][5] << " ";
        }
//         std::cout << "\n";
    }
    
    
    
    
    //loop for every time step
//     std::cout << "inicializacion" << "\n" << std::flush;
    for(int t=0; t<tsteps; t++){
      
//       int rparts=0;
//       for(int i=0; i<xdim; i++){
//       for(int j=0; j<ydim; j++){
// 	for(int k=0; k<NODES; k++){
// 	  rparts+=s->lattice[i][j][k];
// 	}
//       }
//     }
//     
//     std::cout << rparts << " " << double(rparts)/double(6*xdim*ydim) << "\n" << std::flush;
//         senseff=t>(tsteps/2)?senseff-heating:(heating*t); //linear heating
// 	std::cout << t << " " << senseff << "\n" << std::flush;
//         std::cout << "interaccion" << "\n";
        
        //loop for every node
        for(int i=0; i<xdim; i++){
            for(int j=0; j<ydim; j++){
                
                //if node not empty, select post interaction configuration
                if(!empty(s->lattice[i][j])){
//                     std::cout << t << "\n";
                    s->newconfgen(temp,i,j,/*senseff*/sens); //sense for no heating, senseff for heating
                    
                    s->lattice_temp[i][j]=temp;
                    
//                     std::cout << "i= " << i << " " << "j= " << j << " " <<  "c1= " << temp[0] << " " << "c2= " << temp[1] << " " << "c3= " << temp[2] << " " << "c4= " << temp[3] << " " << "c5= " << temp[4] << " " << "c6= " << temp[5] << " " "\n";
                    
                }
                
            }
        }
        
//         std::cout << "interaccion" << "\n" << std::flush;

//         for(int i=0; i<xdim; i++){
//             for(int j=0; j<ydim; j++){
//                 initialize(s->lattice[i][j]);
//             }
//         }
        
        //migration step
        for(int i=0; i<xdim; i++){
            for(int j=0; j<ydim; j++){
//                 if(WELLRNG1024a()<syncp){
                s->lattice[(i+1)%xdim][j][0]=s->lattice_temp[i][j][0];
                s->lattice[(i+1)%xdim][(j-1+ydim)%ydim][1]=s->lattice_temp[i][j][1];
                s->lattice[i][(j-1+ydim)%ydim][2]=s->lattice_temp[i][j][2];
                s->lattice[(i-1+xdim)%xdim][j][3]=s->lattice_temp[i][j][3];
                s->lattice[(i-1+xdim)%xdim][(j+1)%ydim][4]=s->lattice_temp[i][j][4];
                s->lattice[i][(j+1)%ydim][5]=s->lattice_temp[i][j][5];
//                 }else{
//                     s->lattice[i][j][0]+=s->lattice_temp[i][j][0];
//                 s->lattice[i][j][1]+=s->lattice_temp[i][j][1];
//                 s->lattice[i][j][2]+=s->lattice_temp[i][j][2];   //uncomment for asynchronous
//                 s->lattice[i][j][3]+=s->lattice_temp[i][j][3];
//                 s->lattice[i][j][4]+=s->lattice_temp[i][j][4];
//                 s->lattice[i][j][5]+=s->lattice_temp[i][j][5];
//                 }
                for(int k=0; k<NODES; k++){     //clear post interaction configurations
                    s->lattice_temp[i][j][k]=0;
                }
                
            }
        }
        
//         for(int i=1; i<xdim; i++){
//             for(int j=0; j<ydim; j++){
// //                 std::cout << "i= " << i << " " << "j= " << j << " " << "c1= " << s->lattice[i][j][0] << " " << "c2= " << s->lattice[i][j][1] << " " << "c3= " << s->lattice[i][j][2] << " " << "c4= " << s->lattice[i][j][3] << " " << "c5= " << s->lattice[i][j][4] << " " << "c6= " << s->lattice[i][j][5] << " ";
//             }
// //             std::cout << "\n";
//         }
            
//routine if a histogram of particles per node is desired
#ifdef HISTOG
        if(t==tsteps-1){
	  std::stringstream filename;
	  filename << "nemfhistogram_s" << sens << "_d" << dens << ".csv"; 
	  std::ofstream histogram;
	  histogram.open(filename.str(),std::ios::out | std::ios::app);
	  for(int i=0; i<xdim; i++){
	    for(int j=0; j<ydim; j++){
	      
	      if(!empty(s->lattice[i][j])){
		int particles=0;
		for(int k=0; k<NODES; k++){
		  particles+=s->lattice[i][j][k];
		}
		histogram << double(particles)/double(NODES) << ", ";
	      }
	      
	    }
	  }
	  histogram.close();
	}
#endif

        //calculation of polar, nematic, entropy order parameters
	double x1=0, x2=0, y1=0, y2=0,nloc=0,tentrop=0;
	for(int i=0; i<xdim; i++){
            for(int j=0; j<ydim; j++){
                nloc=0;
                if(!empty(s->lattice[i][j])){
		  
		  for(int k=0; k<NODES; k++){
		    
		    x1+=s->lattice[i][j][k]*s->ci[k][0];
		    y1+=s->lattice[i][j][k]*s->ci[k][1];
		    x2+=s->lattice[i][j][k]*(2*sqr(s->ci[k][0])-1);
		    y2+=s->lattice[i][j][k]*2*s->ci[k][0]*s->ci[k][1];
                    nloc+=s->lattice[i][j][k];
		  }
		  tentrop+=(nloc/double(nuparts))*log(nloc/double(nuparts));
		}
	    }
	}
	polaro[t]+=sqrt(sqr(x1)+sqr(y1))/double(nuparts);
	sqpolaro[t]+=(sqr(x1)+sqr(y1))/sqr(double(nuparts));
	nemato[t]+=sqrt(sqr(x2)+sqr(y2))/double(nuparts);
	sqnemato[t]+=(sqr(x2)+sqr(y2))/sqr(double(nuparts));
        entro[t]-=tentrop;
        sqentro[t]+=sqr(tentrop);
        
        if(polmax[t]<sqrt(sqr(x1)+sqr(y1))/double(nuparts)) polmax[t]=sqrt(sqr(x1)+sqr(y1))/double(nuparts);
        if(nemmax[t]<sqrt(sqr(x2)+sqr(y2))/double(nuparts)) nemmax[t]=sqrt(sqr(x2)+sqr(y2))/double(nuparts);
        if(entromax[t]>-1*tentrop || r==0) entromax[t]=-1*tentrop;
        
        
        //band order parameter calculation
//         int ixx,jxx;
        double x3=0,y3=0,nband=0,preband=0,nsite=0;
//         do{
//             ixx=int(round(WELLRNG1024a()*(xdim-1)));
//             jxx=int(round(WELLRNG1024a()*(ydim-1)));
//         }
//         while(empty(s->lattice[ixx][jxx]));
        
        for(int ixx=0; ixx<xdim; ixx++){
            for(int jxx=0; jxx<ydim; jxx++){
                if(!empty(s->lattice[ixx][jxx])){
                    nband=0;
                    x3=0;
                    y3=0;
                    nsite=0;
        for(int k=1; k<5; k++){
            for(int mx=1; mx<NODES; mx++){
                x3+=s->lattice[(ixx+k)%xdim][jxx][mx]*(2*sqr(s->ci[0][0])-1);
                x3+=s->lattice[(ixx+k)%xdim][(jxx-k+ydim)%ydim][mx]*(2*sqr(s->ci[1][0])-1);
                x3+=s->lattice[ixx][(jxx-k+ydim)%ydim][mx]*(2*sqr(s->ci[2][0])-1);
                x3+=s->lattice[(ixx-k+xdim)%xdim][jxx][mx]*(2*sqr(s->ci[3][0])-1);
                x3+=s->lattice[(ixx-k+xdim)%xdim][(jxx+k)%ydim][mx]*(2*sqr(s->ci[4][0])-1);
                x3+=s->lattice[ixx][(jxx+k)%ydim][mx]*(2*sqr(s->ci[5][0])-1);
                
                y3+=s->lattice[(ixx+k)%xdim][jxx][mx]*2*s->ci[0][0]*s->ci[0][1];
                y3+=s->lattice[(ixx+k)%xdim][(jxx-k+ydim)%ydim][mx]*2*s->ci[1][0]*s->ci[1][1];
                y3+=s->lattice[ixx][(jxx-k+ydim)%ydim][mx]*2*s->ci[2][0]*s->ci[2][1];
                y3+=s->lattice[(ixx-k+xdim)%xdim][jxx][mx]*2*s->ci[3][0]*s->ci[3][1];
                y3+=s->lattice[(ixx-k+xdim)%xdim][(jxx+k)%ydim][mx]*2*s->ci[4][0]*s->ci[4][1];
                y3+=s->lattice[ixx][(jxx+k)%ydim][mx]*2*s->ci[5][0]*s->ci[5][1];
                
                nband+=s->lattice[(ixx+k)%xdim][jxx][mx]+s->lattice[(ixx+k)%xdim][(jxx-k+ydim)%ydim][mx]+s->lattice[ixx][(jxx-k+ydim)%ydim][mx]+s->lattice[(ixx-k+xdim)%xdim][jxx][mx]+s->lattice[(ixx-k+xdim)%xdim][(jxx+k)%ydim][mx]+s->lattice[ixx][(jxx+k)%ydim][mx];
                
                nsite+=s->lattice[ixx][jxx][mx];
            }
        }
        if(nband>0) preband+=nsite*sqrt(sqr(x3)+sqr(y3))/(4*double(nband));
            }
            }
        }
        
//         if(nband>0){
        band[t]+=preband/double(nuparts);
        sqband[t]+=sqr(preband/double(nuparts));
//         }
        if(bandmax[t]<preband/double(nuparts)) bandmax[t]=preband/double(nuparts);
//         std::cout << "migracion" << "\n" << std::flush;
        
        //coloring the graphic output of the lattice
#ifdef CIMG
        p.fill(0);
        for(int i=0;i<xdim;i++){
            for(int j=0;j<ydim;j++){
                white[0]=0;     //intiialization with black color (RGB entries)
                white[1]=0;
                white[2]=0;
            }
        }
        
        std::vector<double> val;
        double norma;
        val.resize(3);
        for(int i=0;i<xdim;i++){
            for(int j=0; j<ydim; j++){
                if(!empty(s->lattice[i][j])){
               /*COLOR BY ORIENTATION {*/    
//                     for(int w=0; w<3; w++){
//                         val[w]=0;
//                     }
//                     
//                     norma=0;
//                     
//                     for(int k=0; k<NODES; k++){
//                         if(s->lattice[i][j][k]){
//                             
//                                /*count if there are particles in certain channel/*
//                             norma+=1;
//                            
//                                /*color by color wheel (RGB)*/
//                             if(k<2 || k==5) val[0]+=255;
//                             if((k>0) && (k<4)) val[1]+=255;
//                             if(k>2) val[2]+=255;
//                             
//                         }
//                     }
//                     
                    /*normalization of the RGB entries for maximum saturation! ;) */
//                     for(int w=0; w<3; w++){
// //                         val[w]=val[w]*norma/6.0;
// 		      val[w]=val[w]/norma;
//                     }
                     /*COLOR BY ORIENTATION }*/
		     
 		    /*COLOR BY DENSITY {*/
		    
		    norma=0;
		    
                    /*calculate number of particles in node*/
		    for(int k=0; k<NODES; k++){
		      norma+=s->lattice[i][j][k];
		    }
		    
		    /*maximum brightness at mean density value*/
		    norma=norma/(double)(NODES*dens)*100+155;
		    
                    /*RGB normalization*/
		    val[0]=norma>510?255:norma/2.0;
		    val[1]=norma>510?255:norma/2.0;
		    val[2]=norma>255?255:norma;
		    
		    /*COLOR BY DENSITY {*/
                    
                }else{
                    val[0]=0;
                    val[1]=0;
                    val[2]=0;
                }
                
                white[0]=val[0];
                white[1]=val[1];
                white[2]=val[2];
                
                /*final drawing of lattice*/
                for(int ip=0;ip<4;ip++){
                    for(int jp=0;jp<4;jp++){
	      /* if(df)*/ p.draw_point((const unsigned int)(4*((i+(j/2))%xdim)+ip+2*(j%2)),(const unsigned int)(4*j+jp),white);
                    }
                }
                
            }
        }
        
        dis.display(p);
        /*this for saving the lattice as an image*/
        std::stringstream filename;
         if(t<10)           filename << "cbpol-" << sens << "_d-" << dens << "_t" << "-000" << t << ".bmp";
        if(t>=10&&t<100)   filename << "cbpol-" << sens << "_d-" << dens << "_t" << "-00"  << t << ".bmp";
        if(t>=100&&t<1000) filename << "cbpol-" << sens << "_d-" << dens << "_t" << "-0"  << t << ".bmp";
 
        if(t>=1000)         filename << "cbpol-" << sens << "_d-" << dens << "_t" << "-"    << t << ".bmp";
        if(t==tsteps-1) p.save_bmp(filename.str().c_str());
        
#endif
        
    }
    
}

/*text output of observable*/
// std::cout << "SENS POLAR_ORDER POerr NEMATIC_ORDER NOerr" << "\n";
std::cout << "TIME POLAR_ORDER POerr NEMATIC_ORDER NOerr SPATIAL_ORDER SOerr BAND_ORDER BOerr" << "\n";
// std::cout << "TIME MAX_POLAR_ORDER MAX_NEMATIC_ORDER MAX_SPATIAL_ORDER MAX_BAND_ORDER" << "\n";

for (int t=1; t<tsteps+1; t++){
  
//   if (t-1>(tsteps/2)){
//     std::cout << heating*(2*(tsteps/2)-(t-1)) << " ";
//   }else{
//     std::cout << (heating*(t-1)) << " ";
//   }
  std::cout << double(4*D*t)/double(sqr(v)) << " ";
  std::cout << polaro[t-1]/iters << " ";
  std::cout << sqrt(fabs(sqr(polaro[t-1]/iters)-sqpolaro[t-1]/iters)/iters) << " ";
  std::cout << nemato[t-1]/iters << " ";
  std::cout << sqrt(fabs(sqr(nemato[t-1]/iters)-sqnemato[t-1]/iters)/iters) << " ";
  std::cout << 1.0-(entro[t-1]/iters)/log(xdim*ydim) << " ";
  std::cout << sqrt(fabs(sqr(entro[t-1]/iters)-sqentro[t-1]/iters)/iters)/log(xdim*ydim) << " ";
  std::cout << band[t-1]/iters << " ";
  std::cout << sqrt(fabs(sqr(band[t-1]/iters)-sqband[t-1]/iters)/iters) << "\n";
  
//   std::cout << polmax[t] << " ";
//   std::cout << nemmax[t] << " ";
//   std:: cout << 1.0-(entromax[t]/log(xdim*ydim)) << " ";
//   std::cout << bandmax[t] << "\n";
  
}

// std::cout << sens*dens << " " << polaro[tsteps-1]/iters << " " << sqrt(fabs(sqr(polaro[tsteps-1]/iters)-sqpolaro[tsteps-1]/iters)/iters) << " " << nemato[tsteps-1]/iters << " " << sqrt(fabs(sqr(nemato[tsteps-1]/iters)-sqnemato[tsteps-1]/iters)/iters) << " " << 1.0-(entro[tsteps-1]/iters)/log(xdim*ydim) << " " << sqrt(fabs(sqr(entro[tsteps-1]/iters)-sqentro[tsteps-1]/iters)/iters)/log(xdim*ydim) << " " << band[tsteps-1]/iters << " " << sqrt(fabs(sqr(band[tsteps-1]/iters)-sqband[tsteps-1]/iters)/iters) << "\n" ;


}
