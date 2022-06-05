/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2014 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "CLTool.h"
#include "CLToolRegister.h"
#include "wrapper/Plumed.h"
#include "tools/Vector.h"
#include "tools/Matrix.h"
#include "tools/Random.h"
#include <string>
#include <cstdio>
#include <cmath>
#include <vector>

using namespace std;

namespace PLMD{
namespace cltools{


//+PLUMEDOC TOOLS pesmd
/*
Pesmd allows one to do (biased) Langevin dynamics on a two-dimensional potential energy surface.
The input to pesmd is specified in an input file. 
The directives available are as follows:

\par Examples

You run a Langevin simulation using pesmd with the following command:
\verbatim
plumed pesmd < input
\endverbatim

The following is an example of an input file for a pesmd simulation. This file
instructs pesmd to do 50 steps of Langevin dynamics on a 2D potential energy surface
at a temperature of 0.722 
\verbatim
pesfile potential.dat
temperature 0.722
tstep 0.005
friction 1
nstep 50
posx_init 0.0 
posy_init 0.0 
\endverbatim

If you run the following a description of all the directives that can be used in the
input file will be output.
\verbatim
plumed pesmd --help
\endverbatim

The potential energy surface is read from a file, via the directive "pesfile", 
with the following format. A first header lines containing "#", and two integers that
indicate the number of gridpoints in the x and y directions. The header line is followed
by three columns, with on the first two columns the evenly spaced (x,y) gridpoints and 
the third column containing the potential value. Every block of lines at constant x, 
and y running from its minimum to its maximum, is followed by an empty line (this is the
gnuplot format for contour plots). Derivatives are estimated numerically. 

The plumed machinery is invoked through a plumed file "plumed.dat", in which the dynamics
may be probed or biased via the components of a DISTANCE collective variable applied to
atoms one and two. The following example will print the position of the Langevin particle as
the "d1" distance components to the "colvar.out" file, every 100 dynamics steps.
\verbatim
UNITS LENGTH=nm TIME=fs

d1: DISTANCE ATOMS=1,2 COMPONENTS

PRINT ARG=d1.x,d1.y STRIDE=100   FILE=colvar.out

\endverbatim



*/
//+ENDPLUMEDOC

class PesMD:
public PLMD::CLTool
{
  string description()const{
    return "run pesmd code";
  }


public:
static void registerKeywords( Keywords& keys ){ 
  keys.add("compulsory","nstep","The number of steps of dynamics you want to run");
  keys.add("compulsory","temperature","NVE","the temperature at which you wish to run the simulation in LJ units");
  keys.add("compulsory","friction","off","The friction (in LJ units) for the langevin thermostat that is used to keep the temperature constant");
  keys.add("compulsory","tstep","0.005","the integration timestep in LJ units");
  keys.add("compulsory","pesfile","A data file containing the potential surface to do dynamics on");  
  keys.add("compulsory","posx_init","0","Starting x position");
  keys.add("compulsory","posy_init","0","Starting y position");
  keys.add("compulsory","idum","0","The random number seed");
  keys.add("compulsory","wrapatoms","false","If true, atomic coordinates are written wrapped in minimal cell");
}

PesMD( const CLToolOptions& co ) :
  CLTool(co)
{
  inputdata=ifile;
}

private:

/*------------------------------------------------------------------------*/
  void read_input(double& temperature,
		  double& tstep,
		  double& friction,
		  int&    nstep,
		  string& pesfile,
		  double& posx_init,
		  double& posy_init,
		  bool&   wrapatoms,
		  int&    idum)
{

  // Read everything from input file
  std::string tempstr; parse("temperature",tempstr);
  if( tempstr!="NVE" ) Tools::convert(tempstr,temperature);
  parse("tstep",tstep);
  std::string frictionstr; parse("friction",frictionstr);
  if( tempstr!="NVE" ){
      if(frictionstr=="off"){ fprintf(stderr,"Specify friction for thermostat\n"); exit(1); }
      Tools::convert(frictionstr,friction); 
  } 
  parse("posx_init",posx_init);
  parse("posy_init",posy_init);
  parse("nstep",nstep);
  parse("idum",idum);

  // Read in stuff with sanity checks
  parse("pesfile",pesfile);
  if(pesfile.length()==0){
      fprintf(stderr,"Specify pes file\n");
      exit(1);
  }  
  std::string w;
  parse("wrapatoms",w);
  wrapatoms=false;
  if(w.length()>0 && (w[0]=='T' || w[0]=='t')) wrapatoms=true;
}

/*------------------------------------------------------------------------*/
void read_gridsize(const string & pesfile,vector<int>& ngrid){
// read the size of the grid from pesfile "potential.dat"
  FILE* fp=fopen(pesfile.c_str(),"r");
  if(!fp){
    fprintf(stderr,"ERROR: file %s not found\n",pesfile.c_str());
    exit(1);
  }
  char hash[2];
  fscanf(fp,"%s%1000d%1000d",hash,&ngrid[0],&ngrid[1]);
  fclose(fp);
}

/*------------------------------------------------------------------------*/
void read_potential(const string & pesfile,vector<int>& ngrid,vector<double>& x1pot,
		    vector<double>& x2pot,Matrix<double>& extpot){
  char line[128];
  FILE* fp=fopen(pesfile.c_str(),"r");

  if(!fp){
    fprintf(stderr,"ERROR: file %s not found\n",pesfile.c_str());
    exit(1);
  }
  char hash[2];
  fscanf(fp,"%s%1000d%1000d",hash,&ngrid[0],&ngrid[1]);
  for(int i=0;i<ngrid[0];++i){
    for(int j=0;j<ngrid[1];++j){
      fscanf(fp,"%lf%lf%lf",&x1pot[i],&x2pot[j],&extpot[i][j]);
    }
    fgets(line,128,fp);
  }
  fclose(fp);
}

/*------------------------------------------------------------------------*/
void get_derivatives(vector<int>& ngrid,vector<double> x1pot,vector<double> x2pot,Matrix<double> extpot,Matrix<double>& extdvx1,Matrix<double>& extdvx2,Matrix<double>& extdvxx){

  int i, j, nx1, nx2;

  nx1 = ngrid[0] - 1;
  nx2 = ngrid[1] - 1;
  for(i=0;i<ngrid[0];i++){
    extdvx1[i][0] =  extdvx1[i][nx2] = 0.;
    extdvx2[i][0] =  extdvx2[i][nx2] = 0.;
    extdvxx[i][0] =  extdvxx[i][nx2] = 0.;
  }
  for(j=0;j<ngrid[1];j++){
    extdvx1[0][j] =  extdvx1[nx1][j] = 0.;
    extdvx2[0][j] =  extdvx2[nx1][j] = 0.;
    extdvxx[0][j] =  extdvxx[nx1][j] = 0.;
  }
  for(i=1;i<nx1;i++){
    for(j=1;j<nx2;j++){
      extdvx1[i][j] = (extpot[i+1][j] - extpot[i-1][j]) / (x1pot[i+1] - x1pot[i-1]);
      extdvx2[i][j] = (extpot[i][j+1] - extpot[i][j-1]) / (x2pot[j+1] - x2pot[j-1]);
      extdvxx[i][j] = (extpot[i+1][j+1] - extpot[i+1][j-1] - extpot[i-1][j+1] + extpot[i-1][j-1]) 
        / ((x1pot[i+1] - x1pot[i-1]) * (x2pot[j+1] - x2pot[j-1]));        
    }
  }
}


/*------------------------------------------------------------------------*/
void randomize_velocities(const int natoms,const int ndim,const double temperature,const vector<double>&masses,vector<Vector>& velocities,Random&random){
// randomize the velocities according to the temperature
  for(int iatom=0;iatom<natoms;iatom++) for(int i=0;i<ndim;i++)
      velocities[iatom][i]=sqrt(temperature/masses[iatom])*random.Gaussian();
}

void zero_forces(const int natoms,vector<Vector>& forces){

  for(int i=0;i<natoms;++i){
    forces[i][0] = forces[i][1] = forces[i][2] = 0.0;
  }
}

void pbc(const double cell[3],const Vector & vin,Vector & vout){
// apply periodic boundary condition to a vector
  for(int i=0;i<3;i++){
    vout[i]=vin[i]-floor(vin[i]/cell[i]+0.5)*cell[i];
  }
}

/*------------------------------------------------------------------------*/
void bcucof(std::vector<double>& y,std::vector<double>& y1,std::vector<double>& y2,
	    std::vector<double>& y12, double d1, double d2,
	    Matrix<double>& c){

  static int wt_d[16*16]={ 1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			  0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,
			  -3,0,0,3,0,0,0,0,-2,0,0,-1,0,0,0,0,
			  2,0,0,-2,0,0,0,0,1,0,0,1,0,0,0,0,
			  0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
			  0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,
			  0,0,0,0,-3,0,0,3,0,0,0,0,-2,0,0,-1,
			  0,0,0,0,2,0,0,-2,0,0,0,0,1,0,0,1,
			  -3,3,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0,
			  0,0,0,0,0,0,0,0,-3,3,0,0,-2,-1,0,0,
			  9,-9,9,-9,6,3,-3,-6,6,-6,-3,3,4,2,1,2,
			  -6,6,-6,6,-4,-2,2,4,-3,3,3,-3,-2,-1,-1,-2,
			  2,-2,0,0,1,1,0,0,0,0,0,0,0,0,0,0,
			  0,0,0,0,0,0,0,0,2,-2,0,0,1,1,0,0,
			  -6,6,-6,6,-3,-3,3,3,-4,4,2,-2,-2,-2,-1,-1,
			  4,-4,4,-4,2,2,-2,-2,2,-2,-2,2,1,1,1,1};
  int l=0,k,j,i,wt[16][16];
  double xx,d1d2,cl[16],x[16];

  // This is to set up the coefficient matrix
  for (unsigned i=0;i<16;i++) for (unsigned j=0;j<16;j++){ wt[i][j]=wt_d[l++]; }

  d1d2=d1*d2;
  for (i=1;i<=4;i++) {
    x[i-1]=y[i];
    x[i+3]=y1[i]*d1;
    x[i+7]=y2[i]*d2;
    x[i+11]=y12[i]*d1d2;
  }
  for (i=0;i<=15;i++) {
    xx=0.0;
    for (k=0;k<=15;k++) xx += wt[i][k]*x[k];
    cl[i]=xx;
  }
  l=0;
  for (i=1;i<=4;i++)
    for (j=1;j<=4;j++) c[i][j]=cl[l++];
}
  
/*------------------------------------------------------------------------*/
void bcuint(std::vector<double>& y,std::vector<double>& y1,std::vector<double>& y2,
	    std::vector<double>& y12, double x1l, double x1u, double x2l, double x2u, 
	    double x1, double x2, double & ansy,double & ansy1, double & ansy2){


  int i;
  double t,u,d1,d2;
  Matrix<double> c(5,5);

  //  c=dmatrix(1,4,1,4);  /* allocate a matrix that runs from 1-4, 1-4  */
  d1=x1u-x1l;
  d2=x2u-x2l;
  bcucof(y,y1,y2,y12,d1,d2,c);
  if (x1u == x1l || x2u == x2l) {
    fprintf(stderr,"Bad input in routine bcuint");
    exit(1);
  } 
  t=(x1-x1l)/d1;
  u=(x2-x2l)/d2;
  ansy = ansy2 = ansy1 = 0.0;
  for (i=4;i>=1;i--) {
    ansy=t*ansy+((c[i][4]*u+c[i][3])*u+c[i][2])*u+c[i][1];
    ansy2=t*ansy2+(3.0*c[i][4]*u+2.0*c[i][3])*u+c[i][2];
    ansy1=u*ansy1+(3.0*c[4][i]*t+2.0*c[3][i])*t+c[2][i];
  }
  ansy1 /= d1;
  ansy2 /= d2;
  //  free_dmatrix(c,1,4,1,4);
}

/*------------------------------------------------------------------------*/
double compute_extpot_force(const vector<Vector>& positions,
			    vector<Vector>& forces,vector<int>& ngrid, vector<double>& x1pot,
			    vector<double>& x2pot,Matrix<double>& extpot,Matrix<double>& extdvx1,
			    Matrix<double>& extdvx2,Matrix<double>& extdvxx){


  int i, j;
  std::vector<double> y(5), y1(5), y2(5), y3(5);
  double v, vdx1, vdx2;

  /* compute the force due to the 2D external potential on atom 1 (atom 0 is fixed at the origin) */
  double px1=positions[1][0];
  double px2=positions[1][1];

  /* find in which grid square pcolvar is */
  if((px1 < x1pot[0]) || (px1 > x1pot[ngrid[0]-1]) || (px2 < x2pot[0]) || (px2 > x2pot[ngrid[1]-1])){
    /* pcolvar is off the electrostatic potential grid where the external force is zero */
    return 0.0;
  }

  i = j = 1;
  while(px1 > x1pot[i]) i++;
  while(px2 > x2pot[j]) j++;
  i--;
  j--;

  y[1] = extpot[i][j];
  y[2] = extpot[i+1][j];
  y[3] = extpot[i+1][j+1];
  y[4] = extpot[i][j+1];

  y1[1] = extdvx1[i][j];
  y1[2] = extdvx1[i+1][j];
  y1[3] = extdvx1[i+1][j+1];
  y1[4] = extdvx1[i][j+1];

  y2[1] = extdvx2[i][j];
  y2[2] = extdvx2[i+1][j];
  y2[3] = extdvx2[i+1][j+1];
  y2[4] = extdvx2[i][j+1];

  y3[1] = extdvxx[i][j];
  y3[2] = extdvxx[i+1][j];
  y3[3] = extdvxx[i+1][j+1];
  y3[4] = extdvxx[i][j+1];

  bcuint(y,y1,y2,y3,x1pot[i],x1pot[i+1],x2pot[j],x2pot[j+1],px1,px2,v,vdx1,vdx2);

#define DEBUG_OFF
#ifdef DEBUG
  printf("punt %d (%d,%d): x=%lf y=%lf E=%lf (%lf %lf %lf)\n",1,i,j,x1pot[i],x2pot[j],y[1],y1[1],y2[1],y3[1]);
  printf("punt %d (%d,%d): x=%lf y=%lf E=%lf (%lf %lf %lf)\n",2,i+1,j,x1pot[i+1],x2pot[j],y[2],y1[2],y2[2],y3[2]);
  printf("punt %d (%d,%d): x=%lf y=%lf E=%lf (%lf %lf %lf)\n",3,i+1,j+1,x1pot[i+1],x2pot[j+1],y[3],y1[3],y2[3],y3[3]);
  printf("punt %d (%d,%d): x=%lf y=%lf E=%lf (%lf %lf %lf)\n",4,i,j+1,x1pot[i],x2pot[j+1],y[4],y1[4],y2[4],y3[4]);
  printf("Potential at (%lf,%lf) = %g; Forces are %g, %g\n",px1,px2,v,vdx1,vdx2);
#endif

  forces[1][0]  -= vdx1;
  forces[1][1]  -= vdx2;

  return v;

}


/*------------------------------------------------------------------------*/
void compute_engkin(const int natoms,const vector<double>& masses,const vector<Vector>& velocities,double & engkin)
{
// calculate the kinetic energy from the velocities
  engkin=0.0;
  for(int iatom=0;iatom<natoms;iatom++)for(int k=0;k<3;k++){
    engkin+=0.5*masses[iatom]*velocities[iatom][k]*velocities[iatom][k];
  }
}


/*------------------------------------------------------------------------*/
void thermostat(const int natoms,const int ndim,const vector<double>& masses,const double dt,const double friction,
                const double temperature,vector<Vector>& velocities,double & engint,Random & random){
// Langevin thermostat, implemented as decribed in Bussi and Parrinello, Phys. Rev. E (2007)
// it is a linear combination of old velocities and new, randomly chosen, velocity,
// with proper coefficients
  double c1=exp(-friction*dt);
  for(int iatom=0;iatom<natoms;iatom++){
    double c2=sqrt((1.0-c1*c1)*temperature/masses[iatom]);
    for(int i=0;i<ndim;i++){
      engint+=0.5*masses[iatom]*velocities[iatom][i]*velocities[iatom][i];
      velocities[iatom][i]=c1*velocities[iatom][i]+c2*random.Gaussian();
      engint-=0.5*masses[iatom]*velocities[iatom][i]*velocities[iatom][i];
    }
  }
}



/*------------------------------------------------------------------------*/
virtual int main(FILE* in,FILE*out,PLMD::Communicator& pc){
  int            natoms=2;     // number of atoms
  int            ndim;         // dimensionality of the system (1, 2, or 3)
  vector<int>    ngrid;        // size of the pes grid
  vector<Vector> positions;    // atomic positions
  vector<Vector> velocities;   // velocities
  vector<double> masses;       // masses
  vector<Vector> forces;       // forces
  double         posx_init=0.0;  // starting position on the PES
  double         posy_init=0.0;  // starting position on the PES
  double         cell[3];      // cell size
  double         cell9[3][3];  // cell size


// input parameters
// all of them have a reasonable default value, set in read_input()
  double      tstep;             // simulation timestep
  double      temperature;       // temperature
  double      friction;          // friction for Langevin dynamics (for NVE, use 0)
  int         nstep;             // number of steps
  int         idum;              // seed
  int         plumedWantsToStop; // stop flag
  bool        wrapatoms;         // if true, atomic coordinates are written wrapped in minimal cell
  string      pesfile;           // name of file with potential energy surface on a grid
  string      string;            // a string for parsing

  double engkin;                 // kinetic energy
  double engpot;                 // external potential energy
  double engconf=0.0;            // configurational energy
  double engint;                 // integral for conserved energy in Langevin dynamics

  Random random;                 // random numbers stream

  PLMD::Plumed* plumed=NULL;

// Commenting the next line it is possible to switch-off plumed
  plumed=new PLMD::Plumed;

  if(plumed){
    int s=sizeof(double);
    plumed->cmd("setRealPrecision",&s);
  }

  read_input(temperature,tstep,friction,nstep,pesfile,posx_init,posy_init,wrapatoms,idum);

// number of atoms is read from file inputfile
  ndim=2; 
  ngrid.resize(ndim);
  read_gridsize(pesfile,ngrid);
  vector<double> x1pot(ngrid[0]);
  vector<double> x2pot(ngrid[1]);
  Matrix<double> extpot(ngrid[0],ngrid[1]);
  Matrix<double> extdvxx(ngrid[0],ngrid[1]);
  Matrix<double> extdvx1(ngrid[0],ngrid[1]);
  Matrix<double> extdvx2(ngrid[0],ngrid[1]);
  read_potential(pesfile,ngrid,x1pot,x2pot,extpot);
  get_derivatives(ngrid,x1pot,x2pot,extpot,extdvx1,extdvx2,extdvxx);

// write the parameters in output so they can be checked
  fprintf(out,"%s %f\n","Temperature                      :",temperature);
  fprintf(out,"%s %f\n","Time step                        :",tstep);
  fprintf(out,"%s %f\n","Friction                         :",friction);
  fprintf(out,"%s %d\n","Number of steps                  :",nstep);
  fprintf(out,"%s %d  ","Grid size                        :",ngrid[0]); 
  for(int i=1;i<ndim;++i)   fprintf(out,"x %d",ngrid[i]); fprintf(out,"\n");
  fprintf(out,"Grid range cv%d                   : %lf <> %lf \n",1,x1pot[0],x1pot[ngrid[0]-1]); 
  fprintf(out,"Grid range cv%d                   : %lf <> %lf \n",2,x2pot[0],x2pot[ngrid[1]-1]); 
  fprintf(out,"%s %d\n","Seed                             :",idum);
  fprintf(out,"%s %s\n","Are atoms wrapped on output?     :",(wrapatoms?"T":"F"));

// Setting the seed
  random.setSeed(idum);

// allocation of dynamical arrays
  natoms=2;
  positions.resize(natoms);
  velocities.resize(natoms);
  forces.resize(natoms);
  masses.resize(natoms);

// set masses and positions 
  masses[0] = 0.0;
  masses[1] = 1.0;
  positions[0][0] = positions[0][1] = 0.0;
  positions[1][0] = posx_init; positions[1][1] = posy_init;

// energy integral initialized to 0
  engint=0.0;


// velocities are randomized according to temperature
  randomize_velocities(natoms,ndim,temperature,masses,velocities,random);
  velocities[0][0] = velocities[0][1] = 0.0;
  //  velocities[1][0] = velocities[1][1] = 0.0;

  if(plumed){
    plumed->cmd("setNoVirial");
    plumed->cmd("setNatoms",&natoms);
    plumed->cmd("setMDEngine","pesMD");
    plumed->cmd("setTimestep",&tstep);
    plumed->cmd("setPlumedDat","plumed.dat");
    int pversion=0;
    plumed->cmd("getApiVersion",&pversion);
// setting kbT is only implemented with api>1
// even if not necessary in principle in PesMD (which is part of plumed)
// we leave the check here as a reference
      if(pversion>1){
        plumed->cmd("setKbT",&temperature);
      }
    plumed->cmd("init");
  }


// forces are computed before starting md
  zero_forces(natoms,forces);
  engpot = compute_extpot_force(positions,forces,ngrid,x1pot,x2pot,extpot,extdvx1,extdvx2,extdvxx);
  //  printf("DEBUG main %d p= %lf %lf | v= %lf %lf | f= %lf %lf | e= %lf\n",0,positions[1][0],positions[1][1],velocities[1][1],velocities[1][0],forces[1][0],forces[1][1],engpot);

// remove forces if ndim<3
  if(ndim<3)
  for(int iatom=0;iatom<natoms;++iatom) for(int k=ndim;k<3;++k) forces[iatom][k]=0.0;

// here is the main md loop
// Langevin thermostat is applied before and after a velocity-Verlet integrator
// the overall structure is:
//   thermostat
//   update velocities
//   update positions
//   (eventually recompute neighbour list)
//   compute forces
//   update velocities
//   thermostat
//   (eventually dump output informations)
  for(int istep=0;istep<nstep;istep++){
    thermostat(natoms,ndim,masses,0.5*tstep,friction,temperature,velocities,engint,random);

    for(int k=0;k<3;k++){
      velocities[1][k]+=forces[1][k]*0.5*tstep/masses[1];
      positions[1][k]+=velocities[1][k]*tstep;
    }
    
    zero_forces(natoms,forces);
    engpot = compute_extpot_force(positions,forces,ngrid,x1pot,x2pot,extpot,extdvx1,extdvx2,extdvxx);
    //    printf("DEBUG main %d p= %lf %lf | v= %lf %lf | f= %lf %lf | e= %lf\n",istep,positions[1][0],positions[1][1],velocities[1][1],velocities[1][0],forces[1][0],forces[1][1],engpot);

    if(plumed){
      int istepplusone=istep+1;
      plumedWantsToStop=0;
      for(int i=0;i<3;i++)for(int k=0;k<3;k++) cell9[i][k]=0.0;
      for(int i=0;i<3;i++) cell9[i][i]=cell[i];
      plumed->cmd("setStep",&istepplusone);
      plumed->cmd("setMasses",&masses[0]);
      plumed->cmd("setForces",&forces[0]);
      plumed->cmd("setEnergy",&engconf);
      plumed->cmd("setPositions",&positions[0]);
      plumed->cmd("setBox",cell9);
      plumed->cmd("setStopFlag",&plumedWantsToStop);
      plumed->cmd("calc");
      if(plumedWantsToStop) nstep=istep;
    }
// remove forces if ndim<3
   if(ndim<3)
    for(int iatom=0;iatom<natoms;++iatom) for(int k=ndim;k<3;++k) forces[iatom][k]=0.0;

   for(int k=0;k<3;k++){
     velocities[1][k]+=forces[1][k]*0.5*tstep/masses[1];
   }

   thermostat(natoms,ndim,masses,0.5*tstep,friction,temperature,velocities,engint,random);

// kinetic energy is calculated
   compute_engkin(natoms,masses,velocities,engkin);

  }

  if(plumed) delete plumed;

  return 0;
}


};

PLUMED_REGISTER_CLTOOL(PesMD,"pesmd")

}
}




