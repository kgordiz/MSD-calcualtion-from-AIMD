#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <unistd.h>
#include <math.h>

using namespace std;

int main(){
  int m,n,dt,Natoms,Ntimesteps,indx,indx_prev,indx_dt;
  double dx,dy,dz,x_size,y_size,z_size,Sum_x,Sum_y,Sum_z,**r_real,**r_cartesian,r_direct[3],**MSD,a[3][3];

  Natoms=23; // Number of mobile atoms (O)
  Ntimesteps=274785;

  a[0][0] = 7.875519;
  a[0][1] = -0.000031;
  a[0][2] = 0.031486;

  a[1][0] = -0.000165;
  a[1][1] = 7.871420;
  a[1][2] = 0.000215;

  a[2][0] = 0.031519;
  a[2][1] = 0.000075;
  a[2][2] = 7.863395;

  x_size = a[0][0];
  y_size = a[1][1];
  z_size = a[2][2];

  r_real = new double* [Ntimesteps*Natoms];
  for (m=0;m<Ntimesteps*Natoms;m++) r_real[m] = new double [3];

  r_cartesian = new double* [Ntimesteps*Natoms];
  for (m=0;m<Ntimesteps*Natoms;m++) r_cartesian[m] = new double [3];

  MSD = new double* [Ntimesteps];
  for (m=0;m<Ntimesteps;m++) MSD[m] = new double [3];

  ifstream readfile ("r.txt");
  ofstream writefile ("MSD_vs_dt_40_1100C.txt");
  writefile.precision(10);
  readfile.precision(10);

  for (m=0;m<Ntimesteps;m++){  //This block calculates r_cartesian from r_direct obtained from XDATCAR
    for (n=0;n<Natoms;n++){
      readfile>>r_direct[0];
      readfile>>r_direct[1];
      readfile>>r_direct[2];
      indx = m*Natoms + n;
      r_cartesian[indx][0]=r_direct[0]*a[0][0]+r_direct[1]*a[1][0]+r_direct[2]*a[2][0];
      r_cartesian[indx][1]=r_direct[0]*a[0][1]+r_direct[1]*a[1][1]+r_direct[2]*a[2][1];
      r_cartesian[indx][2]=r_direct[0]*a[0][2]+r_direct[1]*a[1][2]+r_direct[2]*a[2][2];
    }
  }

  for (m=0;m<Ntimesteps;m++){  //This block calculates r_real from r_cartesian (cancelling the PBC effects)
    if (m==0){  // positions at time step zero
      for (n=0;n<Natoms;n++){
        indx = m*Natoms + n;
        r_real[indx][0] = r_cartesian[indx][0];
        r_real[indx][1] = r_cartesian[indx][1];
        r_real[indx][2] = r_cartesian[indx][2];
      }
    }
    else {
      for (n=0;n<Natoms;n++){
        indx_prev = (m-1)*Natoms + n;
        indx = m*Natoms + n;
        dx = r_cartesian[indx][0] - r_cartesian[indx_prev][0];
        dy = r_cartesian[indx][1] - r_cartesian[indx_prev][1];
        dz = r_cartesian[indx][2] - r_cartesian[indx_prev][2];

        if (dx > x_size*0.5) dx -= x_size;
        else if (dx <= -x_size*0.5) dx += x_size;
        if (dy > y_size*0.5) dy -= y_size;
        else if (dy <= -y_size*0.5) dy += y_size;
        if (dz > z_size*0.5) dz -= z_size;
        else if (dz <= -z_size*0.5) dz += z_size;

        r_real[indx][0] = r_real[indx_prev][0] + dx;
        r_real[indx][1] = r_real[indx_prev][1] + dy;
        r_real[indx][2] = r_real[indx_prev][2] + dz;
      }
    }
  }

  // Loop over all delta_t possibilities
  for (dt=0;dt<Ntimesteps;dt++){
    Sum_x = 0;
    Sum_y = 0;
    Sum_z = 0;
    for (m=0;m<Ntimesteps-dt;m++){
      for (n=0;n<Natoms;n++){
        indx = m*Natoms + n;
        indx_dt = (m+dt)*Natoms + n;
        Sum_x += pow((r_real[indx_dt][0]-r_real[indx][0]),2);
        Sum_y += pow((r_real[indx_dt][1]-r_real[indx][1]),2);
        Sum_z += pow((r_real[indx_dt][2]-r_real[indx][2]),2);
      }
    }
    MSD[dt][0] = Sum_x/(2.0*Natoms*3*(Ntimesteps-dt)); // MSD in x-direction
    MSD[dt][1] = Sum_y/(2.0*Natoms*3*(Ntimesteps-dt)); // MSD in y-direction
    MSD[dt][2] = Sum_z/(2.0*Natoms*3*(Ntimesteps-dt)); // MSD in z-direction
    writefile<<dt<<" "<<MSD[dt][0]<<" "<<MSD[dt][1]<<" "<<MSD[dt][2]<<endl;
  }

  readfile.close();
  writefile.close();
  return 78;
}

