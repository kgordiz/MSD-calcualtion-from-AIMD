#include <stdlib.h>
#include <fstream>
#include <iostream>

using namespace std;

int main(){
int m,n,Natoms,Natoms_O,Ntimesteps;
string line;

Natoms=39;
Natoms_O=23;

Ntimesteps=274785;

ofstream writefile ("r.txt");
ifstream readfile ("XDATCAR_40atoms_1100C");
//writefile.precision(10);
//readfile.precision(10);
if (!readfile.is_open()) {
	cout<<"Unable to open the file!"<<endl;
	exit(1);}
else cout<<"The file could be opened successfully!"<<endl;

for (m=0;m<7;m++){
  getline(readfile,line); //drop the first 7 lines
}
for (m=0;m<Ntimesteps;m++){
  getline(readfile,line); //Drop the configuration number line
  for (n=0;n<Natoms-Natoms_O;n++){ //Drop the lines belonging to other atoms
    getline(readfile,line);
  }
  for (n=0;n<Natoms_O;n++){ //Get the lines belonging to O
    getline(readfile,line);
    writefile<<line<<endl;
  }
}

readfile.close();
writefile.close();

return 21;
}
