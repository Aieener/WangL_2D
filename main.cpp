// main.cpp

#include "MC.h"
#include <iostream>
using namespace std;

int main()
{
	stringstream t;
	double start = clock();

	// ======================= MCRUN & Plotting the final config ===============================
	array<double,10000>  wf;
	vector<HR> R;
	MC m(1E11L,8,64,64,13);
	wf = m.MCRUN();
	// m.MCWEIGHT();
	// ======================= end of simulation, print out the time =======
	double end = clock();
	cout <<"This simulation takes "<< (double(end-start)/CLOCKS_PER_SEC)<<" sec."<<endl;
	t <<"This simulation takes "<< (double(end-start)/CLOCKS_PER_SEC) <<" sec."<<endl;
	ofstream myfile ("time.txt");
	string datat = t.str();
	myfile << datat;
	myfile.close();	
	return 0;
}