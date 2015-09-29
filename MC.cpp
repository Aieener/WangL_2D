/*
* S2LG.cpp
* Simulation of 2-D Rods (one specie algorithm)By WangL
* Author: Yuding Ai
* Date: 2015.09.24
* *************************** MC implementation ********************************
* This simulation follows Wang-Laudau sampling.
* ******************************************************************************
*/


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "square.h"
#include "cells.h"
#include "MC.h"
#include "hardrods.h"
#include <cstdlib>
#include <cmath>
#include <time.h>
#include <vector>
#include "histogram.h"
#include <array>
using namespace std;


MC::MC(long int ST, int LEN,int C, int R, double Z)
{
	r = R;
	c = C;
	length = LEN;
	step = ST;
	z = Z;
	nh=nv=dh=dv=ah=av=0;
}

// vector<HR> MC::getRodlist() 
// {
// 	return Rodlist;
// }

double MC::getTho() const
{
	double tho;	
	tho = double(length*(av+ah-dv-dh))/double(r*c);
	return tho;
}

double MC::getQ() const
{
	double Q;	
	Q = (nv - nh)/(nv + nh);
	return Q;
}

double MC::getAaccp() const
{
	double A;
	A = z*(double(r*c))/(double(av+ah-dv-dh+1.0)*double(length));
	return A;
}

double MC::getDaccp() const
{
	double D;
	D = (double(av+ah-dv-dh)*double(length))/(z*(double(r*c)));
	return D;
}

double MC::getNh() const
{
	return nh;
}
double MC::getNv() const
{
	return nv;
}

// void MC::setRodlist(std::vector<HR> RodL)
// {
// 	Rodlist = RodL;
// }

void MC::Add(Cells &s,double &prob,double &aaccp,array<double,10000> &WF)
{
	int x,y,o; // pick a random position and orientation for the HR to be added;
	double proba;
	x = rand()%c;
	y = rand()%r;
	o = rand()%2;// change it to 1 for  lattice gas case

	if(s.getSquare(x,y).isEmpty()) // if it's open, try to do Addition;
	{
		HR rod(x,y,length,o);

		//======================== Vertical ===============================

		if(o == 0)
		{
			int counter = 0;

			for (int j = 0; j < length-1; j++)
			{					
				// check if the vertical space is wide open
				if(s.getSquare(x,(y+j+1)%r).isOccupied())
				{
					counter++;
				}					
			}
			if(counter == 0)
			{
				proba = min(1.0,aaccp*(exp(WF[int(nv+1)] - WF[int(nv)])));
				if(prob <= proba)
				{					
					// Do addition;
					// push the new rod into the Rodlist;
					VRodlist.push_back(rod);
					av++;
					nv++;// accumulate the # of ver rod;
					// update new N, E and new config;
					for (int i = 0; i < length; i++)
					{	
						s.getSquare(x,(y+i)%r).setStatus(1);
					}
				}
			}								
		}

		else 
		{
        //======================= Horizontal  ============================
			int counter = 0;
			for (int j = 0; j< length-1 ; j++)
			{
				// check if the horizontal space is wide open
				if(s.getSquare((x+1+j)%c,y).isOccupied())
				{
					counter++;
				}							
			}
			if (counter == 0)
			{
				proba = min(1.0,aaccp);
				if(prob <= proba)
				{
					//Do addition;
					//push the new rod into the Rodlist;
					HRodlist.push_back(rod);
					ah++;
					nh++;// accumulate the # of hor rod;

					// update new N, E and new config;
					for (int i = 0; i < length; i++)
					{
						s.getSquare((x+i)%c,y).setStatus(1);
					}
				}	
			}							
		}
    }
}

void MC::Del(Cells &s,double &prob,double &daccp,double &size,array<double,10000> &WF)
{
	double probd;

	if(nv+nh > 0)// make sure there are rod;
	{
		int indx; // pick a random index from the Rodlist;
		indx = rand()%int(nv+nh);

		//remove Rodlist[indx];
		int x,y;// the position of the target on the cells;

		if (indx < nv) // vertical
		{
			x = VRodlist[indx].getX();
			y = VRodlist[indx].getY();

			probd = min(1.0,daccp*(exp(WF[int(nv-1)] - WF[int(nv)])));
			if(prob <= probd)
			{
				// --------------------- it's a vertical rod -----------------------			
				for(int i = 0; i<VRodlist[indx].getLength(); i++)
				{
					// update the new config of cells
					s.getSquare(x,(y+i)%r).setStatus(0);
				}
				// remove the target rod from the vector Rodlist;
				VRodlist.erase(VRodlist.begin() + indx);
				nv--;// substract the # of ver rod;
				dv++;					
			}
		}

		else
		{
			x = HRodlist[indx - nv].getX();
			y = HRodlist[indx - nv].getY();

			probd = min(1.0,daccp);
			if(prob <= probd)
			{
			// --------------------- it's a Horizontal rod -----------------------
				for(int i = 0; i<HRodlist[indx-nv].getLength(); i++)
				{
					// update the new config of cells
					s.getSquare((x+i)%c,y).setStatus(0);
				}
				// remove the target rod from the vector Rodlist;
				HRodlist.erase(HRodlist.begin() + indx - nv);
				nh--;// substract the # of hor rod;
				dh++;

			}				
		}										
	}
}


array<double,10000>  MC::MCRUN()
{
	Cells s(c,r); //  setting the lattice;
    
    //==========================================================  declare instance variables ============================================================= //
	stringstream sh;
	sh.precision(20);
	double addordel;           // the prob to decide either add or del;
	// double probd,proba;      // the acceptance prob of addition and deletion; 
	double prob;               // the prob to decide either accept add/del;
	double aaccp,daccp;      // the acceptance probabilities: 
	double V = double(r*c);    // the total lattice size
	double K = double(length); //
    array<double,10000> WF;
    double g = 1;
		
	srand(time(NULL));
	long int i = 0;
	Histogram histotal(0,0.8*V/K,4); // take 80% of the full range.
	Histogram histotalacc(0,0.8*V/K,4); // take 80% of the full range.


	// int av = 0; 
	// int ah = 0; // an interger keep track of the times of we reset the histogram

	// =============================================================Start MC runs ======================================================================== //
	while (g>=1E-7)
	{
		i++;
		// generate a random probability to decide either add or del;
		addordel = rand()%2;
		double size = nv+nh;
		prob = ((double) rand() / (RAND_MAX)); 

		aaccp = V/((size+1.0)*K);
		daccp = (size*K)/V;	

		// probd = min(1.0,daccp);
		// proba = min(1.0,aaccp);

        // ===========================Addition ===================================
		if(addordel == 0) 
		{
			if(size <= 0.8*V/K) // make sure does not go beyond the histogram
			{
				//Do Addition;
				Add(s,prob,aaccp,WF);
			}
	        WF[int(size)] -= g;

		}

		// ============================Deletion=============================
		else
		{
			if (size != 0) // make sure there are rods to be del;
			{
				//Do deletion;
				Del(s,prob,daccp,size,WF);
	            WF[int(size)] -= g;
			}			
		}

		// ======================= Record the datas =============================================

		//================================================ load data into histogram ===================================================================//
		histotal.record(nv);
		histotalacc.record(nv);
	
		//************************************************ check if the current histogram nv is "flat enough" *********************************************//
		double hisvmin,hisvmean;
		hisvmin = double(histotal.Minave().first);
		hisvmean = double(histotal.Minave().second);

		if(hisvmin/hisvmean >=0.8)
		{   
			g = 0.5*g;
			histotal.reset();
		}

		if (i%100000 == 0) // print out result in the terminal
		{
			cout <<"g= "<<g<<"  Min = "<<hisvmin<<" Mean = "<< hisvmean <<" # of Ver Rod: "<<nv <<" # of Hor Rod: "<<nh <<  "      WF[0] = "<< WF[0]<<" WF[0.2*full#] = "<< WF[0.2*0.8*V/K]<<" WF[0.5*full#] = "<< WF[0.5*0.8*V/K]<<"  WF[0.7*full#] = "<< WF[0.7*0.8*V/K]<<endl;
			// cout <<"g= "<<g<<"  Min = "<<hisvmin<<" Mean = "<< hisvmean <<" # of Ver Rod: "<<nv <<" # of Hor Rod: "<<nh <<  "      WFV[0] = "<< WFV[0]<<" WFV[50] = "<< WFV[50]<<" WFV[70] = "<< WFV[70]<<"  WFV[100] = "<< WFV[100]<<endl;
		}        	        
	}
   

	for(int i = 0; i< r*c+1; i++)
	{
		sh<<WF[i]<<endl;
	}

	ofstream myfile ("Weight_function.txt");
	myfile.precision(20);
	string data2 = sh.str();
	myfile << data2;
	myfile.close();
	histotalacc.plot(0);

	return WF;    	
}



void MC::plot(const vector<HR>& VRodlist, const vector<HR>& HRodlist)
{
	stringstream stv,sth;

	for (int i = 0; i < VRodlist.size(); i++)
	{
		stv<< VRodlist[i].getX() << "   "<< VRodlist[i].getY()<<endl;
	}

	ofstream myfilev ("2dplotv.txt");
	string datav = stv.str();
	myfilev << datav;
	myfilev.close();

	for (int j = 0; j < HRodlist.size(); j++)
	{
		sth<< HRodlist[j].getX() << "   "<< HRodlist[j].getY() <<endl;
	}

	ofstream myfileh ("2dploth.txt");
	string datah = sth.str();
	myfileh << datah;
	myfileh.close();
}



int main()
{
	double start = clock();

	// ======================= MCRUN & Plotting the final config ===============================
	array<double,10000>  wf;
	vector<HR> R;
	MC m(1E8L,8,128,128,1);
	wf = m.MCRUN();
	// ======================= end of simulation, print out the time =======
	double end = clock();
	cout <<"This simulation takes "<< (double(end-start)/CLOCKS_PER_SEC)<<endl;
	return 0;
}

