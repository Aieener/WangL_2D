/*
* S2LG.cpp
* Simulation of 2-D lattice gas By WangL
* Author: Yuding Ai
* Date: 2015.06.24
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
#include <deque>
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

deque<HR> MC::getRodlist() 
{
	return Rodlist;
}

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

void MC::setRodlist(std::deque<HR> RodL)
{
	Rodlist = RodL;
}

int MC::Add(Cells &s,double &prob,double &probav, double &probah)
{
	int x,y,o; // pick a random position and orientation for the HR to be added;
	x = rand()%c;
	y = rand()%r;
	o = rand()%2;

	if(s.getSquare(x,y).isEmpty()) // if it's open, try to do Addition;
	{
		HR rod(x,y,length,o);
		//======================== Vertical inside boundary===============================
		if(o == 0)
		{
			if(prob <= probav)
			{
				if(y + length <= r)
				{
					// the vertical case
					int counter = 0;

					for (int j = 0; j < length-1; j++)
					{
						// check if the vertical space is wide open
						if(s.getSquare(x,y+j+1).isOccupied())
						{
							counter++;
						}						
					}

					if (counter == 0)
					{
						// Do addition;
						// push the new rod into the Rodlist;
						Rodlist.push_front(rod);
						av++;
						nv++;// accumulate the # of ver rod;
						// update new N, E and new config;
						for (int i = 0; i < length; i++)
						{	
							s.getSquare(x,y+i).setStatus(1);
						}
						return 1;		
					}

					else{
						return 0; // forbiden move
					}						
				}
				//========================== Vertical Apply peiodic boundary ===================

				else 
				{ 
					// the vertical case apply periodic BC
					int counter2 = 0;
					for (int j = 0; j <r-y-1; j++)
					{
						// check if the vertical space is wide open
						if(s.getSquare(x,y+j+1).isOccupied())
						{
							counter2++;
						}
					}
					// if (counter2 == 0)
					// {
						for (int i = 0; i < y+length-r; i++)
						{
							// check if the vertical space is wide open
							if(s.getSquare(x,i).isOccupied())
							{								
								counter2++;
							}
						}
					// }

					if (counter2 == 0)
					{
						// Do addition;
						// push the new rod into the Rodlist;
						Rodlist.push_front(rod);	
						av++;
						nv++;// accumulate the # of ver rod;
						
						for (int j = 0; j <r-y; j++)
						{
							s.getSquare(x,y+j).setStatus(1);
						}
						for (int i = 0; i < y+length-r; i++)
						{
							s.getSquare(x,i).setStatus(1);
						}
						return 1;							
					}
					else{
						return 0;// forbiden move
					}					
				}
			}
			else{
				return 1; // the case that reject the move but not forbiden!
			}			
		}

		else 
		{
        //======================= Horizontal inside boundary ============================
			if(prob <= probah)
			{
				if( x + length <= c)
				{
					int counter3 = 0;
					for (int j = 0; j< length-1 ; j++)
					{
						// check if the horizontal space is wide open
						if(s.getSquare(x+1+j,y).isOccupied())
						{
							counter3++;
						}							
					}
					if (counter3 == 0)
					{
						//Do addition;
						//push the new rod into the Rodlist;
						Rodlist.push_back(rod);
						ah++;
						nh++;// accumulate the # of hor rod;

						// update new N, E and new config;
						for (int i = 0; i < length; i++)
						{
							s.getSquare(x+i,y).setStatus(1);
						}
						return 1;
					}
					else{
						return 0;
					}					
				}
				//======================= Horizontal periodic boundary ============================

				else
				{ 		
					// the Horizontal case apply periodic BC
					int counter4 = 0;
					for (int j = 0; j <c-x-1; j++)
					{
						// check if the Horizontal space is wide open
						if(s.getSquare(x+j+1,y).isOccupied())
						{
							counter4++;
						}
					}
					// if (counter4 == 0)
					// {
						for (int i = 0; i < x + length-c; i++)
						{
							// check if the Horizontal space is wide open
							if(s.getSquare(i,y).isOccupied())
							{
								counter4++;
							}
						}						
					// }

					if (counter4 == 0)
					{
						// Do addition;
						// push the new rod into the Rodlist;
						Rodlist.push_back(rod);	
						ah++;
						nh++;// accumulate the # of hor rod;
						
						for (int j = 0; j <c-x; j++)
						{
							s.getSquare(x+j,y).setStatus(1);
						}
						for (int i = 0; i < x+length-c; i++)
						{
							s.getSquare(i,y).setStatus(1);
						}
						return 1;							
					}

					else{
						return 0; // forbiden move
					}					
				}				    												
			}
			else{
				return 1; // the case that reject the move but not forbiden!
			}			
		}

    }
    else 
    {    	
    	return 0; // the forbiden move
    }
}


void MC::Del(Cells &s,double &prob,double &probdv, double &probdh,double &size)
{
	//Do Del;
	int DE; //pick a random config of rod to delete with 50% 50% chance for eachl;
	DE = rand()%2;

	if(DE == 0) // delete Vertical rod; which means delete indx from Rodlist[0,nv-1]
	{
		if(Rodlist[0].getOrientation()==0)// make sure there are Vertical rod;
		{
			int indx; // pick a random index from the Rodlist;
			indx = rand()%int(nv);

			//remove Rodlist[indx];
			int x,y;// the position of the target on the cells;
			x = Rodlist[indx].getX();
			y = Rodlist[indx].getY();

			if(prob <= probdv)
			{
			// --------------------- it's a vertical rod -----------------------
			// ============== the case rod is inside the Boundary ==============
				if(y + length <= r)
				{					
					for(int i = 0; i<Rodlist[indx].getLength(); i++)
					{
						// update the new config of cells
						s.getSquare(x,y + i).setStatus(0);
					}
					// remove the target rod from the deque Rodlist;
					Rodlist.erase(Rodlist.begin() + indx);
					nv--;// substract the # of ver rod;
					dv++;
				}
				
				else
				{
					// ==============the case apply periodic Boundary============
					for (int j = 0; j <r-y; j++)
					{
						s.getSquare(x,y+j).setStatus(0);
					}
					for (int i = 0; i < y+length-r; i++)
					{
						s.getSquare(x,i).setStatus(0);
					}
					Rodlist.erase(Rodlist.begin() + indx);
					nv--;// substract the # of ver rod;
					dv++;
				}
			}										
		}
	}

	else
	{
		if(Rodlist[size-1].getOrientation()==1)// make sure there are Hor rod;
		{
			int indx;
			indx = rand()%int(nh) + int(nv); // redefine indx from Rodlist[nv,nv+nh-1] 

			//remove Rodlist[indx];
			int x,y;// the position of the target on the cells;
			x = Rodlist[indx].getX();
			y = Rodlist[indx].getY();
			// --------------------- it's a Horizontal rod -----------------------
			if(prob <= probdh)
			{
                // ==============the case rod is inside the Boundary============
				if(x + length <= c)
				{					
					for(int i = 0; i<Rodlist[indx].getLength(); i++)
					{
						// update the new config of cells
						s.getSquare(x+i,y).setStatus(0);
					}
					// remove the target rod from the deque Rodlist;
					Rodlist.erase(Rodlist.begin() + indx);
					nh--;// substract the # of hor rod;
					dh++;
				}

				else
				{
					// ==============the case apply periodic Boundary============
					for (int j = 0; j <c-x; j++)
					{
						s.getSquare(x+j,y).setStatus(0);
					}
					for (int i = 0; i < x+length-c; i++)
					{

						s.getSquare(i,y).setStatus(0);
					}

					Rodlist.erase(Rodlist.begin() + indx);
					nh--;// substract the # of hor rod;
					dh++;
				}
			}
		}
	}
}


array<double,10000>  MC::MCRUN(int o)
{
	Cells s(c,r); //  setting the lattice;
	// ******************  if there is an initial state:************************** //
	// Rodlist = s.Initial(length,753,1);
	// int k = 0;
	// for(int i = 0; i < Rodlist.size();i++)
	// {
	// 	if (Rodlist[i].getOrientation() == 0)
	// 	{
	// 		k++;
	// 	}
	// }
	// nv = av = k;
	// nh = ah = Rodlist.size() - k;
	// ******************  finish setting initial state************************** //
    
    //==========================================================  declare instance variables ============================================================= //
	stringstream sh;
	sh.precision(20);
	double addordel;           // the prob to decide either add or del;
	double probah,probav;      // the acceptance prob of addition; proba = min(1.0,aaccp);
	double probdh,probdv;      // the acceptance prob of deletion; probd = min(1.0,daccp);
	double prob;               // the prob to decide either accept add/del;
	double aaccph,aaccpv;      // the acceptance probabilities: 
	double daccph,daccpv;      // the acceptance probabilities: 
	double V = double(r*c);    // the total lattice size
	double K = double(length); //
    // double WF[400] = {};
    array<double,10000> WFV;
	array<double,10000> WFH;
    double g = 1;
		
	srand(time(NULL));
	long int i = 0;
	Histogram histotalv(0,0.8*512,4); // take 80% of the full range.
	Histogram histotalh(0,0.8*512,4);

	Histogram hisv(0, 0.8*512, 4);
	Histogram hish(0, 0.8*512, 4);

	int av = 0; 
	int ah = 0; // an interger keep track of the times of we reset the histogram

	// =============================================================Start MC runs ======================================================================== //
	while (g>=1E-7)
	// while (i<step)
	{
		i++;
		// generate a random probability to decide either add or del;
		addordel = rand()%2;
		double size = nv+nh;
		prob = ((double) rand() / (RAND_MAX)); 

		if(o == 0) // Generating WFV
		{
			aaccph = (z*(V/2.0)/((nh+1.0)*K));//M.S.Shell page10/22
			aaccpv = (1*(V/2.0)/((nv+1.0)*K))*(exp(WFV[int(nv+1)] - WFV[int(nv)]));

			daccph = (((nh)*K)/(z*(V/2.0)));
			daccpv = (((nv)*K)/(1*(V/2.0)))*(exp(WFV[int(nv-1)] - WFV[int(nv)]));	
		}

		else // generating WFH
		{
			aaccph = (1*(V/2.0)/((nh+1.0)*K))*(exp(WFH[int(nh+1)] - WFH[int(nh)]));//M.S.Shell page10/22
			aaccpv = (z*(V/2.0)/((nv+1.0)*K));

			daccph = (((nh)*K)/(1*(V/2.0)))*(exp(WFH[int(nh-1)] - WFH[int(nh)]));
			daccpv = (((nv)*K)/(z*(V/2.0)));
		}

		// aaccpv = (z*(V/2.0)/((nv+1.0)*K))*(exp(WFV[int(nv+1)] - WFV[int(nv)]));
		// daccpv = (((nv)*K)/(z*(V/2.0)))*(exp(WFV[int(nv-1)] - WFV[int(nv)]));
		// aaccph = (z*(V/2.0)/((nh+1.0)*K))*(exp(WFH[int(nh+1)] - WFH[int(nh)]));//M.S.Shell page10/22
		// daccph = (((nh)*K)/(z*(V/2.0)))*(exp(WFH[int(nh-1)] - WFH[int(nh)]));

		probdh = min(1.0,daccph);
		probdv = min(1.0,daccpv);
		probah = min(1.0,aaccph);
		probav = min(1.0,aaccpv);	



        // ===========================Addition ===================================
		if(addordel == 0) 
		{
			if(nv <= 0.8*512) // make sure does not go beyond the histogram
			{
				//Do Addition;
				Add(s,prob,probav,probah);

			}
	        WFH[int(nh)] -= g;
			WFV[int(nv)] -= g;	

		}

		// ============================Deletion=============================
		else
		{
			if (size != 0) // make sure there are rods to be del;
			{
				//Do deletion;
				Del(s,prob,probdv,probdh,size);
	            WFH[int(nh)] -= g;
				WFV[int(nv)] -= g;
			}			
		}

		// ======================= Record the datas =============================================

		//================================================ load data into histogram ===================================================================//
		hisv.record(nv);
		histotalv.record(nv);

		hish.record(nh);
		histotalh.record(nh);
	
		//************************************************ check if the current histogram nv is "flat enough" *********************************************//
		double hisvmin,hisvmean;
		hisvmin = double(hisv.Minave().first);
		hisvmean = double(hisv.Minave().second);

		double hishmin,hishmean;
		hishmin = double(hish.Minave().first);
		hishmean = double(hish.Minave().second);

        if (o == 0) // we are generating the weighting function for N of vertical rod
        {
			if(hisvmin/hisvmean >=0.8)
			{   
				av++;
				g = 0.5*g;
				// his.plot(a);
				hisv.reset();
			}

			if (i%100000 == 0) // print out result in the terminal
			{
				cout <<"g= "<<g<<"  Min = "<<hisvmin<<" Mean = "<< hisvmean <<" # of Ver Rod: "<<nv <<" # of Hor Rod: "<<nh <<  "      WFV[0] = "<< WFV[0]<<" WFV[200] = "<< WFV[200]<<" WFV[300] = "<< WFV[300]<<"  WFV[400] = "<< WFV[400]<<endl;
				// cout <<"g= "<<g<<"  Min = "<<hisvmin<<" Mean = "<< hisvmean <<" # of Ver Rod: "<<nv <<" # of Hor Rod: "<<nh <<  "      WFV[0] = "<< WFV[0]<<" WFV[50] = "<< WFV[50]<<" WFV[70] = "<< WFV[70]<<"  WFV[100] = "<< WFV[100]<<endl;
			}        	
        }

        else // we are generating the weighting function for N of Hor rod
        {
			if(hishmin/hishmean >=0.8)
			{   
				ah++;
				g = 0.5*g;
				// his.plot(a);
				hish.reset();
			}

			if (i%100000 == 0) // print out result in the terminal
			{
				cout <<"g= "<<g<<"  Min = "<<hishmin<<" Mean = "<< hishmean <<" # of Ver Rod: "<<nv <<" # of Hor Rod: "<<nh <<  "      WFH[0] = "<< WFH[0]<<"WFH[150] = "<< WFH[150]<<"WFH[300] = "<< WFH[300]<<"  WFH[560] = "<< WFH[560]<<endl;
			} 
        }
	}
   
    if (o == 0) // we are generating the weighting function for N of vertical rod
    {
		for(int i = 0; i< r*c+1; i++)
		{
			sh<<WFV[i]<<endl;
		}

		ofstream myfile ("VWeight_function.dat");
		myfile.precision(20);
		string data2 = sh.str();
		myfile << data2;
		myfile.close();
		histotalv.plot(0);

		return WFV;    	
    }
    else
    {
    	for(int i = 0; i< r*c+1; i++)
		{
			sh<<WFH[i]<<endl;
		}

		ofstream myfile2 ("HWeight_function.dat");
		myfile2.precision(20);
		string data3 = sh.str();
		myfile2 << data3;
		myfile2.close();
		histotalh.plot(0);

		return WFH;  
    }    
}



void MC::plot(const deque<HR>& Rodlist)
{

	FILE* gnuplot = popen("gnuplot -persistent","w");
	// fprintf(gnuplot, "set terminal png \n  set output 'RvsNi.png'\n");
	fprintf(gnuplot, "set grid\n f(x)=0\n set xrange [0:%d]\nset yrange [0:%d]\n",r,c);
	fprintf(gnuplot, "set xtics 0,1,%d\n set ytics 0,1,%d\n set format x\"\"\n set format y\"\"\n", r-1,c-1);
	for(int i = 0; i< (nh+nv); i++)
	{
		double x = Rodlist[i].getX();
		double y = Rodlist[i].getY();
		double ori = Rodlist[i].getOrientation();

		if(ori ==1)
		{
			fprintf(gnuplot, "set object %d rect from %lf,%lf to %lf,%lf front fc rgb \"blue\" fillstyle solid 1.0\n", i + 1, x, y, x + Rodlist[i].getLength(), y + 1);
		}
		else
		{
			fprintf(gnuplot, "set object %d rect from %lf,%lf to %lf,%lf front fc rgb \"red\" fillstyle solid 1.0\n", i + 1, x, y, x + 1, y + Rodlist[i].getLength());
		}

	}
	fprintf(gnuplot, "plot f(x) ls 0 notitle\n");
	fflush(gnuplot);
	pclose(gnuplot);
}

void MC::wfplot(array<double,10000> wf) const
{
	// stringstream sh;
	// Histogram his(0, r*c, 16);
	// for(int i = 0; i<(1028/16);i++)
	// {
	// 	cout << exp(-wf[i]) << endl;
	// 	his.record(wf[i]);
	// }
	// his.plot(0);


	// for(int i = 0; i< 100; i++)
	// {
	// 	sh<< i*10 + 0.5*10<< "  "<< wf[i] << endl;
	// }

	// ofstream myfile ("mchis.dat");
	// string data2 = sh.str();
	// myfile << data2;
	// myfile.close();
	// FILE* gnuplot = popen("gnuplot -persistent","w");
	// fprintf(gnuplot, "plot \'his.dat\' using 1:2 smooth freq with boxes \n");
	// fflush(gnuplot);
	// pclose(gnuplot);
}


void MC::MCRUNCHECK()
{
	Cells s(c,r);

	// ******************  if there is an initial state:************************** //
	// Rodlist = s.Initial(length,753,1);
	// int k = 0;
	// for(int i = 0; i < Rodlist.size();i++)
	// {
	// 	if (Rodlist[i].getOrientation() == 0)
	// 	{
	// 		k++;
	// 	}
	// }
	// nv = av = k;
	// nh = ah = Rodlist.size() - k;
	// ******************  finish setting initial state************************** //

	stringstream st;

	double addordel; // the prob to decide either add or del;
	double probah,probav; // the acceptance prob of addition; proba = min(1.0,aaccp);
	double probdh,probdv; // the acceptance prob of deletion; probd = min(1.0,daccp);
	double prob; // the prob to decide either accept add/del;
	double aaccph,aaccpv; 
	double daccph,daccpv;
	double Q; // the fraction of hor and ver particle;
	double tho; // the density 
	double AD;// addition and deletion fraction
    array<double,101> WF;
	double V = double(r*c);    // the total lattice size
	double K = double(length); //
    //================= recording dat file into a sstring

    
    //================= recording WF ==============
    ifstream myfilewf ("VWeight_function.dat");
    double record;
    string number;
    for ( int w = 0; w<104; w++)
    {
    	getline(myfilewf,number);
    	record = stod(number);
    	WF[w] = record;
    	cout << WF[w]<<endl;
    }
    myfilewf.close();

		
	srand(time(NULL));
	long int i = 0;
	Histogram his(0, r*c, 1); // the histogram of nv

	//================================Start my MC simulation=================================
	while (i<step)
	// while(false)
	{
		i++;
		// generate a random probability to decide either add or del;
		addordel = rand()%2;
		double size = av+ah-dv-dh;
		
		prob = ((double) rand() / (RAND_MAX)); 
		tho = K*size/V;

		aaccph = (z*(V/2.0)/((nh+1.0)*K))*(exp(WF[int(nh+1)] - WF[int(nh)]));//M.S.Shell page10/22
		aaccpv = (z*(V/2.0)/((nv+1.0)*K))*(exp(WF[int(nv+1)] - WF[int(nv)]));// LATTICE GAS CASE

		daccph = (((nh)*K)/(z*(V/2.0)))*(exp(WF[int(nh-1)] - WF[int(nh)]));
		daccpv = (((nv)*K)/(z*(V/2.0)))*(exp(WF[int(nv-1)] - WF[int(nv)]));// LATTICE GAS CASE

		probdh = min(1.0,daccph);
		probdv = min(1.0,daccpv);
		probah = min(1.0,aaccph);
		probav = min(1.0,aaccpv);

		//******************* The sturcture of my deque list of HR ***********************
		// the Vertical rod is always push in the front
		// the Horizontal rod is always push in the back
		// the index of last vertical rod in the list can be found by index[nv-1]
		// *******************************************************************************

        // ===========================Addition ===================================
		if(addordel == 0) 
		{
			//Do Addition;
			Add(s,prob,probav,probah);
		}

		// ============================Deletion=============================
		else
		{
			if (size != 0) // make sure there are rods to be del;
			{
				//Do deletion;
				Del(s,prob,probdv, probdh,size);
			}			
		}

		// ======================= Record the datas =============================================		
        Q = (nv - nh)/(nh + nv);
		AD = (av+ah-dv-dh)/(av+ah+dv+dh);

		if (i%(step/10000) == 0)
		{
			his.record(nv);
			st << i << "         " << Q <<"        "<< nv << "          "<< nh << "         "<< tho << "         "<< AD<< "         "<< endl;
			cout <<"Process: "<< ((10000*i)/step)/100.00 <<"%"<<"    "<<"SIZE: "<<av+ah-dv-dh<<"    "<<"# of Ver Rod: "<<nv<<"    "<<"# of Hor Rod: "<< nh <<"   "<<"Qis "<<Q <<"   "<<"tho is: "<<tho << endl;
		}
	}
	// Record the data into a txt file
	ofstream myfile3 ("dataplot.dat");
	string data = st.str();
	myfile3 << data;
	myfile3.close();
	his.plot(0);
}


int main()
{
	double start = clock();

	// ======================= MCRUN & Plotting the final config ===============================
	array<double,10000>  wf;
	deque<HR> R;
	MC m(1E8L,8,64,64,14.02555);
	wf = m.MCRUN(0);
	R= m.getRodlist();
	m.plot(R);
	m.wfplot(wf);

	// ======================= check the wf to see if it's flat ================================
	// MC mc(1E8L,8,32,32,14);
	// deque<HR> RC;
	// mc.MCRUNCHECK();
	// RC= mc.getRodlist();
	// mc.plot(RC);

	// ======================= end of simulation, print out the time =======
	double end = clock();
	cout <<"This simulation takes "<< (double(end-start)/CLOCKS_PER_SEC)<<endl;
	return 0;
}

