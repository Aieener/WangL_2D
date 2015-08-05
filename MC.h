#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "square.h"
#include "cells.h"
#include "hardrods.h"
#include <cstdlib>
#include <cmath>
#include <time.h>
#include <deque>
#include <array>
using namespace std;

#ifndef MC_H
#define MC_H

class MC
{
    private:
    	//data members;
    	std::deque<HR> Rodlist; // the list storage the Rods;
    	int r,c;
    	int length;
    	long int step;
    	double z; 
    	double nh,nv,dh,dv,ah,av;
        

    public:
    	MC(long int ST, int LEN,int C, int R, double Z);

    	// ********* Getters********//
    	deque<HR> getRodlist();
    	double getTho() const;
    	double getQ() const;
    	double getAaccp() const;
    	double getDaccp() const;
        double getNh() const;
        double getNv() const;
        // ******** Setters ******//
        void setRodlist(std::deque<HR> RodL);



    	// ******** Other Functianality *******//
        int Add(Cells &s,double &prob,double &probav, double &probah);
        void Del(Cells &s,double &prob,double &probdv, double &probdh,double &size);
    	array<double,10000>  MCRUN(int o);
        void wfplot(array<double,10000> wf) const;
        void Zvs_();
        void MCRUNCHECK();

    	void plot(const deque<HR>& Rodlist);

};

#endif /* MC_H */