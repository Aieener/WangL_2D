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
#include <vector>
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

vector<HR> MC::getRodlist() 
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

void MC::setRodlist(std::vector<HR> RodL)
{
    Rodlist = RodL;
}

void MC::Add(Cells &s,double &prob,double &probav, double &probah)
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
                    // cout <<"counter = "<< counter<< endl;
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
                    if (counter2 == 0)
                    {
                        for (int i = 0; i < y+length-r; i++)
                        {
                            // check if the vertical space is wide open
                            if(s.getSquare(x,i).isOccupied())
                            {                               
                                counter2++;
                            }
                        }
                    }

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
                    }
                }
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
                    if (counter4 == 0)
                    {
                        for (int i = 0; i < x + length-c; i++)
                        {
                            // check if the Horizontal space is wide open
                            if(s.getSquare(i,y).isOccupied())
                            {
                                counter4++;
                            }
                        }                       
                    }

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
                    }
                }                                                                   
            }
        }
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
                    // remove the target rod from the vector Rodlist;
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
                    // remove the target rod from the vector Rodlist;
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


void MC::MCRUN()
{
    Cells s(c,r);

    // ******************  if there is an initial state:************************** //
    // Rodlist = s.Initial(length,753,1);
    // int k = 0;
    // for(int i = 0; i < Rodlist.size();i++)
    // {
    //  if (Rodlist[i].getOrientation() == 0)
    //  {
    //      k++;
    //  }
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
    double size;
        
    srand(time(NULL));
    long int i = 0;

    //================================Start my MC simulation=================================
    while (i<step)
    {
        i++;
        // generate a random probability to decide either add or del;
        addordel = rand()%2;
        size = av+ah-dv-dh;

        // *****************define the probabilities ***********************************//
        prob = ((double) rand() / (RAND_MAX)); 
        tho = double(length*size)/double(r*c);
        aaccph = z*(double(r*c)/2.0)/(double(ah-dh+1.0)*double(length));
        aaccpv = z*(double(r*c)/2.0)/(double(av-dv+1.0)*double(length));

        daccph = (double(ah-dh)*double(length))/(z*(double(r*c)/2.0));
        daccpv = (double(av-dv)*double(length))/(z*(double(r*c)/2.0));

        probdh = min(1.0,daccph);
        probdv = min(1.0,daccpv);
        probah = min(1.0,aaccph);
        probav = min(1.0,aaccpv);

        //******************* The sturcture of my vector list of HR ***********************
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
            st << i << "         " << Q <<"        "<< nv << "          "<< nh << "         "<< tho << "         "<< AD<< "         "<< endl;
            cout <<"Process: "<< ((10000*i)/step)/100.00 <<"%"<<"    "<<"SIZE: "<<av+ah-dv-dh<<"    "<<"# of Ver Rod: "<<nv<<"    "<<"# of Hor Rod: "<< nh <<"   "<<"Qis "<<Q <<"   "<<"tho is: "<<tho << endl;
        }
    }
    // Record the data into a txt file
    ofstream myfile3 ("dataplot.dat");
    string data = st.str();
    myfile3 << data;
    myfile3.close();
}


void MC::plot(const vector<HR>& Rodlist)
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

void Zvs_()
{
    double z = 0;
    double H,V,tho,Q,miubeta,cmiubeta;  
    stringstream st;
    ofstream myfile("dataNvsZ.dat");
    for (int i =0; i < 500; i++)
    {
        MC m(1E5,1,60,60,z);
        z = double (double(10*i)/500.0);
        m.MCRUN();
        H =  m.getNh();
        V = m.getNv();
        tho = m.getTho();
        Q = (H - V)/(H + V);
        miubeta = -log(z);
        cmiubeta = (log(tho) - log(1-tho));
        st << z <<"         " << H << "             "<< V<< "             "<<tho<< "             "<<Q<< "             "<< -miubeta<< "             "<< cmiubeta << endl;
    }
    string data = st.str();
    myfile << data;
    myfile.close();
}


int main()
{
    double start = clock();

    // ======================= Plotting the final config ========================
    vector<HR> R;
    MC m(1E8L,8,120,120,35);
    m.MCRUN();
    R= m.getRodlist();
    m.plot(R);
    // ======================= get data for N vs z ========================
    // Zvs_();
    // ======================= end of simulation, print out the time =======
    double end = clock();
    cout <<"This simulation takes "<< (double(end-start)/CLOCKS_PER_SEC)<<endl;
    return 0;
}

