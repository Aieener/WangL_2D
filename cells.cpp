// cells.cpp
// 2-D lattice gas
// Author: Yuding Ai
// Date: 2015.06.05

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "cells.h"
#include "hardrods.h"
#include "MC.h"
#include <cstdlib>
#include <cmath>
using namespace std;
extern const string EXC_INVALID_DESC = "Not a valid description of cells!";

Cells::Cells()
{
	numCols = 1;
	numRows = 1;
	//Allocate memory
	arr = new Square*[numRows]; // a pointer to an array of pointers that points to square;
	for(int i = 0; i < numRows; i++)
	{
		arr [i] = new Square[numCols]; // each pointer to a square;
	}
		
	for(int j = 0 ; j<numRows*numCols; j++)
	{
		int ri, ci; // declare rowindex ri and colindex ci;
		ri = j/numCols;
	    ci = j%numCols;

		// Assign values
		Square temp = Square(ri,ci);
		arr[ri][ci] = temp;
	}
}
Cells::Cells(int x,int y)
{
	numCols = x; // x index

	numRows = y; // y index


	if (numRows < 0 || numCols < 0 || numRows*numCols == 1)
	{
		throw EXC_INVALID_DESC;
	}

	//Allocate memory
	arr = new Square*[numRows]; // a pointer to an array of pointers that points to square;
	for(int i = 0; i < numRows; i++)
	{
		arr [i] = new Square[numCols]; // each pointer to a square;
	}
		
	for(int j = 0 ; j<numRows*numCols; j++)
	{
		int ri, ci; // declare rowindex ri and colindex ci;
		ri = j/numCols;
	    ci = j%numCols;

		// Assign values
		Square temp = Square(ri,ci);
		arr[ri][ci] = temp;
	}

}

/*
*  Destructor
*/
Cells::~Cells()
{
	for (int i = 0; i < numRows; ++i)
	{
		delete [] arr[i];
	}    
    delete [] arr;
}


// *** Getters *** //

int Cells::getNumRows() const
{
	return numRows;
}
int Cells::getNumCols() const
{
	return numCols;
}

// *** Other Functionality *** //
/*
*  @function      getSquare
*  @param         int i - row index, in range 0 through numRows - 1
*  @param         int j - column index, in range 0 through numCols - 1
*  @return        Square& - return a reference to the (i,j)-Square stored
*                 in this cells
*/
Square& Cells::getSquare(int x, int y) const
{
	return arr[y][x];
}

/*
*  @function      neighbors
*  @param         Square sq
*  @return        deque <Square *>, a deque of pointers to the Squares
*                 in this cells which are neighbors of `sq`
*
*  Note:  Fill the deque in this order:  North (up one), East (right one),
*  South (down one), and West (left one).
*
*  Note:  If there is NO neighbor in a given direction, the deque will
*  have size < 4.  Corner squares have just two neighbors, and side
*  squares have three neihbors.
*/
deque <Square *> Cells::neighbors(Square sq) const
{
	deque<Square *> v1; // declare deque of pointers to the Squares.

	if (sq.getY() > 0 && sq.getY() < numRows-1)
	{
		if(sq.getX() > 0 && sq.getX() < numCols-1)
		{
		    // in this situation, the square has 4 neighbors:	

            // find the 4 Square pointers;
            Square * up = &arr[sq.getY()-1][sq.getX()];
			Square * right =  &arr[sq.getY()][sq.getX()+1];
			Square * down =  &arr[sq.getY()+1][sq.getX()];
			Square * left =  &arr[sq.getY()][sq.getX()-1];

			//push those square pointers to the deque;			
			v1.push_back(up);
			v1.push_back(right);
			v1.push_back(down);
			v1.push_back(left);			
		}

		else if(sq.getX() == 0)
		{
			// in this situation, the square has no left neighbors
			// allocate the 3 Square pointers;
			Square * up = &arr[sq.getY()-1][sq.getX()];
			Square * right = &arr[sq.getY()][sq.getX()+1];
			Square * down = &arr[sq.getY()+1][sq.getX()];

			//push those square pointers to the deque;			
			v1.push_back(up);
			v1.push_back(right);
			v1.push_back(down);
		}

		else if(sq.getX() == numCols-1)
		{
			// in this situation, the square has no right neighbors

			// allocate the 3 Square pointers;
			Square * up = &arr[sq.getY()-1][sq.getX()];
			Square * down = &arr[sq.getY()+1][sq.getX()];
			Square * left = &arr[sq.getY()][sq.getX()-1];
			//push those square pointers to the deque;			
			v1.push_back(up);
			v1.push_back(down);
			v1.push_back(left);
		}

	}

	else if (sq.getY()== 0)
	{
		if(sq.getX() > 0 && sq.getX() < numCols-1)
		{
			// in this situation, the square has no up neighbors

            // allocate the 3 Square pointers;
			Square * right = &arr[sq.getY()][sq.getX()+1];
			Square * down = &arr[sq.getY()+1][sq.getX()];
			Square * left = &arr[sq.getY()][sq.getX()-1];

			//push those square pointers to the deque;			
			v1.push_back(right);
			v1.push_back(down);
			v1.push_back(left);
		}

		else if(sq.getX() == 0)
		{
			// in this situation, the square has no up and left neighbors

			Square * right =  &arr[sq.getY()][sq.getX()+1];
			Square * down = &arr[sq.getY()+1][sq.getX()];
			v1.push_back(right);
			v1.push_back(down);
		}

		else if(sq.getX() == numCols-1)
		{
			// in this situation, the square has no up and right neighbors

			Square * down = &arr[sq.getY()+1][sq.getX()];
			Square * left = &arr[sq.getY()][sq.getX()-1];

			v1.push_back(down);
			v1.push_back(left);
		}
	}

	else if (sq.getY()== numRows-1)
	{
		if(sq.getX() > 0 && sq.getX() < numCols-1)
		{
			// in this situation, the square has no down neighbors

            // allocate the 3 Square pointers;
			Square * up = &arr[sq.getY()-1][sq.getX()];
			Square * right = &arr[sq.getY()][sq.getX()+1];
			Square * left = &arr[sq.getY()][sq.getX()-1];

			//push those square pointers to the deque;			
			v1.push_back(up);
			v1.push_back(right);
			v1.push_back(left);
		}

		else if(sq.getX() == 0)
		{
			// in this situation, the square has no down and left neighbors

			Square * up = &arr[sq.getY()-1][sq.getX()];
			Square * right = &arr[sq.getY()][sq.getX()+1]; 

			v1.push_back(up);
			v1.push_back(right);						
		}

		else if(sq.getX() == numCols-1)
		{
			// in this situation, the square has no down and right neighbors

			Square * up = &arr[sq.getY()-1][sq.getX()];
			Square * left = &arr[sq.getY()][sq.getX()-1];

			//push those square pointers to the deque;			
			v1.push_back(up);
			v1.push_back(left);

		}
	}
	
	return v1;
}

deque<HR> Cells::Initial(int len, int n, double st)
{
	// define an accumulator to indicate # of particles that sit into the cells;
	int accv = 0;
	int acch = 0;

	int c = numCols;
	int r = numRows;
	deque<HR> RodL;


	// make the random initial config;
	srand(time(0));	
	// st is chosen between [0,1]
	// when st == 0; do ver;
	// when st == 1; do hor;
	// whtn st == a; do a*ver and (1-a)*hor;
	while(accv < n*(1-st))
	{
		int x,y; // pick a random position and orientation for the HR to be added;
		x = rand()%c;
		y = rand()%r;	
		if(getSquare(x,y).isEmpty()) // if it's open, try to do Addition;
		{
			HR rod(x,y,len,0);

			//======================== Vertical inside boundary===============================
			if(y + len <= r)
			{
				// the vertical case
				int counter = 0;

				for (int j = 0; j < len-1; j++)
				{
					// check if the vertical space is wide open
					if(getSquare(x,y+j+1).isOccupied())
					{
						counter++;
					}
				}
				if (counter == 0)
				{
					// Do addition;
					// push the new rod into the Rodlist;
					RodL.push_front(rod);
					accv++;
					// update new N, E and new config;
					for (int i = 0; i < len; i++)
					{	
						getSquare(x,y+i).setStatus(1);
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
					if(getSquare(x,y+j+1).isOccupied())
					{
						counter2++;

					}
				}
				if (counter2 == 0)
				{
					for (int i = 0; i < y+len-r; i++)
					{
						// check if the vertical space is wide open
						if(getSquare(x,i).isOccupied())
						{								
							counter2++;
						}
					}
				}

				if (counter2 == 0)
				{
					// Do addition;
					// push the new rod into the Rodlist;
					RodL.push_front(rod);	
					accv++;
					
					for (int j = 0; j <r-y; j++)
					{
						getSquare(x,y+j).setStatus(1);
					}
					for (int i = 0; i < y+len-r; i++)
					{
						getSquare(x,i).setStatus(1);
					}							
				}
			}
		}
	}

	while(acch < n*st) 
	{
    //======================= Horizontal inside boundary ============================
		int x,y; // pick a random position and orientation for the HR to be added;
		x = rand()%c;
		y = rand()%r;
		if(getSquare(x,y).isEmpty())
		{
			HR rod(x,y,len,1);

			if( x + len <= c)
			{
				int counter3 = 0;
				for (int j = 0; j< len-1 ; j++)
				{
					// check if the horizontal space is wide open
					if(getSquare(x+1+j,y).isOccupied())
					{
						counter3++;
					}							
				}
				if (counter3 == 0)
				{
					//Do addition;
					//push the new rod into the Rodlist;
					RodL.push_back(rod);
					acch++;

					// update new N, E and new config;
					for (int i = 0; i < len; i++)
					{
						getSquare(x+i,y).setStatus(1);
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
					if(getSquare(x+j+1,y).isOccupied())
					{
						counter4++;

					}
				}
				if (counter4 == 0)
				{
					for (int i = 0; i < x + len-c; i++)
					{
						// check if the Horizontal space is wide open
						if(getSquare(i,y).isOccupied())
						{
							counter4++;
						}
					}						
				}

				if (counter4 == 0)
				{
					// Do addition;
					// push the new rod into the Rodlist;
					RodL.push_back(rod);	
					acch++;
					
					for (int j = 0; j <c-x; j++)
					{
						getSquare(x+j,y).setStatus(1);
					}
					for (int i = 0; i < x+len-c; i++)
					{
						getSquare(i,y).setStatus(1);
					}							
				}
			}				    												
		}
	}
	cout << acch << endl;
	cout << accv << endl;
	return RodL;
}

string Cells::toString() const
{
	stringstream in;
	string temp;
	string mg ="\n";
	for(int i = 0; i < numRows; i++)
	{
		for(int j = 0; j < numCols; j++)
		{
			char sta = getSquare(j,i).getStatusChar();
			mg = mg + sta + " ";
		}
		mg = mg + "\n";
	}

	return mg;

}


/*
*  output stream operator, for printing the "text" version of the cells
*
*  Hint:  Call toString() on the given cells and put that result into the
*  output stream.
*/
std::ostream& operator<< (std::ostream &out, const Cells &cells)
{
	out<< cells.toString();
	return out;
}
