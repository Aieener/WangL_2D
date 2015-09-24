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
*  @return        vector <Square *>, a vector of pointers to the Squares
*                 in this cells which are neighbors of `sq`
*
*  Note:  Fill the vector in this order:  North (up one), East (right one),
*  South (down one), and West (left one).
*
*  Note:  If there is NO neighbor in a given direction, the vector will
*  have size < 4.  Corner squares have just two neighbors, and side
*  squares have three neihbors.
*/
vector <Square *> Cells::neighbors(Square sq) const
{
	vector<Square *> v1; // declare vector of pointers to the Squares.

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

			//push those square pointers to the vector;			
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

			//push those square pointers to the vector;			
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
			//push those square pointers to the vector;			
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

			//push those square pointers to the vector;			
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

			//push those square pointers to the vector;			
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

			//push those square pointers to the vector;			
			v1.push_back(up);
			v1.push_back(left);

		}
	}
	
	return v1;
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
