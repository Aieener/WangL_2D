// square.cpp
// 2-D lattice gas
// Author: Yuding Ai
// Date: 2015.06.05

#include <iostream>
#include <sstream>
#include <string>
#include "square.h"
using namespace std;
#include <cstdlib>

/*!
*  Constructor with 3 parameters
*  @param      int s - status code, should be one of:
*              0 (Unvisited), 1 (Occupied)
*  @param      int r - y index
*  @param      int c - xumn index
*
*/
Square::Square(int x, int y)
{

	status = 0; // set default status to EMPTY;
	Y = y;
	X = x;
}

// *** Setters *** //
/*
*  @function         setStatus
*  @param            int newStatus - an int in the range 0 and 1
*/
void Square::setStatus(int newStatus)
{
	status = newStatus;
}

// *** Getters *** //
int Square::getStatus() const
{
	return status;
}
int Square::getY() const
{
	return Y;
}
int Square::getX() const
{
	return X;
}

// *** Other Functionality *** //

/*
*  @function      isEmpty
*  @return        bool - return true iff this is Empty 
*/
bool Square::isEmpty() const
{
	if(status == EMPTY)
	{
		return true;
	}
	return false;
}
/*
*  @function      isOccupied
*  @return        bool - return true iff this is Occupied
*/
bool Square::isOccupied() const
{
	if(status == OCCUPIED)
	{
		return true;
	}
	return false;
}

/*
*  @function      getStatusChar
*  @return        char - return the char corresponding to the status, for
*                 printing the cell (Empty = '_', Occupied = '#'.)
*
*/
char Square::getStatusChar() const
{
	char sout;
	
	switch(status)
	{
		case 0:
		sout = '_';
		break;
		case 1:
		sout = '#';
	}
	return sout;
}

/*
*  @function      toString
*  @return        string - return a string for reporting the coordinates
*                 of the square
*
*  Example: if the square is in y  and x 3, return "[3,5]"
*/
string Square::toString() const
{
	
	stringstream rs, cs;
	rs << X;
	cs << Y;
	string ys = rs.str();
	string xs = cs.str();

	return "[" + ys + "," + xs + "]";
}


