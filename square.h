// square.h
// 2-D lattice gas
// Author: Yuding Ai
// Date: 2015.06.05

#ifndef SQUARE_H
#define SQUARE_H

#include <iostream>
#include <string>

class Square
{
	private:
		// instance variables
		int status;
		int Y;
		int X;

	public:

		// class constants
		static const int EMPTY = 0; // char for printing : _
		static const int OCCUPIED = 1;  // char for printing : # 

		/*!
		*  Default constructor
		*/   
		Square() {status = Y = X = 0;}

		/*!
		*  Constructor with 3 parameters
		*  @param      int s - status code, should be one of:
		*              0 (Unvisited), 1 (Occupied)
		*  @param      int r - y index
		*  @param      int c - x index
		*
		*/
		Square(int x, int y);

		// *** Setters *** //
		/*
		*  @function         setStatus
		*  @param            int newStatus - an int in the range 0 and 1
		*/
		void setStatus(int newStatus);


		// *** Getters *** //
		int getStatus() const;
		int getY() const;
		int getX() const;


		// ***  Other Functionality ***//

		/*
		*  @function      isEmpty
		*  @return        bool - return true iff this is Empty 
		*/
		bool isEmpty() const;

		/*
		*  @function      isOccupied
		*  @return        bool - return true iff this is Occupied
		*/
		bool isOccupied() const;

		/*
		*  @function      getStatusChar
		*  @return        char - return the char corresponding to the status, for
		*                 printing the cell (Empty = '_', Occupied = '#'.)
		*
		*/
		char getStatusChar() const;

		/*
		*  @function      toString
		*  @return        string - return a string for reporting the coordinates
		*                 of the square
		*
		*  Example: if the square is in y 5 and x 3, return "[3,5]"
		*/
		std::string toString() const;

};

#endif /* SQUARE_H */
