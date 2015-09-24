// cells.h	
// 2-D lattice gas
// Author: Yuding Ai
// Date: 2015.06.05

#ifndef CELLS_H
#define CELLS_H
#include <vector>
#include <vector>
#include "hardrods.h"
#include "square.h"
using namespace std;

// Rule 1 : one particle --> one cells
// Rule 2 : for nearest neighboring, energy is lowered by epsilon (ignore for now);
class Cells
{
	private:
		Square ** arr;  // pointer to a 2-dimensional dynamic array of Squares
		int numRows;    // number of rows in the cells
		int numCols;    // number of cols in the cells



		// // disabled copy constructor (not needed, don't implement)
		// Cells(const Cells &other);

		// // disabled assignment operator (not needed, don't implement)
		// Cells& operator=(const Cells &other);
	public:
		Cells();

		/*
		*  Constructor builds a cells from an input stream
		*  @param      int r - row index
		*  @param      int c - column index
		*  initialize all states to 0;
		*/
		Cells(int r,int c);

		/*
		*  Destructor
		*/
		~Cells();

		// *** Getters *** //
		int getNumRows() const;
		int getNumCols() const;



		// *** Other Functionality *** //
		/*
		*  @function      getSquare
		*  @param         int i - row index, in range 0 through numRows - 1
		*  @param         int j - column index, in range 0 through numCols - 1
		*  @return        Square& - return a reference to the (i,j)-Square stored
		*                 in this cells
		*/
		Square& getSquare(int i, int j) const;

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

		vector <Square *> neighbors(Square sq) const;

		/*
		*  @function      toString
		*  @return        a string for printing the "text" version of the cells
		*
		*  For example, the text version might look like
		*        _ _ _ _
		*        _ # # _
		*        # _ _ _
		*/
		std::string toString() const;



	};


	/*
	*  output stream operator, for printing the "text" version of the cells
	*/
	std::ostream& operator<< (std::ostream &out, const Cells &cells);


#endif /* CELLS_H */