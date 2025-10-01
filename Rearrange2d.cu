// This file contains a number of rearrangemeent algorithms that can be possibly used
// on a 2D array. They are all deprecated except for Balance&Compress. The way the algorithm
// works is outlined in the drive, and the functions themselves have explanations in the
// header file.

#include "Rearrange2d.h"

using namespace std;

vector<int> one_d_rearrange(vector<bool> start, vector<bool> ending){
	//cout << "right before ordering config" << endl;
    vector<int> endingConfig;
    vector<int> orderedEndingSites;
    for(int i = 0;i<ending.size();i++){
        if(ending[i]){
            orderedEndingSites.push_back(i);
        }
    }
	//cout << "after ordering" << endl;
    int j = 0;
    for(int i = 0;i<start.size();i++){
        if(start[i] && j<orderedEndingSites.size()){
            endingConfig.push_back(orderedEndingSites[j]);
            j ++;
        }else{
            endingConfig.push_back(-1);
        }
    }
	//cout << "after config" << endl;
    return endingConfig;
}

void printArray(vector<vector<bool>> Array){
	cout << endl;
	for(int s1 = 0; s1<Array.size(); s1++){
		for(int s2 = 0; s2<Array[0].size(); s2++){
			cout << Array[s1][s2] << " ";
		}
	cout << endl;
	}
}

//Select a column of the array and return it as a vector
vector<bool> ColumnAt(vector<vector<bool>> Array, int dim){
    vector<bool> new_col;
    for(int k=0; k<Array.size(); k++)
        if(Array[k][dim]){
            new_col.push_back(true);
        }else{
            new_col.push_back(false);
        }
  return new_col;
}

//sum the atoms in the array by row, 0th index is the total atoms in the array
vector<int> RowSum(vector<vector<bool>> Array){
    vector<int> RowTotal(Array.size()+1);
    for(int i=0; i < Array.size(); i++){
        for(int j=0; j < Array[0].size(); j++){
            RowTotal[i+1] += int(Array[i][j]);
            RowTotal[0] += int(Array[i][j]);
        }
    }
    return RowTotal;
}

//same as rowsum but for the columns
vector<int> ColSum(vector<vector<bool>> Array){
    int index = Array[0].size() +1;
    vector<int> ColTotals(Array[0].size() + 1);
    for(int i = 0; i<Array[0].size(); i++){
        for(int j = 0; j<Array.size(); j++){
            ColTotals[i+1] += int(Array[j][i]);
            ColTotals[0] += int(Array[j][i]);
        }
    }
    return ColTotals;
}

// Sum on a row, from a particular position to the end in the direction dir
int RowSumAt(vector<vector<bool>> Array, int index, int position, bool direction){
	if(position >= Array[0].size()){
		position = Array[0].size();
	}
    int total = 0;
	if(direction){
    	for(int j=position; j < Array[0].size(); j++){
        	total += int(Array[index][j]);
		}
	}else{
	    for(int j=0; j < position +1; j++){
	        total += int(Array[index][j]);
		}
    }
    return total;
}

// Sum on a col, from a particular position to the end in the direction dir
int ColSumAt(vector<vector<bool>> Array, int index, int position, bool direction){ 
	if(position >= Array.size()){
		position = Array.size();
	}
    int total = 0;
	if(direction){
    	for(int i = position; i<Array.size(); i++){
        	total += int(Array[i][index]);
		}
    }else{
    	for(int i = 0; i < position + 1; i++){
        	total += int(Array[i][index]);
		}
	}

    return total;
}

vector<RearrangementMove> rearrangeBalanceCompress(vector<vector<bool>> &array, vector<vector<bool>> &target){
  vector<RearrangementMove> moves;
  int numrows = array.size();
  int numcols = array[0].size();

  vector<vector<int>> end;
  vector<int> temp;
  int gcounter = 0;
  int movecounter = 0;

//Balance the array to guarantee that there are enough atoms per row to meet the requirements of the target array.
	for (int i = 0; i < numrows-1; i++){
		//Try to source atoms from the above row if there are not enough atoms to make the target in this row
		int depletion = RowSumAt(target,i,0,1) > RowSumAt(array,i,0,1);
		bool flag = false;
		if (depletion > 0){
			for(int j = 0;j < numcols;j++){
				if (depletion == 0){
					break;
				}
			    if((array[i+1][j]) && !(array[i][j])){
			      	temp.push_back(1);
					array[i+1][j] = 0;
					array[i][j]= 1;
					depletion = depletion - 1;
					flag = true;
			    }else{
					temp.push_back(-1);
				}
			}
			if(flag){
				end.push_back(temp);
  				moves.push_back(RearrangementMove());
  				moves[gcounter].row = true; // Rows are movements along Lab X (AOD Y)
  				moves[gcounter].startingConfig = array[i];
  				moves[gcounter].endingConfig = end[movecounter]; //Kevin changed here to below line to force all traps to ramp up and high
  				moves[gcounter].dim = i+1;
				moves[gcounter].dir = -1;
				moves[gcounter].compression = false;
				movecounter++;
				gcounter++;
				flag = false;
			}
    		temp = {};
		}

//		Try to source atoms from two rows above if there are not enough atoms to make the target in this row
		if ((depletion > 0) && (i < numrows - 2)){
			flag = false;
			for(int j = 0;j < numcols;j++){
				if (depletion == 0){
					break;
				}
			    if((array[i+2][j]) && !(array[i+1][j]) && !(array[i][j])){
			      	temp.push_back(1);
					array[i+2][j] = 0;
					array[i][j]= 1;
					depletion = depletion - 1;
					flag = true;
			    }else{
					temp.push_back(-1);
				}
			}

			if(flag){
				end.push_back(temp);
  				moves.push_back(RearrangementMove());
  				moves[gcounter].row = true;
  				moves[gcounter].startingConfig = array[i];
  				moves[gcounter].endingConfig = end[movecounter]; //Kevin changed here to below line to force all traps to ramp up and high
  				moves[gcounter].dim = i+2;
				moves[gcounter].dir = -1;
				moves[gcounter].compression = false;
				gcounter++;

  				moves.push_back(RearrangementMove());
  				moves[gcounter].row = true;
  				moves[gcounter].startingConfig = array[i];
  				moves[gcounter].endingConfig = end[movecounter]; //Kevin changed here to below line to force all traps to ramp up and high
  				moves[gcounter].dim = i+1;
				moves[gcounter].dir = -1;
				moves[gcounter].compression = false;
				movecounter++;
				gcounter++;

				flag = false;
			}
    		temp = {};
		}


//		Try to remove atoms from this row and give them to the above row if there is space
		if (RowSumAt(array,i,0,1) > RowSumAt(target,i,0,1) ){
			int excess = RowSumAt(array,i,0,1) - RowSumAt(target,i,0,1);
			bool flag2 = false;
			//Try to source atoms from the above row if there are not enough atoms to make the target in this row
			for(int j = 0;j < numcols;j++){
				if (excess == 0){
					break;
				}
				if( (array[i][j]) && !(array[i+1][j])){
					temp.push_back(1);
					array[i+1][j] = 1;
					array[i][j]= 0;
					excess = excess - 1;
					flag2 = true;
				}else{
					temp.push_back(-1);
				}
			}
			if(flag2){
				end.push_back(temp);
				moves.push_back(RearrangementMove());
				moves[gcounter].row = true;
				moves[gcounter].startingConfig = array[i];
				moves[gcounter].endingConfig = end[movecounter]; //Kevin changed here to below line to force all traps to ramp up and high
				moves[gcounter].dim = i;
				moves[gcounter].dir = 1;
				moves[gcounter].compression = false;
				movecounter++;
				gcounter++;
				flag2 = false;
			}
			temp = {};

		}
	}

//// COMPRESSION PART OF ALGORITHM ////

for (int i = 0; i< numrows; i++){
	moves.push_back(RearrangementMove());
	moves[gcounter].row = true;
	moves[gcounter].startingConfig = array[i];
	moves[gcounter].endingConfig = one_d_rearrange(array[i],target[i]); //Kevin changed here to below line to force all traps to ramp up and high
	moves[gcounter].dim = i;
	moves[gcounter].dir = 1;
	moves[gcounter].compression = true;
	gcounter++;
  }

  return moves;
}

vector<RearrangementMove> make1DTrainOnly(vector<vector<bool>> &array, vector<int> rbRows, vector<int> csRows){
	// Takes a 2D array, tries to arrange >=1 atom/row, and moves columns to create a 1D array in the middle
	// Uses only train-type moves, so that we don't need to load compressive moves
	// printArray(array);

	// Defining variables for later
	vector<RearrangementMove> moves;
	int numRows = array.size();
	int numCols = array[0].size();
	int midCol = (numCols-1)/2; // Assuming array with odd x-dimension, this is the middle column
	int idxRbCs;
	bool isCs;
	int nextRow;
	int prevRow;
	int step;
	vector<vector<int>> end;
	int gCounter = 0; // Can maybe get rid of these counters if we instead use .back()
	int moveCounter = 0;

	// If a row is empty, try to take an atom from nearby
	for (int i = 0; i < numRows; i++) {
		if (RowSumAt(array,i,0,1) == 0) {

			// Identify whether this row is Rb/Cs, and find the neighboring rows of the same element
			idxRbCs = std::find(rbRows.begin(), rbRows.end(), i) - rbRows.begin(); // Get location of i in rbRows, if it's in there
			if (idxRbCs >= rbRows.size()) { // Index is NOT a Rb row, must be Cs
				idxRbCs = std::find(csRows.begin(), csRows.end(), i) - csRows.begin();
				isCs = true; // for row-type moves, we'll use this value to determine the species
				if (0 < idxRbCs) {
					prevRow = csRows[idxRbCs-1]; // prevRow is the next row down of the same element, in this case Cs
				} else {
					prevRow = -1; // Use -1 to indicate nonexistant row, e.g. the 0th row has no previous row
				}
				if (idxRbCs < csRows.size()-1) {
					nextRow = csRows[idxRbCs+1];
				} else {
					nextRow = -1;
				}
			} else { // Index is a Rb row
				isCs = false;
				if (0 < idxRbCs) {
					prevRow = rbRows[idxRbCs-1];
				} else {
					prevRow = -1;
				}
				if (idxRbCs < rbRows.size()-1) {
					nextRow = rbRows[idxRbCs+1];
				} else {
					nextRow = -1;
				}
			}


			// If the next-lowest same-element row has >1 atom, take from there first
			if ((prevRow != -1) && (RowSumAt(array,prevRow,0,1) > 1)) {
				// Identify which atom should be moved, update array
				vector<int> whoShouldMove(numCols, -1); // default -1 value means "do nothing"
				for (int j=0; j<numCols; j++) {
					if (array[prevRow][j]) { // prevRow has an atom, let's take it
						whoShouldMove[j] = 1;
						array[prevRow][j] = 0;
						array[i][j] = 1;
						break; // Only need to move one atom
					}
				}

				// Define the rearrangement moves to bring the atom over
				end.push_back(whoShouldMove);
				moves.push_back(RearrangementMove());
				moves[gCounter].row = true; // Rows are movements along Lab X (AOD Y)
				moves[gCounter].startingConfig = array[i];
				moves[gCounter].endingConfig = end[moveCounter];
				moves[gCounter].dim = prevRow; // dim = "which row/col are we starting with"
				moves[gCounter].dir = i-prevRow; // dir = "how far should we move" (e.g. (final row/col index)=dim+dir)
				moves[gCounter].compression = false;
				moves[gCounter].isCs = isCs;
				gCounter++;
				moveCounter++;
			}

			// Else, if the next-highest same-element row has any atom, take from there
			else if ((nextRow != -1) && (RowSumAt(array,nextRow,0,1) > 0)) {
				// Identify which atom should be moved, update array
				vector<int> whoShouldMove(numCols, -1); // default -1 value means "do nothing"
				for (int j=0; j<numCols; j++){
					if (array[nextRow][j]) { // nextRow has an atom, let's take it
						whoShouldMove[j] = 1;
						array[nextRow][j] = 0;
						array[i][j] = 1;
						break; // Only need to move one atom
					}
				}

				// Define the rearrangement moves to bring the atom over
				end.push_back(whoShouldMove);
				moves.push_back(RearrangementMove());
				moves[gCounter].row = true; // Rows are movements along Lab X (AOD Y)
				moves[gCounter].startingConfig = array[i];
				moves[gCounter].endingConfig = end[moveCounter];
				moves[gCounter].dim = nextRow;  // dim = "which row/col are we starting with"
				moves[gCounter].dir = i-nextRow; // dir = "how far should we move" (e.g. (final row/col index)=dim+dir)
				moves[gCounter].compression = false;
				moves[gCounter].isCs = isCs;
				gCounter++;
				moveCounter++;
			}
		}
	}


	// Move left-of-middle columns towards the middle to create a single filled 1D column
	// Assume columns are staggered like Cs-Rb-Cs-Rb...; each move should have a distance 2, to go from e.g. Cs->Cs
	for (int j=0; j < midCol; j++) {
		// For column moves, the Cs columns are {1,3,4,6,8}
		isCs = (j%2 == 1);

		// Move columns 2 at a time, except for the column right next to the middle
		if (j == midCol-1) {
			step = 1;
		} else {step = 2;}

		// Identify which atoms should be moved, update array
		vector<int> whoShouldMove(numRows, -1); // default -1 value means "do nothing"
		for (int i=0; i<numRows; i++){
			if (array[i][j] && !array[i][j+step]) { // Can move an atom to the right
				whoShouldMove[i] = 1;
				array[i][j] = 0;
				array[i][j+step] = 1;
			}
		}

		// Define the rearrangement moves to bring the atoms over
		end.push_back(whoShouldMove);
		moves.push_back(RearrangementMove());
		moves[gCounter].row = false; // Rows are movements along Lab X (AOD Y)
		moves[gCounter].startingConfig = ColumnAt(array,j+step); // Tells us which atoms in this column should be moved
		moves[gCounter].endingConfig = end[moveCounter];
		moves[gCounter].dim = j;
		moves[gCounter].dir = step;
		moves[gCounter].compression = false;
		moves[gCounter].isCs = isCs;
		gCounter++;
		moveCounter++;
	}

	// Move right-of-middle columns towards the middle to create a single filled 1D column
	// Assume columns are staggered like Cs-Rb-Cs-Rb...; each move should have a distance 2, to go from e.g. Cs->Cs
	for (int j=numCols-1; j > midCol; j--) {
		// For column moves, the Cs columns are {1,3,4,6,8}
		isCs = (j%2 == 0);

		// Move columns 2 at a time, except for the column right next to the middle
		if (j == midCol+1) {
			step = 1;
		} else {step = 2;}

		// Identify which atoms should be moved, update array
		vector<int> whoShouldMove(numRows, -1); // default -1 value means "do nothing"
		for (int i=0; i<numRows; i++) {
			if (array[i][j] && !array[i][j-step]) { // Can move an atom to the left
				whoShouldMove[i] = 1; //i;
				array[i][j] = 0;
				array[i][j-step] = 1;
			}
		}

		// Define the rearrangement moves to bring the atoms over
		end.push_back(whoShouldMove);
		moves.push_back(RearrangementMove());
		moves[gCounter].row = false; // Rows are movements along Lab X (AOD Y)
		moves[gCounter].startingConfig = ColumnAt(array,j-step); // Tells us which atoms in this column should be moved
		moves[gCounter].endingConfig = end[moveCounter];
		moves[gCounter].dim = j;
		moves[gCounter].dir = -step;
		moves[gCounter].compression = false;
		moves[gCounter].isCs = isCs;
		gCounter++;
		moveCounter++;
	}

	return moves;
}

vector<RearrangementMove> make1DTrainOnly_yeet(vector<vector<bool>> &array, vector<int> rbRows, vector<int> csRows, vector<bool> target){
	// Same as make1DTrainOnly, but after filling the middle column, we eject unwanted atoms

	// Defining variables for later
	vector<RearrangementMove> moves;
	int numRows = array.size();
	int numCols = array[0].size();
	int midCol = (numCols-1)/2; // Assuming array with odd x-dimension, this is the middle column
	int idxRbCs;
	bool isCs;
	int nextRow;
	int prevRow;
	int step;
	vector<vector<int>> end;
	int gCounter = 0; // Can maybe get rid of these counters if we instead use .back()
	int moveCounter = 0;

	// If a row is empty, try to take an atom from nearby
	for (int i = 0; i < numRows; i++) {
		if (RowSumAt(array,i,0,1) == 0) {

			// Identify whether this row is Rb/Cs, and find the neighboring rows of the same element
			idxRbCs = std::find(rbRows.begin(), rbRows.end(), i) - rbRows.begin(); // Get location of i in rbRows, if it's in there
			if (idxRbCs >= rbRows.size()) { // Index is NOT a Rb row, must be Cs
				idxRbCs = std::find(csRows.begin(), csRows.end(), i) - csRows.begin();
				isCs = true; // for row-type moves, we'll use this value to determine the species
				if (0 < idxRbCs) {
					prevRow = csRows[idxRbCs-1]; // prevRow is the next row down of the same element, in this case Cs
				} else {
					prevRow = -1; // Use -1 to indicate nonexistant row, e.g. the 0th row has no previous row
				}
				if (idxRbCs < csRows.size()-1) {
					nextRow = csRows[idxRbCs+1];
				} else {
					nextRow = -1;
				}
			} else { // Index is a Rb row
				isCs = false;
				if (0 < idxRbCs) {
					prevRow = rbRows[idxRbCs-1];
				} else {
					prevRow = -1;
				}
				if (idxRbCs < rbRows.size()-1) {
					nextRow = rbRows[idxRbCs+1];
				} else {
					nextRow = -1;
				}
			}

			// If the next-lowest same-element row has >1 atom, take from there first
			if ((prevRow != -1) && (RowSumAt(array,prevRow,0,1) > 1)) {
				// Identify which atom should be moved, update array
				vector<int> whoShouldMove(numCols, -1); // default -1 value means "do nothing"
				for (int j=0; j<numCols; j++) {
					if (array[prevRow][j]) { // prevRow has an atom, let's take it
						whoShouldMove[j] = 1;
						array[prevRow][j] = 0;
						array[i][j] = 1;
						break; // Only need to move one atom
					}
				}

				// Define the rearrangement moves to bring the atom over
				end.push_back(whoShouldMove);
				moves.push_back(RearrangementMove());
				moves[gCounter].row = true; // Rows are movements along Lab X (AOD Y)
				moves[gCounter].startingConfig = array[i];
				moves[gCounter].endingConfig = end[moveCounter];
				moves[gCounter].dim = prevRow; // dim = "which row/col are we starting with"
				moves[gCounter].dir = i-prevRow; // dir = "how far should we move" (e.g. (final row/col index)=dim+dir)
				moves[gCounter].compression = false;
				moves[gCounter].isCs = isCs;
				gCounter++;
				moveCounter++;
			}

			// Else, if the next-highest same-element row has any atom, take from there
			else if ((nextRow != -1) && (RowSumAt(array,nextRow,0,1) > 0)) {
				// Identify which atom should be moved, update array
				vector<int> whoShouldMove(numCols, -1); // default -1 value means "do nothing"
				for (int j=0; j<numCols; j++){
					if (array[nextRow][j]) { // nextRow has an atom, let's take it
						whoShouldMove[j] = 1;
						array[nextRow][j] = 0;
						array[i][j] = 1;
						break; // Only need to move one atom
					}
				}

				// Define the rearrangement moves to bring the atom over
				end.push_back(whoShouldMove);
				moves.push_back(RearrangementMove());
				moves[gCounter].row = true; // Rows are movements along Lab X (AOD Y)
				moves[gCounter].startingConfig = array[i];
				moves[gCounter].endingConfig = end[moveCounter];
				moves[gCounter].dim = nextRow;  // dim = "which row/col are we starting with"
				moves[gCounter].dir = i-nextRow; // dir = "how far should we move" (e.g. (final row/col index)=dim+dir)
				moves[gCounter].compression = false;
				moves[gCounter].isCs = isCs;
				gCounter++;
				moveCounter++;
			}
		}
	}


	// Move left-of-middle columns towards the middle to create a single filled 1D column
	// Assume columns are staggered like Cs-Rb-Cs-Rb...; each move should have a distance 2, to go from e.g. Cs->Cs
	for (int j=0; j < midCol; j++) {
		// For column moves, the Cs columns are {1,3,4,6,8}
		isCs = (j%2 == 1);

		// Move columns 2 at a time, except for the column right next to the middle
		if (j == midCol-1) {
			step = 1;
		} else {step = 2;}

		// Identify which atoms should be moved, update array
		vector<int> whoShouldMove(numRows, -1); // default -1 value means "do nothing"
		for (int i=0; i<numRows; i++){
			if (array[i][j] && !array[i][j+step]) { // Can move an atom to the right
				whoShouldMove[i] = 1;
				array[i][j] = 0;
				array[i][j+step] = 1;
			}
		}

		// Define the rearrangement moves to bring the atoms over
		end.push_back(whoShouldMove);
		moves.push_back(RearrangementMove());
		moves[gCounter].row = false; // Rows are movements along Lab X (AOD Y)
		moves[gCounter].startingConfig = ColumnAt(array,j+step); // Tells us which atoms in this column should be moved
		moves[gCounter].endingConfig = end[moveCounter];
		moves[gCounter].dim = j;
		moves[gCounter].dir = step;
		moves[gCounter].compression = false;
		moves[gCounter].isCs = isCs;
		gCounter++;
		moveCounter++;
	}

	// Move right-of-middle columns towards the middle to create a single filled 1D column
	// Assume columns are staggered like Cs-Rb-Cs-Rb...; each move should have a distance 2, to go from e.g. Cs->Cs
	for (int j=numCols-1; j > midCol; j--) {
		// For column moves, the Cs columns are {1,3,4,6,8}
		isCs = (j%2 == 0);

		// Move columns 2 at a time, except for the column right next to the middle
		if (j == midCol+1) {
			step = 1;
		} else {step = 2;}

		// Identify which atoms should be moved, update array
		vector<int> whoShouldMove(numRows, -1); // default -1 value means "do nothing"
		for (int i=0; i<numRows; i++) {
			if (array[i][j] && !array[i][j-step]) { // Can move an atom to the left
				whoShouldMove[i] = 1; //i;
				array[i][j] = 0;
				array[i][j-step] = 1;
			}
		}

		// Define the rearrangement moves to bring the atoms over
		end.push_back(whoShouldMove);
		moves.push_back(RearrangementMove());
		moves[gCounter].row = false; // Rows are movements along Lab X (AOD Y)
		moves[gCounter].startingConfig = ColumnAt(array,j-step); // Tells us which atoms in this column should be moved
		moves[gCounter].endingConfig = end[moveCounter];
		moves[gCounter].dim = j;
		moves[gCounter].dir = -step;
		moves[gCounter].compression = false;
		moves[gCounter].isCs = isCs;
		gCounter++;
		moveCounter++;
	}

	// Eject any atoms that are not in the target array
	for (int j=0; j<2; j++) { // Do two moves, one for Rb one for Cs
		isCs = (j%2 == 0);

		// Yeet one column over for Rb, two for Cs
		step = 2-j;

		// Identify which atoms should be moved, update array
		vector<int> whoShouldMove(numRows, -1); // default -1 value means "do nothing"
		for (int i=j; i<numRows; i+=2) { // Start at j, do every other row (i.e. all Cs or all Rb)
			if (!target[i]) { // There's an atom in the middle column, but it's not in the target
				whoShouldMove[i] = 1;
				array[i][midCol] = 0;
				array[i][midCol-step] = 1;
			} 
		}

		// Define the rearrangement moves to bring the atoms over
		end.push_back(whoShouldMove);
		moves.push_back(RearrangementMove());
		moves[gCounter].row = false; // Rows are movements along Lab X (AOD Y)
		moves[gCounter].startingConfig = target;//ColumnAt(array,midCol); // Tells us which atoms in this column should be moved
		moves[gCounter].endingConfig = end[moveCounter];
		moves[gCounter].dim = midCol;
		moves[gCounter].dir = -step;
		moves[gCounter].compression = false;
		moves[gCounter].isCs = isCs;
		gCounter++;
		moveCounter++;
	}

	return moves;
}

vector<RearrangementMove> make1DTrainOnly_yeet_split(vector<vector<bool>> &array, vector<int> rbRows, vector<int> csRows, vector<bool> target){
	// Same as make1DTrainOnly_yeet, but the column moves are split between the top and bottom half of the array

	// Defining variables for later
	vector<RearrangementMove> moves;
	int numRows = array.size();
	int numCols = array[0].size();
	int midCol = (numCols-1)/2; // Assuming array with odd x-dimension, this is the middle column
	int idxRbCs;
	bool isCs;
	int nextRow;
	int prevRow;
	int step;
	vector<vector<int>> end;
	int gCounter = 0; // Can maybe get rid of these counters if we instead use .back()
	int moveCounter = 0;

	// If a row is empty, try to take an atom from nearby
	for (int i = 0; i < numRows; i++) {
		if (RowSumAt(array,i,0,1) == 0) {

			// Identify whether this row is Rb/Cs, and find the neighboring rows of the same element
			idxRbCs = std::find(rbRows.begin(), rbRows.end(), i) - rbRows.begin(); // Get location of i in rbRows, if it's in there
			if (idxRbCs >= rbRows.size()) { // Index is NOT a Rb row, must be Cs
				idxRbCs = std::find(csRows.begin(), csRows.end(), i) - csRows.begin();
				isCs = true; // for row-type moves, we'll use this value to determine the species
				if (0 < idxRbCs) {
					prevRow = csRows[idxRbCs-1]; // prevRow is the next row down of the same element, in this case Cs
				} else {
					prevRow = -1; // Use -1 to indicate nonexistant row, e.g. the 0th row has no previous row
				}
				if (idxRbCs < csRows.size()-1) {
					nextRow = csRows[idxRbCs+1];
				} else {
					nextRow = -1;
				}
			} else { // Index is a Rb row
				isCs = false;
				if (0 < idxRbCs) {
					prevRow = rbRows[idxRbCs-1];
				} else {
					prevRow = -1;
				}
				if (idxRbCs < rbRows.size()-1) {
					nextRow = rbRows[idxRbCs+1];
				} else {
					nextRow = -1;
				}
			}

			// If the next-lowest same-element row has >1 atom, take from there first
			if ((prevRow != -1) && (RowSumAt(array,prevRow,0,1) > 1)) {
				// Identify which atom should be moved, update array
				vector<int> whoShouldMove(numCols, -1); // default -1 value means "do nothing"
				for (int j=0; j<numCols; j++) {
					if (array[prevRow][j]) { // prevRow has an atom, let's take it
						whoShouldMove[j] = 1;
						array[prevRow][j] = 0;
						array[i][j] = 1;
						break; // Only need to move one atom
					}
				}

				// Define the rearrangement moves to bring the atom over
				end.push_back(whoShouldMove);
				moves.push_back(RearrangementMove());
				moves[gCounter].row = true; // Rows are movements along Lab X (AOD Y)
				moves[gCounter].startingConfig = array[i];
				moves[gCounter].endingConfig = end[moveCounter];
				moves[gCounter].dim = prevRow; // dim = "which row/col are we starting with"
				moves[gCounter].dir = i-prevRow; // dir = "how far should we move" (e.g. (final row/col index)=dim+dir)
				moves[gCounter].compression = false;
				moves[gCounter].isCs = isCs;
				gCounter++;
				moveCounter++;
			}

			// Else, if the next-highest same-element row has any atom, take from there
			else if ((nextRow != -1) && (RowSumAt(array,nextRow,0,1) > 0)) {
				// Identify which atom should be moved, update array
				vector<int> whoShouldMove(numCols, -1); // default -1 value means "do nothing"
				for (int j=0; j<numCols; j++){
					if (array[nextRow][j]) { // nextRow has an atom, let's take it
						whoShouldMove[j] = 1;
						array[nextRow][j] = 0;
						array[i][j] = 1;
						break; // Only need to move one atom
					}
				}

				// Define the rearrangement moves to bring the atom over
				end.push_back(whoShouldMove);
				moves.push_back(RearrangementMove());
				moves[gCounter].row = true; // Rows are movements along Lab X (AOD Y)
				moves[gCounter].startingConfig = array[i];
				moves[gCounter].endingConfig = end[moveCounter];
				moves[gCounter].dim = nextRow;  // dim = "which row/col are we starting with"
				moves[gCounter].dir = i-nextRow; // dir = "how far should we move" (e.g. (final row/col index)=dim+dir)
				moves[gCounter].compression = false;
				moves[gCounter].isCs = isCs;
				gCounter++;
				moveCounter++;
			}
		}
	}

	// Move right-of-middle columns towards the middle to create a single filled 1D column
	// Assume columns are staggered like Cs-Rb-Cs-Rb...; each move should have a distance 2, to go from e.g. Cs->Cs
	for (int j=numCols-1; j > midCol; j--) {
		// For column moves, the Cs columns are {1,3,4,6,8}
		isCs = (j%2 == 0);

		// Move columns 2 at a time, except for the column right next to the middle
		if (j == midCol+1) {
			step = 1;
		} else {step = 2;}

		// Identify which atoms should be moved, update array
		vector<int> whoShouldMove(numRows, -1); // default -1 value means "do nothing"
		for (int i=0; i<numRows/2; i++) {
			if (array[i][j] && !array[i][j-step]) { // Can move an atom to the left
				whoShouldMove[i] = 1; //i;
				array[i][j] = 0;
				array[i][j-step] = 1;
			}
		}

		// Define the rearrangement moves to bring the atoms over
		end.push_back(whoShouldMove);
		moves.push_back(RearrangementMove());
		moves[gCounter].row = false; // Rows are movements along Lab X (AOD Y)
		moves[gCounter].startingConfig = ColumnAt(array,j-step); // Tells us which atoms in this column should be moved
		moves[gCounter].endingConfig = end[moveCounter];
		moves[gCounter].dim = j;
		moves[gCounter].dir = -step;
		moves[gCounter].compression = false;
		moves[gCounter].isCs = isCs;
		gCounter++;
		moveCounter++;

		// Identify which atoms should be moved, update array
		vector<int> whoShouldMove2(numRows, -1); // default -1 value means "do nothing"
		for (int i=numRows/2; i<numRows; i++) {
			if (array[i][j] && !array[i][j-step]) { // Can move an atom to the left
				whoShouldMove2[i] = 1; //i;
				array[i][j] = 0;
				array[i][j-step] = 1;
			}
		}

		// Define the rearrangement moves to bring the atoms over
		end.push_back(whoShouldMove2);
		moves.push_back(RearrangementMove());
		moves[gCounter].row = false; // Rows are movements along Lab X (AOD Y)
		moves[gCounter].startingConfig = ColumnAt(array,j-step); // Tells us which atoms in this column should be moved
		moves[gCounter].endingConfig = end[moveCounter];
		moves[gCounter].dim = j;
		moves[gCounter].dir = -step;
		moves[gCounter].compression = false;
		moves[gCounter].isCs = isCs;
		gCounter++;
		moveCounter++;
	}

	// Move left-of-middle columns towards the middle to create a single filled 1D column
	// Assume columns are staggered like Cs-Rb-Cs-Rb...; each move should have a distance 2, to go from e.g. Cs->Cs
	for (int j=0; j < midCol; j++) {
		// For column moves, the Cs columns are {1,3,4,6,8}
		isCs = (j%2 == 1);

		// Move columns 2 at a time, except for the column right next to the middle
		if (j == midCol-1) {
			step = 1;
		} else {step = 2;}

		// Identify which atoms should be moved, update array
		vector<int> whoShouldMove(numRows, -1); // default -1 value means "do nothing"
		for (int i=0; i<numRows/2; i++){
			if (array[i][j] && !array[i][j+step]) { // Can move an atom to the right
				whoShouldMove[i] = 1;
				array[i][j] = 0;
				array[i][j+step] = 1;
			}
		}

		// Define the rearrangement moves to bring the atoms over
		end.push_back(whoShouldMove);
		moves.push_back(RearrangementMove());
		moves[gCounter].row = false; // Rows are movements along Lab X (AOD Y)
		moves[gCounter].startingConfig = ColumnAt(array,j+step); // Tells us which atoms in this column should be moved
		moves[gCounter].endingConfig = end[moveCounter];
		moves[gCounter].dim = j;
		moves[gCounter].dir = step;
		moves[gCounter].compression = false;
		moves[gCounter].isCs = isCs;
		gCounter++;
		moveCounter++;

		// Identify which atoms should be moved, update array
		vector<int> whoShouldMove2(numRows, -1); // default -1 value means "do nothing"
		for (int i=numRows/2; i<numRows; i++){
			if (array[i][j] && !array[i][j+step]) { // Can move an atom to the right
				whoShouldMove2[i] = 1;
				array[i][j] = 0;
				array[i][j+step] = 1;
			}
		}

		// Define the rearrangement moves to bring the atoms over
		end.push_back(whoShouldMove2);
		moves.push_back(RearrangementMove());
		moves[gCounter].row = false; // Rows are movements along Lab X (AOD Y)
		moves[gCounter].startingConfig = ColumnAt(array,j+step); // Tells us which atoms in this column should be moved
		moves[gCounter].endingConfig = end[moveCounter];
		moves[gCounter].dim = j;
		moves[gCounter].dir = step;
		moves[gCounter].compression = false;
		moves[gCounter].isCs = isCs;
		gCounter++;
		moveCounter++;
	}

	// Eject any atoms that are not in the target array
	for (int j=0; j<2; j++) { // Do two moves, one for Rb one for Cs
		isCs = (j%2 == 0);

		// Yeet one column over for Rb, two for Cs
		step = 2-j;

		// Identify which atoms should be moved, update array
		vector<int> whoShouldMove(numRows, -1); // default -1 value means "do nothing"
		for (int i=j; i<numRows/2; i+=2) { // Start at j, do every other row (i.e. all Cs or all Rb)
			if (!target[i]) { // There's an atom in the middle column, but it's not in the target
				whoShouldMove[i] = 1;
				array[i][midCol] = 0;
				array[i][midCol-step] = 1;
			} 
		}

		// Define the rearrangement moves to bring the atoms over
		end.push_back(whoShouldMove);
		moves.push_back(RearrangementMove());
		moves[gCounter].row = false; // Rows are movements along Lab X (AOD Y)
		moves[gCounter].startingConfig = target;//ColumnAt(array,midCol); // Tells us which atoms in this column should be moved
		moves[gCounter].endingConfig = end[moveCounter];
		moves[gCounter].dim = midCol;
		moves[gCounter].dir = -step;
		moves[gCounter].compression = false;
		moves[gCounter].isCs = isCs;
		gCounter++;
		moveCounter++;

		// Identify which atoms should be moved, update array
		vector<int> whoShouldMove2(numRows, -1); // default -1 value means "do nothing"
		for (int i=j+numRows/2; i<numRows; i+=2) { // Start at j, do every other row (i.e. all Cs or all Rb)
			if (!target[i]) { // There's an atom in the middle column, but it's not in the target
				whoShouldMove2[i] = 1;
				array[i][midCol] = 0;
				array[i][midCol-step] = 1;
			} 
		}

		// Define the rearrangement moves to bring the atoms over
		end.push_back(whoShouldMove2);
		moves.push_back(RearrangementMove());
		moves[gCounter].row = false; // Rows are movements along Lab X (AOD Y)
		moves[gCounter].startingConfig = target;//ColumnAt(array,midCol); // Tells us which atoms in this column should be moved
		moves[gCounter].endingConfig = end[moveCounter];
		moves[gCounter].dim = midCol;
		moves[gCounter].dir = -step;
		moves[gCounter].compression = false;
		moves[gCounter].isCs = isCs;
		gCounter++;
		moveCounter++;
	}

	return moves;
}

vector<RearrangementMove> make1DTrainOnly_yeet_split_csOnly(vector<vector<bool>> &array, vector<int> rbRows, vector<int> csRows, vector<bool> target){
	// For test purposes

	// Defining variables for later
	vector<RearrangementMove> moves;
	int numRows = array.size();
	int numCols = array[0].size();
	int midCol = (numCols-1)/2; // Assuming array with odd x-dimension, this is the middle column
	int idxRbCs;
	bool isCs;
	int nextRow;
	int prevRow;
	int step;
	vector<vector<int>> end;
	int gCounter = 0; // Can maybe get rid of these counters if we instead use .back()
	int moveCounter = 0;

	// If a row is empty, try to take an atom from nearby
	for (int i = 0; i < numRows; i++) {
		if (RowSumAt(array,i,0,1) == 0) {

			// Identify whether this row is Rb/Cs, and find the neighboring rows of the same element
			idxRbCs = std::find(rbRows.begin(), rbRows.end(), i) - rbRows.begin(); // Get location of i in rbRows, if it's in there
			if (idxRbCs >= rbRows.size()) { // Index is NOT a Rb row, must be Cs
				idxRbCs = std::find(csRows.begin(), csRows.end(), i) - csRows.begin();
				isCs = true; // for row-type moves, we'll use this value to determine the species
				if (0 < idxRbCs) {
					prevRow = csRows[idxRbCs-1]; // prevRow is the next row down of the same element, in this case Cs
				} else {
					prevRow = -1; // Use -1 to indicate nonexistant row, e.g. the 0th row has no previous row
				}
				if (idxRbCs < csRows.size()-1) {
					nextRow = csRows[idxRbCs+1];
				} else {
					nextRow = -1;
				}
			} else { // Index is a Rb row
				isCs = false;
				if (0 < idxRbCs) {
					prevRow = rbRows[idxRbCs-1];
				} else {
					prevRow = -1;
				}
				if (idxRbCs < rbRows.size()-1) {
					nextRow = rbRows[idxRbCs+1];
				} else {
					nextRow = -1;
				}
			}

			// If the next-lowest same-element row has >1 atom, take from there first
			if ((prevRow != -1) && (RowSumAt(array,prevRow,0,1) > 1) && isCs) {
				// Identify which atom should be moved, update array
				vector<int> whoShouldMove(numCols, -1); // default -1 value means "do nothing"
				for (int j=0; j<numCols; j++) {
					if (array[prevRow][j]) { // prevRow has an atom, let's take it
						whoShouldMove[j] = 1;
						array[prevRow][j] = 0;
						array[i][j] = 1;
						break; // Only need to move one atom
					}
				}

				// Define the rearrangement moves to bring the atom over
				end.push_back(whoShouldMove);
				moves.push_back(RearrangementMove());
				moves[gCounter].row = true; // Rows are movements along Lab X (AOD Y)
				moves[gCounter].startingConfig = array[i];
				moves[gCounter].endingConfig = end[moveCounter];
				moves[gCounter].dim = prevRow; // dim = "which row/col are we starting with"
				moves[gCounter].dir = i-prevRow; // dir = "how far should we move" (e.g. (final row/col index)=dim+dir)
				moves[gCounter].compression = false;
				moves[gCounter].isCs = isCs;
				gCounter++;
				moveCounter++;
			}

			// Else, if the next-highest same-element row has any atom, take from there
			else if ((nextRow != -1) && (RowSumAt(array,nextRow,0,1) > 0) && isCs) {
				// Identify which atom should be moved, update array
				vector<int> whoShouldMove(numCols, -1); // default -1 value means "do nothing"
				for (int j=0; j<numCols; j++){
					if (array[nextRow][j]) { // nextRow has an atom, let's take it
						whoShouldMove[j] = 1;
						array[nextRow][j] = 0;
						array[i][j] = 1;
						break; // Only need to move one atom
					}
				}

				// Define the rearrangement moves to bring the atom over
				end.push_back(whoShouldMove);
				moves.push_back(RearrangementMove());
				moves[gCounter].row = true; // Rows are movements along Lab X (AOD Y)
				moves[gCounter].startingConfig = array[i];
				moves[gCounter].endingConfig = end[moveCounter];
				moves[gCounter].dim = nextRow;  // dim = "which row/col are we starting with"
				moves[gCounter].dir = i-nextRow; // dir = "how far should we move" (e.g. (final row/col index)=dim+dir)
				moves[gCounter].compression = false;
				moves[gCounter].isCs = isCs;
				gCounter++;
				moveCounter++;
			}
		}
	}

	// Move right-of-middle columns towards the middle to create a single filled 1D column
	// Assume columns are staggered like Cs-Rb-Cs-Rb...; each move should have a distance 2, to go from e.g. Cs->Cs
	for (int j=numCols-1; j > midCol; j--) {
		// For column moves, the Cs columns are {1,3,4,6,8}
		isCs = (j%2 == 0);

		if (isCs) {
			// Move columns 2 at a time, except for the column right next to the middle
			if (j == midCol+1) {
				step = 1;
			} else {step = 2;}

			// Identify which atoms should be moved, update array
			vector<int> whoShouldMove(numRows, -1); // default -1 value means "do nothing"
			for (int i=0; i<numRows/2; i++) {
				if (array[i][j] && !array[i][j-step]) { // Can move an atom to the left
					whoShouldMove[i] = 1; //i;
					array[i][j] = 0;
					array[i][j-step] = 1;
				}
			}

			// Define the rearrangement moves to bring the atoms over
			end.push_back(whoShouldMove);
			moves.push_back(RearrangementMove());
			moves[gCounter].row = false; // Rows are movements along Lab X (AOD Y)
			moves[gCounter].startingConfig = ColumnAt(array,j-step); // Tells us which atoms in this column should be moved
			moves[gCounter].endingConfig = end[moveCounter];
			moves[gCounter].dim = j;
			moves[gCounter].dir = -step;
			moves[gCounter].compression = false;
			moves[gCounter].isCs = isCs;
			gCounter++;
			moveCounter++;

			// Identify which atoms should be moved, update array
			vector<int> whoShouldMove2(numRows, -1); // default -1 value means "do nothing"
			for (int i=numRows/2; i<numRows; i++) {
				if (array[i][j] && !array[i][j-step]) { // Can move an atom to the left
					whoShouldMove2[i] = 1; //i;
					array[i][j] = 0;
					array[i][j-step] = 1;
				}
			}

			// Define the rearrangement moves to bring the atoms over
			end.push_back(whoShouldMove2);
			moves.push_back(RearrangementMove());
			moves[gCounter].row = false; // Rows are movements along Lab X (AOD Y)
			moves[gCounter].startingConfig = ColumnAt(array,j-step); // Tells us which atoms in this column should be moved
			moves[gCounter].endingConfig = end[moveCounter];
			moves[gCounter].dim = j;
			moves[gCounter].dir = -step;
			moves[gCounter].compression = false;
			moves[gCounter].isCs = isCs;
			gCounter++;
			moveCounter++;
		}
	}

	// Move left-of-middle columns towards the middle to create a single filled 1D column
	// Assume columns are staggered like Cs-Rb-Cs-Rb...; each move should have a distance 2, to go from e.g. Cs->Cs
	for (int j=0; j < midCol; j++) {
		// For column moves, the Cs columns are {1,3,4,6,8}
		isCs = (j%2 == 1);


		if (isCs) {
			// Move columns 2 at a time, except for the column right next to the middle
			if (j == midCol-1) {
				step = 1;
			} else {step = 2;}

			// Identify which atoms should be moved, update array
			vector<int> whoShouldMove(numRows, -1); // default -1 value means "do nothing"
			for (int i=0; i<numRows/2; i++){
				if (array[i][j] && !array[i][j+step]) { // Can move an atom to the right
					whoShouldMove[i] = 1;
					array[i][j] = 0;
					array[i][j+step] = 1;
				}
			}

			// Define the rearrangement moves to bring the atoms over
			end.push_back(whoShouldMove);
			moves.push_back(RearrangementMove());
			moves[gCounter].row = false; // Rows are movements along Lab X (AOD Y)
			moves[gCounter].startingConfig = ColumnAt(array,j+step); // Tells us which atoms in this column should be moved
			moves[gCounter].endingConfig = end[moveCounter];
			moves[gCounter].dim = j;
			moves[gCounter].dir = step;
			moves[gCounter].compression = false;
			moves[gCounter].isCs = isCs;
			gCounter++;
			moveCounter++;

			// Identify which atoms should be moved, update array
			vector<int> whoShouldMove2(numRows, -1); // default -1 value means "do nothing"
			for (int i=numRows/2; i<numRows; i++){
				if (array[i][j] && !array[i][j+step]) { // Can move an atom to the right
					whoShouldMove2[i] = 1;
					array[i][j] = 0;
					array[i][j+step] = 1;
				}
			}

			// Define the rearrangement moves to bring the atoms over
			end.push_back(whoShouldMove2);
			moves.push_back(RearrangementMove());
			moves[gCounter].row = false; // Rows are movements along Lab X (AOD Y)
			moves[gCounter].startingConfig = ColumnAt(array,j+step); // Tells us which atoms in this column should be moved
			moves[gCounter].endingConfig = end[moveCounter];
			moves[gCounter].dim = j;
			moves[gCounter].dir = step;
			moves[gCounter].compression = false;
			moves[gCounter].isCs = isCs;
			gCounter++;
			moveCounter++;
		}
	}

	// Eject any atoms that are not in the target array
	for (int j=0; j<2; j++) { // Do two moves, one for Rb one for Cs
		isCs = (j%2 == 0);

		if (isCs) {
			// Yeet one column over for Rb, two for Cs
			step = 2-j;

			// Identify which atoms should be moved, update array
			vector<int> whoShouldMove(numRows, -1); // default -1 value means "do nothing"
			for (int i=j; i<numRows/2; i+=2) { // Start at j, do every other row (i.e. all Cs or all Rb)
				if (!target[i]) { // There's an atom in the middle column, but it's not in the target
					whoShouldMove[i] = 1;
					array[i][midCol] = 0;
					array[i][midCol-step] = 1;
				} 
			}

			// Define the rearrangement moves to bring the atoms over
			end.push_back(whoShouldMove);
			moves.push_back(RearrangementMove());
			moves[gCounter].row = false; // Rows are movements along Lab X (AOD Y)
			moves[gCounter].startingConfig = target;//ColumnAt(array,midCol); // Tells us which atoms in this column should be moved
			moves[gCounter].endingConfig = end[moveCounter];
			moves[gCounter].dim = midCol;
			moves[gCounter].dir = -step;
			moves[gCounter].compression = false;
			moves[gCounter].isCs = isCs;
			gCounter++;
			moveCounter++;

			// Identify which atoms should be moved, update array
			vector<int> whoShouldMove2(numRows, -1); // default -1 value means "do nothing"
			for (int i=j+numRows/2; i<numRows; i+=2) { // Start at j, do every other row (i.e. all Cs or all Rb)
				if (!target[i]) { // There's an atom in the middle column, but it's not in the target
					whoShouldMove2[i] = 1;
					array[i][midCol] = 0;
					array[i][midCol-step] = 1;
				} 
			}

			// Define the rearrangement moves to bring the atoms over
			end.push_back(whoShouldMove2);
			moves.push_back(RearrangementMove());
			moves[gCounter].row = false; // Rows are movements along Lab X (AOD Y)
			moves[gCounter].startingConfig = target;//ColumnAt(array,midCol); // Tells us which atoms in this column should be moved
			moves[gCounter].endingConfig = end[moveCounter];
			moves[gCounter].dim = midCol;
			moves[gCounter].dir = -step;
			moves[gCounter].compression = false;
			moves[gCounter].isCs = isCs;
			gCounter++;
			moveCounter++;
		}
	}

	return moves;
}

vector<RearrangementMove> rearrangeCleanUpY(vector<vector<bool>> &array, vector<vector<bool>> &target){
	//printArray(array);
	vector<RearrangementMove> moves;
	int numrows = array.size();
	int numcols = array[0].size();

	//int numAtoms = 0;
	vector<vector<int>> end;
	vector<int> temp;
	int gcounter = 0;
	//int movecounter = 0;

	for (int i = 0; i< numcols; i++){
		moves.push_back(RearrangementMove());
		moves[gcounter].row = false;
		moves[gcounter].startingConfig = ColumnAt(array,i);
		moves[gcounter].endingConfig = one_d_rearrange(ColumnAt(array,i),ColumnAt(target,i)); //Kevin changed here to below line to force all traps to ramp up and high
		moves[gcounter].dim = i;
		moves[gcounter].dir = 1;
		moves[gcounter].compression = true;
		moves[gcounter].isCs = (i%2 == 1); // assuming Cs columns are 1,3,5,7,9...
		gcounter++;
	}

	return moves;
}

vector<RearrangementMove> rearrangeCleanUpX(vector<vector<bool>> &array, vector<vector<bool>> &target){
  vector<RearrangementMove> moves;
  int numrows = array.size();
  int numcols = array[0].size();

  //int numAtoms = 0;
  vector<vector<int>> end;
  vector<int> temp;
  int gcounter = 0;
  //int movecounter = 0;


for (int i = 0; i< numrows; i++){
	moves.push_back(RearrangementMove());
	moves[gcounter].row = true;
	moves[gcounter].startingConfig = array[i];
	moves[gcounter].endingConfig = one_d_rearrange(array[i],target[i]); //Kevin changed here to below line to force all traps to ramp up and high
	moves[gcounter].dim = i;
	moves[gcounter].dir = 1;
	moves[gcounter].compression = true;
	moves[gcounter].isCs = (i%2 == 1); // assuming Cs rows are 1,3,5,7,9...
	gcounter++;
  }



  return moves;
}
