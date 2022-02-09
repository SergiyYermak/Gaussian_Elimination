/*	This algorithm transforms an n*n(+1) augmented matrix M representing a linear system into its reduced form.
	At each step, M refers to the current state of the matrix, not the original state.

	A. Set the row i equal to 0.
	B. Set the column j equal to 0. We will loop though columns 0 to n-1.
	C. Find the row k with k >= i for which M_kj has the largest absolute value.
		If no such row exists for which M_kj != 0, then skip to step H.
	D. If k != i, then exchange rows k and i using elementry row operation (a) under Def. 3.3
	E. Multiply row i by 1/M_ij. This sets the (i, j) entry of M to 1 using elementary row operation (b).
	F. For each row r, where 1 <= r <= n and r != i, add -Mrj times row i to row r.
		This step clears each entry above and below row i in column j to 0 using elementary row operation (c).
	G. Increment i.
	F. If j<n, increment j and loop to step C.

	Definition 3.3 Elementary row operation
	(a) Exchange two rows.
	(b) Multiply a row by a nonzero scalar.
	(c) Add a multiple of one row to another row.

*/
#include <iostream>
#include <fstream> //read and write file
#include <iomanip> //format output to the screen
#include <cmath> //for absolute value functions

int main()
{
	//---------- Set up ----------//

	//Variables
	int rows, columns = 0;
	double *b; //RHS
	double **A; //LHS
	std::ifstream fin;
	std::ofstream fout;

	//Open file containing the matrix
	fin.open("TestMatrix.txt");

	//Find dimensions
	fin >> rows >> columns;
	//std::cout << "DEBUG: rows = " << rows << " columns = " << columns << std::endl;

	//Finish constructing matrix
	b = new double[rows]; //RHS

	A = new double *[rows]; //LHS
	for (int i = 0; i < rows; i++)
		A[i] = new double[columns];

	//Fill LHS
	for (int row = 0; row < rows; row++)
	{
		for (int col = 0; col < columns; col++)
		{
			fin >> A[row][col];
		}
	}
	//Fil RHS
	for (int row = 0; row < rows; row++)
	{
		fin >> b[row];
	}
	fin.close();

	//---------- Gaussian Elimination ----------//
	//A. Set the row i equal to 0.
	int i = 0;
	//B. Set the column j equal to 0. We will loop though columns 0 to n-1.
	int j = 0;
	//F.
	while(j < columns)
    {

        //C. Find the row k with k >= i for which M_kj has the largest absolute value.
        //   If no such row exists for which M_kj != 0, then skip to step H.
        int k = i;
        double largestAbsValue = 0.0;
        for(int i = 0; i < rows; i++)
        {
            if(abs(A[i][j]) > largestAbsValue)
                k = i;
        }

        //D. If k != i, then exchange rows k and i using elementry row operation (a) under Def. 3.3
        if(k != i)
        {
            double *temp = new double[rows];
            //RHS
            temp[0] = b[i];
            b[i] = b[k];
            b[k] = temp[0];
            for(int col = 0; col < columns; col++)
            {
                //LHS
                temp[col] = A[i][col];
                A[i][col] = A[k][col];
                A[k][col] = temp[col];
            }
        }

        //E. Multiply row i by 1/M_ij. This sets the (i, j) entry of M to 1 using elementary row operation (b).
        double M_ij = 1/A[i][j];
        //RHS
        b[i] *= M_ij;
        for(int col = 0; col < columns; col++)
        {
            //LHS
            A[i][col] *= M_ij;
        }

        //F. For each row r, where 0 <= r <= n-1 and r != i, add -Mrj times row i to row r.
		//   This step clears each entry above and below row i in column j to 0 using elementary row operation (c).
		for(int r = 0; r < rows; r++)
        {
            if(r != i)
            {
                double M_rj = A[r][j];
                //RHS
                b[r] += (-M_rj * b[i]);
                for(int col = 0; col < columns; col++)
                {
                    //LHS
                    A[r][col] = A[r][col] + (-M_rj * A[i][col]);

                }
            }
        }

        //G. Increment i
        i++;
        //F. If j<n, increment j and loop to step C. (done in the while loop)
        j++;
    }

	//---------- Display ----------//
	for (int row = 0; row < rows; row++)
	{
		for (int col = 0; col < columns; col++)
		{
			std::cout << A[row][col] << " ";
		}
		std::cout << " | " << b[row] << std::endl;
	}

	//---------- Write to file ----------//
	fout.open("solution.txt");
	for (int row = 0; row < rows; row++)
	{
		for (int col = 0; col < columns; col++)
		{
			fout << A[row][col] << " ";
		}
		fout << " | " << b[row] << std::endl;
	}

	return 0;
}

