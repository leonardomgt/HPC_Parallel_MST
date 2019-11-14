
#include <math.h>
#include <string.h>
#include <stdio.h>

struct Matrix
{
	double *values;		// all Non-zero values in the matrix
	int *col_ind;		// column index for each Non-zero value
	int *row_ptr_begin; // first value index for each row
	int *row_ptr_end;   // last value index for each row

	int n_rows;	// number of rows in the matrix
	int n_cols;	// number of cols in the matrix
	int non_zeros; // number of non-zero values in the matrix
};

// Cleaner way to load the matrix from file
bool load_matrix(const char *filename, const int max_n_elements, const int max_n_rows, Matrix &matrix);

// Returns the value (row, column) of the given matrix
double get_value(const Matrix matrix, int row, int column);

// Sets the position (row, column) of the given matrix with the given value
void set_value(Matrix &matrix, int row, int column, double value);

// Searches through the given column and returns the highest absolute value
int findPivotRow(const Matrix matrix, int column);

// Swaps the row r1 and r2 in the given matrix.
void changeRows(Matrix &matrix, int r1, int r2);

// Multiplies matrix by vector, returning the result
void multiplyMatrixByVector(const Matrix matrix, const double vector[], double result[]);

// Print matrix and vector with meta data
void printMatrix(const Matrix matrix);
void printVector(double vector[], int size);
