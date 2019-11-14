
#include "sparse.h"
#include "matrix.h"

bool load_matrix(const char *filename, const int max_n_elements, const int max_n_rows, Matrix &matrix)
{
	return load_matrix_market(filename, max_n_elements, max_n_rows,
							  matrix.non_zeros, matrix.n_rows, matrix.n_cols,
							  matrix.values, matrix.col_ind, matrix.row_ptr_begin, matrix.row_ptr_end);
}

double get_value(Matrix matrix, int row, int column)
{
	for (int idx = matrix.row_ptr_begin[row]; idx <= matrix.row_ptr_end[row]; ++idx)
	{
		if (matrix.col_ind[idx] == column)
			return matrix.values[idx];
	}

	return 0;
}

void set_value(Matrix &matrix, int row, int column, double value)
{
	for (int idx = matrix.row_ptr_begin[row]; idx <= matrix.row_ptr_end[row]; ++idx)
	{
		// If the OLD value is Non-zero
		if (matrix.col_ind[idx] == column)
		{
			// If the NEW value is Non-zero
			if (value != 0)
			{
				matrix.values[idx] = value;
				return;
			}
			// If the NEW value is ZERO - Remove value from sparse matrix
			else
			{
				for (int i = idx; i < matrix.non_zeros - 1; i++)
				{
					matrix.values[i] = matrix.values[i + 1];
					matrix.col_ind[i] = matrix.col_ind[i + 1];
				}
				matrix.non_zeros--;

				matrix.row_ptr_end[row]--;
				for (int i = row + 1; i < matrix.n_rows; i++)
				{
					matrix.row_ptr_begin[i]--;
					matrix.row_ptr_end[i]--;
				}
				return;
			}
		}
	}

	// OLD value is ZERO - Didn't find the value

	// If NEW value is ZERO - Nothing to be done, return.
	if (value == 0)
		return;

	// If NEW value is Non-zero - Add value to sparse matrix
	for (int idx = matrix.row_ptr_begin[row]; idx <= matrix.row_ptr_end[row]; idx++)
	{
		if (matrix.col_ind[idx] >= column)
		{

			for (int i = matrix.non_zeros - 1; i >= idx; i--)
			{
				matrix.values[i + 1] = matrix.values[i];
				matrix.col_ind[i + 1] = matrix.col_ind[i];
			}
			matrix.values[idx] = value;
			matrix.col_ind[idx] = column;

			matrix.non_zeros++;

			matrix.row_ptr_end[row]++;
			for (int i = row + 1; i < matrix.n_rows; i++)
			{
				matrix.row_ptr_begin[i]++;
				matrix.row_ptr_end[i]++;
			}
			return;
		}
	}

	// Case which new value is the last Non-zero of the row.
	int idx = matrix.row_ptr_begin[row + 1];

	for (int i = matrix.non_zeros - 1; i >= idx; i--)
	{
		matrix.values[i + 1] = matrix.values[i];
		matrix.col_ind[i + 1] = matrix.col_ind[i];
	}
	matrix.values[idx] = value;
	matrix.col_ind[idx] = column;

	matrix.non_zeros++;

	matrix.row_ptr_end[row]++;
	for (int i = row + 1; i < matrix.n_rows; i++)
	{
		matrix.row_ptr_begin[i]++;
		matrix.row_ptr_end[i]++;
	}
	return;
}

int findPivotRow(const Matrix matrix, int column)
{
	double max = get_value(matrix, column, column);

	int rowMax = column;

	for (int j = column + 1; j < matrix.n_rows; j++) // Iterate through rows
	{
		double value = get_value(matrix, j, column);

		if (fabs(value) > fabs(max))
		{
			max = value;
			rowMax = j;
		}
	}

	return rowMax;
}

void changeRows(Matrix &matrix, int r1, int r2)
{
	if (r1 > r2)
	{
		int rt = r1;
		r1 = r2;
		r2 = rt;
	}
	else if (r1 == r2)
	{
		return;
	}

	int idx_begin_1 = matrix.row_ptr_begin[r1];
	int idx_begin_2 = matrix.row_ptr_begin[r2];
	int idx_end_1 = matrix.row_ptr_end[r1];
	int idx_end_2 = matrix.row_ptr_end[r2];

	int n_values_r1 = idx_end_1 - idx_begin_1 + 1;
	int n_values_r2 = idx_end_2 - idx_begin_2 + 1;
	int n_values_between = idx_begin_2 - idx_end_1 - 1;

	double *values1 = new double[n_values_r1];
	double *values2 = new double[n_values_r2];
	double *valuesBewtween = new double[n_values_between];
	int *colInd1 = new int[n_values_r1];
	int *colInd2 = new int[n_values_r2];
	int *colIndBewtween = new int[n_values_between];

	memcpy(values1, &(matrix.values[idx_begin_1]), n_values_r1 * sizeof(double));
	memcpy(values2, &(matrix.values[idx_begin_2]), n_values_r2 * sizeof(double));
	memcpy(valuesBewtween, &(matrix.values[idx_end_1 + 1]), n_values_between * sizeof(double));

	memcpy(colInd1, &(matrix.col_ind[idx_begin_1]), n_values_r1 * sizeof(int));
	memcpy(colInd2, &(matrix.col_ind[idx_begin_2]), n_values_r2 * sizeof(int));
	memcpy(colIndBewtween, &(matrix.col_ind[idx_end_1 + 1]), n_values_between * sizeof(int));

	memcpy(&(matrix.values[idx_begin_1]), values2, n_values_r2 * sizeof(double));
	memcpy(&(matrix.values[idx_begin_1 + n_values_r2]), valuesBewtween, n_values_between * sizeof(double));
	memcpy(&(matrix.values[idx_begin_1 + n_values_r2 + n_values_between]), values1, n_values_r1 * sizeof(double));

	memcpy(&(matrix.col_ind[idx_begin_1]), colInd2, n_values_r2 * sizeof(int));
	memcpy(&(matrix.col_ind[idx_begin_1 + n_values_r2]), colIndBewtween, n_values_between * sizeof(int));
	memcpy(&(matrix.col_ind[idx_begin_1 + n_values_r2 + n_values_between]), colInd1, n_values_r1 * sizeof(int));

	matrix.row_ptr_end[r1] += n_values_r2 - n_values_r1;
	for (int ri = r1 + 1; ri < r2; ri++)
	{
		matrix.row_ptr_begin[ri] += n_values_r2 - n_values_r1;
		matrix.row_ptr_end[ri] += n_values_r2 - n_values_r1;
	}
	matrix.row_ptr_begin[r2] += n_values_r2 - n_values_r1;
}

void changeRows(double vector[], int r1, int r2)
{
	int temp = vector[r1];
	vector[r1] = vector[r2];
	vector[r2] = temp;
}

void multiplyMatrixByVector(const Matrix matrix, const double vector[], double result[])
{
	for (int i = 0; i < matrix.n_rows; i++)
	{
		result[i] = 0;

		for (int j = 0; j < matrix.n_cols; j++)
		{
			result[i] += vector[j] * get_value(matrix, i, j);
		}
	}
}

void printMatrix(Matrix matrix)
{
	printf("Matrix (%d, %d) with %d Non-zero values:\n\n", matrix.n_rows, matrix.n_cols, matrix.non_zeros);
	for (int i = 0; i < matrix.n_rows; i++)
	{
		for (int j = 0; j < matrix.n_cols; j++)
		{
			printf("%f, ", get_value(matrix, i, j));
		}
		printf("\n");
	}
	printf("\n");
}

void printVector(double vector[], int size)
{
	printf("Vector of length: %d\n", size);

	for (int i = 0; i < size; i++)
	{
		printf("[ %f ]\n", vector[i]);
	}
	printf("\n");
}