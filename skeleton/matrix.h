#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <vector>

const int INVALID_VALUE = -1;

struct Element
{
	int row, col;
	double val;

	Element() : row(INVALID_VALUE), col(INVALID_VALUE) {}
	Element(int row, int col, double val) : row(row), col(col), val(val) {}
	bool isInvalid() { return row == INVALID_VALUE || col == INVALID_VALUE; }
};

bool read_matrix_market(const char *filename,
						std::vector<Element> &elements,
						int &n_rows, int &n_cols);

void dump_nonzeros(const int n_rows,
				   const double values[],
				   const int col_ind[],
				   const int row_ptr_begin[],
				   const int row_ptr_end[]);

bool load_matrix_market(const char *filename,
						const int max_n_elements,
						const int max_n_rows,
						int &nnz,
						int &n_rows,
						int &n_cols,
						double values[],
						int col_ind[],
						int row_ptr_begin[],
						int row_ptr_end[]);

#endif /* __MATRIX_H__ */
