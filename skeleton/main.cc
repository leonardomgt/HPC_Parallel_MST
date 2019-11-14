#include "matrix.h"
#include "sparse.h"
#include <mpi.h>
#include <vector>
#include <limits>
#include <algorithm>

using namespace std;

/* Global variables holding the matrix data. To complete this assignment
 * you are requested to only use arrays and access these arrays with
 * subscripts. Do not use pointers.
 */

const int max_n_elements = 131072;
const int max_n_rows = 16384;

vector<vector<int>> initialSetupMST(int n_nodes);
void printGroupSets(vector<vector<int>> group_sets);

int main(int argc, char **argv)
{

	// * Initalize MPI Library

	MPI_Init(&argc, &argv);

	int n_procs;
	MPI_Comm_size(MPI_COMM_WORLD, &n_procs);

	int id;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);

	printf("\n\nN_PROCESSOS: %d\n", n_procs);
	printf("ID: %d\n\n", id);

	// * Check arguments.
	if (argc != 2)
	{
		fprintf(stderr, "usage: %s <filename>\n", argv[0]);
		return -1;
	}

	struct Matrix matrix = {
		.values = new double[max_n_elements],
		.col_ind = new int[max_n_elements],
		.row_ptr_begin = new int[max_n_rows],
		.row_ptr_end = new int[max_n_rows],
		.n_rows = 0,
		.n_cols = 0,
		.non_zeros = 0};

	bool ok = load_matrix(argv[1], max_n_elements, max_n_rows, matrix);
	if (!ok)
	{
		fprintf(stderr, "failed to load matrix.\n");
		return -1;
	}

	// For debugging, can be removed when implementation is finished.
	dump_nonzeros(matrix.n_rows, matrix.values, matrix.col_ind, matrix.row_ptr_begin, matrix.row_ptr_end);
	printMatrix(matrix);

	vector<vector<int>> group_sets = initialSetupMST(matrix.n_rows);

	printGroupSets(group_sets);

	for (size_t i = 0; i < group_sets.size(); i++)
	{
		vector<int> group_set = group_sets.at(i);

		double minWeight = numeric_limits<double>::max();
		int minSrc = -1, minDest = -1;

		for (size_t j = 0; j < group_set.size(); j++)
		{
			int row = group_set.at(j);

			for (int idx = matrix.row_ptr_begin[row]; idx <= matrix.row_ptr_end[row]; ++idx)
			{
				if (matrix.values[idx] < minWeight)
				{
					minWeight = matrix.values[idx];
					minSrc = row;
					minDest = matrix.col_ind[idx];
				}

				printf("Value (%d, %d): %f\n", row, matrix.col_ind[idx], matrix.values[idx]);
			}
		}
		printf("MinWeight (%d, %d): %f\n\n", minSrc, minDest, minWeight);
	}

	vector<int> a = {1, 3, 5, 6, 9};
	vector<int> b = {2, 3, 6, 7, 9};
	vector<int> r(a.size() + b.size());

	vector<int>::iterator it = set_union(a.begin(), a.end(), b.begin(), b.end(), r.begin());
	r.resize(it - r.begin());

	printf("Size: %d\n", r.size());
	for (size_t i = 0; i < r.size(); i++)
	{
		printf("%d: %d\n", i, r.at(i));
	}

	MPI_Finalize();
	return 0;
}

vector<vector<int>> initialSetupMST(int n_nodes)
{
	vector<vector<int>> result;

	for (int i = 0; i < n_nodes; i += 2)
	{
		vector<int> node;
		node.push_back(i);
		node.push_back(i + 1);
		result.push_back(node);
	}

	return result;
}

void printGroupSets(vector<vector<int>> group_sets)
{
	for (size_t i = 0; i < group_sets.size(); i++)
	{
		printf("Set %lu: {", i);
		for (size_t j = 0; j < group_sets.at(i).size(); j++)
		{
			printf("%d, ", group_sets.at(i).at(j));
		}
		printf("}\n");
	}
	printf("\n");
}
/*
int main(int argc, char *argv[])
{
	int n, myid, numprocs, i;
	double PI25DT = 3.141592653589793238462643;
	double mypi, pi, h, sum, x;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	for (size_t j = 1; j < 1000000; j *= 2)
	{
		// if (myid == 0)
		// {
		// 	printf("Enter number of intervals: (0 quits)");
		// 	scanf("%d", &n);
		// }
		n = j;
		MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if (n == 0)
			break;
		h = 1.0 / (double)n;
		sum = 0.0;
		for (i = myid + 1; i <= n; i += numprocs)
		{
			x = h * ((double)i - 0.5);
			sum += 4.0 / (1.0 + x * x);
		}
		mypi = h * sum;
		MPI_Reduce(&mypi, &pi, 1, MPI_DOUBLE, MPI_SUM, 0,
				   MPI_COMM_WORLD);
		if (myid == 0)
			printf("pi = approximately %.16f, Error is %.16f\n",
				   pi, fabs(pi - PI25DT));
	}
	MPI_Finalize();
	return 0;
}
*/