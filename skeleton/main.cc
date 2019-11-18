#include "matrix.h"
#include "graph.h"
#include <mpi.h>
#include <vector>
#include <limits>
#include <algorithm>
#include <unordered_map>
#include <time.h>

using namespace std;

/* Global variables holding the matrix data. To complete this assignment
 * you are requested to only use arrays and access these arrays with
 * subscripts. Do not use pointers.
 */

const int max_n_elements = 131072;
const int max_n_rows = 16384;

// vector<vector<int>> initialSetupMST(int n_nodes);
// void printGroupSets(vector<vector<int>> group_sets);
Graph mst_sequential(const Graph &graph);
int find_root(unordered_map<int, int> &mst_i, int vertex);
bool contract_subtrees(unordered_map<int, int> &mst_i, int root1, int root2);

int main(int argc, char **argv)
{

	// // * Initalize MPI Library

	// MPI_Init(&argc, &argv);

	// int n_procs;
	// MPI_Comm_size(MPI_COMM_WORLD, &n_procs);

	// int id;
	// MPI_Comm_rank(MPI_COMM_WORLD, &id);

	// printf("\n\nN_PROCESSOS: %d\n", n_procs);
	// printf("ID: %d\n\n", id);

	// * Check arguments
	if (argc != 2)
	{
		fprintf(stderr, "usage: %s <filename>\n", argv[0]);
		return -1;
	}

	Graph graph;

	clock_t begin = clock();

	bool ok = load_graph(argv[1], graph);
	if (!ok)
	{
		fprintf(stderr, "failed to load matrix.\n");
		return -1;
	}

	clock_t end = clock();
	double time_read = (double)(end - begin) / CLOCKS_PER_SEC;

	printf("Reading time: %f sec\n", time_read);

	print_graph(graph);

	begin = clock();

	Graph mst = mst_sequential(graph);

	end = clock();
	double time_execute = (double)(end - begin) / CLOCKS_PER_SEC;

	printf("Execution time: %f sec\n", time_execute);

	print_graph(mst);
	printf("\nTotal MST weight: %f\n", total_mst_weight(mst));

	// vector<vector<int>> group_sets = initialSetupMST(matrix.n_rows);

	// printGroupSets(group_sets);

	// for (size_t i = 0; i < group_sets.size(); i++)
	// {
	// 	vector<int> group_set = group_sets.at(i);

	// 	double minWeight = numeric_limits<double>::max();
	// 	int minSrc = -1, minDest = -1;

	// 	for (size_t j = 0; j < group_set.size(); j++)
	// 	{
	// 		int row = group_set.at(j);

	// 		for (int idx = matrix.row_ptr_begin[row]; idx <= matrix.row_ptr_end[row]; ++idx)
	// 		{
	// 			if (matrix.values[idx] < minWeight)
	// 			{
	// 				minWeight = matrix.values[idx];
	// 				minSrc = row;
	// 				minDest = matrix.col_ind[idx];
	// 			}

	// 			printf("Value (%d, %d): %f\n", row, matrix.col_ind[idx], matrix.values[idx]);
	// 		}
	// 	}
	// 	printf("MinWeight (%d, %d): %f\n\n", minSrc, minDest, minWeight);
	// }

	// vector<int> a = {1, 3, 5, 6, 9};
	// vector<int> b = {2, 3, 6, 7, 9};
	// vector<int> r(a.size() + b.size());

	// vector<int>::iterator it = set_union(a.begin(), a.end(), b.begin(), b.end(), r.begin());
	// r.resize(it - r.begin());

	// printf("Size: %lu\n", r.size());
	// for (size_t i = 0; i < r.size(); i++)
	// {
	// 	printf("%lu: %d\n", i, r.at(i));
	// }

	// MPI_Finalize();
	return 0;
}

Graph mst_sequential(const Graph &graph)
{
	Graph mst_result(graph.n_nodes);

	// * Maps <original_node, parent> as a tree (mst_i stands for intermediate minimum spanning tree)
	unordered_map<int, int> mst_i;

	// * Initialize mst_i with N subtrees, one for each original node
	for (int i = 0; i < graph.n_nodes; i++)
	{
		mst_i.emplace(i, i);
	}

	int nSubtrees = graph.n_nodes;
	int last_nSubtrees = 0;

	while (nSubtrees > 1 && last_nSubtrees != nSubtrees)
	{
		last_nSubtrees = nSubtrees;

		// * Maps <super_node, minimum_edge>
		unordered_map<int, Edge> minimum_edge;

		// * Iterate all the edges
		for (size_t i = 0; i < graph.edges.size(); i++)
		{
			Edge edge = graph.edges.at(i);

			int rootSrc = find_root(mst_i, edge.row);
			int rootDest = find_root(mst_i, edge.col);

			// * If source and destination have the same root, they're already contracted. Continue.
			if (rootSrc == rootDest)
				continue;

			auto it = minimum_edge.find(rootSrc);

			if (it == minimum_edge.end())
			{
				minimum_edge.emplace(rootSrc, edge);
			}
			else
			{
				if (it->second.val > edge.val)
				{
					it->second = edge;
				}
			}
		}

		// for (auto it = mst_i.begin(); it != mst_i.end(); ++it)
		// 	printf("SN(%d) = %d\n", it->first, it->second);

		// for (auto it = minimum_edge.begin(); it != minimum_edge.end(); ++it)
		// 	printf("ME(%d) = Edge(%d, %d, %f)\n", it->first, it->second.row, it->second.col, it->second.val);

		for (auto it = minimum_edge.begin(); it != minimum_edge.end(); ++it)
		{
			int rootSrc = find_root(mst_i, it->second.row);
			int rootDest = find_root(mst_i, it->second.col);

			if (rootSrc != rootDest)
			{
				mst_result.edges.push_back(it->second);
				if (contract_subtrees(mst_i, rootSrc, rootDest))
					nSubtrees--;
			}
		}

		// for (auto it = mst_i.begin(); it != mst_i.end(); ++it)
		// 	printf("SN(%d) = %d\n", it->first, it->second);

		// printf("\nlast_nSubtrees: %d\n", last_nSubtrees);
		// printf("Subtrees: %d\n", nSubtrees);
	}

	return mst_result;
}

int find_root(unordered_map<int, int> &mst_i, int vertex)
{
	if (mst_i.at(vertex) == vertex)
		return vertex;
	else
		return find_root(mst_i, mst_i.at(vertex));
}

bool contract_subtrees(unordered_map<int, int> &mst_i, int root1, int root2)
{
	auto it1 = mst_i.find(root1);
	auto it2 = mst_i.find(root2);

	if (it1 != mst_i.end() && it2 != mst_i.end())
	{
		it1->second = root2;
		return true;
	}

	return false;
}

// vector<vector<int>> initialSetupMST(int n_nodes)
// {
// 	vector<vector<int>> result;

// 	for (int i = 0; i < n_nodes; i += 2)
// 	{
// 		vector<int> node;
// 		node.push_back(i);
// 		node.push_back(i + 1);
// 		result.push_back(node);
// 	}

// 	return result;
// }

// void printGroupSets(vector<vector<int>> group_sets)
// {
// 	for (size_t i = 0; i < group_sets.size(); i++)
// 	{
// 		printf("Set %lu: {", i);
// 		for (size_t j = 0; j < group_sets.at(i).size(); j++)
// 		{
// 			printf("%d, ", group_sets.at(i).at(j));
// 		}
// 		printf("}\n");
// 	}
// 	printf("\n");
// }
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