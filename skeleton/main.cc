#include "matrix.h"
#include "graph.h"
#include <mpi.h>
#include <vector>
#include <limits>
#include <algorithm>
#include <unordered_map>
#include <time.h>
#include <math.h>

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
Graph mst_parallel(const Graph &graph);
int find_root(unordered_map<int, int> &mst_i, int vertex);
int find_root(int *mst_forest, int vertex);
bool contract_subtrees(unordered_map<int, int> &mst_i, int root1, int root2);
void contract_subtrees(int *mst_forest, int root1, int root2);
MPI_Datatype create_edge_MPI_datatype();
void add_edge_to_graph(Graph &graph, Edge edge);
void merge_minimum_edges(Edge *edges, Edge *edges_rcv, int n_nodes);

int main(int argc, char **argv)
{
	// * Check arguments
	if (argc != 2)
	{
		fprintf(stderr, "usage: %s <filename>\n", argv[0]);
		exit(EXIT_FAILURE);
	}

	// * Initalize MPI Library
	MPI_Init(&argc, &argv);

	int n_procs;
	MPI_Comm_size(MPI_COMM_WORLD, &n_procs);

	int id;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);

	printf("\n\nN_PROCESSOS: %d\n", n_procs);
	printf("ID: %d\n\n", id);

	Graph graph;

	if (id == 0)
	{
		double start_reading = MPI_Wtime();

		bool ok = load_graph(argv[1], graph);

		double end_reading = MPI_Wtime();

		printf("Reading complete in %f s\n", end_reading - start_reading);

		if (!ok)
		{
			fprintf(stderr, "Failed to load graph.\n");
			return -1;
		}

		// print_graph(graph);
	}
	double start_exec = MPI_Wtime();

	// Graph mst = mst_sequential(graph);
	Graph mst = mst_parallel(graph);

	double end_exec = MPI_Wtime();

	if (id == 0)
	{
		printf("Execution complete in %f s\n", start_exec - end_exec);
		// print_graph(mst);
		print_graph_info(mst);
	}

	MPI_Finalize();
	return 0;
}

Graph mst_parallel(const Graph &graph)
{
	int id;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);

	int n_procs;
	MPI_Comm_size(MPI_COMM_WORLD, &n_procs);

	int n_nodes;
	int n_edges;

	if (id == 0)
	{
		n_nodes = graph.n_nodes;
		n_edges = graph.n_edges;
	}

	// * Broadcast the number of edges and nodes in the original graph.
	MPI_Bcast(&n_nodes, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&n_edges, 1, MPI_INT, 0, MPI_COMM_WORLD);

	Graph mst_result(n_nodes, n_nodes - 1);

	// * Maps <original_node, parent> as a tree (mst_i stands for intermediate minimum spanning tree)

	int *mst_forest = new int[n_nodes];

	if (id == 0)
	{
		for (int i = 0; i < n_nodes; i++)
		{
			mst_forest[i] = i;
		}
	}

	int edgesPerProcess = (n_edges + n_procs - 1) / n_procs;

	Edge *processPart = new Edge[edgesPerProcess];

	MPI_Datatype MPI_EDGE = create_edge_MPI_datatype();

	MPI_Scatter(
		graph.edges,
		edgesPerProcess, MPI_EDGE,
		processPart,
		edgesPerProcess, MPI_EDGE,
		0, MPI_COMM_WORLD);

	// * Adjust edgesPerProcess for the last process.
	if (id == n_procs - 1)
	{
		edgesPerProcess = n_edges % edgesPerProcess == 0 ? edgesPerProcess : n_edges % edgesPerProcess;
	}

	// * Process 0 initializes nSubtrees and broadcast to the others.
	int nSubtrees = graph.n_nodes;
	MPI_Bcast(&nSubtrees, 1, MPI_INT, 0, MPI_COMM_WORLD);

	int last_nSubtrees = 0;

	while (nSubtrees > 1 && last_nSubtrees != nSubtrees)
	{
		last_nSubtrees = nSubtrees;

		// * Broadcast the mst_forest.
		MPI_Bcast(mst_forest, n_nodes, MPI_INT, 0, MPI_COMM_WORLD);

		// * Maps <super_node, minimum_edge>
		Edge *minimum_edge = new Edge[n_nodes];

		// * Iterate through edges
		for (int i = 0; i < edgesPerProcess; i++)
		{
			Edge edge = processPart[i];

			int rootSrc = find_root(mst_forest, edge.row);
			int rootDest = find_root(mst_forest, edge.col);

			// * If source and destination have the same root, they're already contracted. Continue.
			if (rootSrc == rootDest)
				continue;

			if (minimum_edge[rootSrc].isInvalid() || edge.val < minimum_edge[rootSrc].val)
			{
				minimum_edge[rootSrc] = edge;
			}
		}

		// for (int i = 0; i < n_nodes; i++)
		// {
		// 	Edge edge = minimum_edge[i];
		// 	printf("min_edge i = %d: (%d, %d, %f)\n", i, edge.row, edge.col, edge.val);
		// }

		Edge *minimum_edge_rcv = new Edge[n_nodes];
		MPI_Status status;

		for (size_t k = 1; k <= ceil(log2(n_procs)); k++)
		{
			if (id % (int)pow(2, k) == pow(2, k - 1))
			{
				// printf("\nId %d sending to %f\n", id, id - pow(2, k - 1));
				MPI_Send(minimum_edge, n_nodes, MPI_EDGE, id - pow(2, k - 1), k, MPI_COMM_WORLD);
			}
			else if (id % (int)pow(2, k) == 0 && id + pow(2, k - 1) < n_procs)
			{
				// printf("\nId %d receiving from %f\n", id, id + pow(2, k - 1));
				MPI_Recv(minimum_edge_rcv, n_nodes, MPI_EDGE, id + pow(2, k - 1), k, MPI_COMM_WORLD, &status);

				// for (int i = 0; i < n_nodes; i++)
				// {
				// 	Edge edge = minimum_edge_rcv[i];
				// 	printf("min_edge_rcv i = %d: (%d, %d, %f)\n", i, edge.row, edge.col, edge.val);
				// }
				merge_minimum_edges(minimum_edge, minimum_edge_rcv, n_nodes);
			}
		}

		if (id == 0)
		{
			for (int i = 0; i < n_nodes; i++)
			{
				Edge edge = minimum_edge[i];
				// printf("MIN_EDGE FINAL i = %d: E(%d, %d, %f)\n", i, edge.row, edge.col, edge.val);

				int rootSrc = find_root(mst_forest, edge.row);
				int rootDest = find_root(mst_forest, edge.col);

				if (rootSrc != rootDest)
				{
					add_edge_to_graph(mst_result, edge);

					contract_subtrees(mst_forest, rootSrc, rootDest);
					nSubtrees--;
				}
			}
		}

		// * Broadcast nSubtrees
		MPI_Bcast(&nSubtrees, 1, MPI_INT, 0, MPI_COMM_WORLD);
	}

	return mst_result;
}

void merge_minimum_edges(Edge *edges, Edge *edges_rcv, int n_nodes)
{
	for (int i = 0; i < n_nodes; i++)
	{
		Edge edge = edges[i];
		Edge edge_rcv = edges_rcv[i];

		if (edge_rcv.isInvalid())
			continue;

		else if (edge.isInvalid() || edge_rcv.val < edge.val)
		{
			edges[i] = edge_rcv;
		}
	}
}

MPI_Datatype create_edge_MPI_datatype()
{
	int blockcount[3] = {1, 1, 1};
	MPI_Aint offsets[3] = {offsetof(Edge, row), offsetof(Edge, col), offsetof(Edge, val)};
	MPI_Datatype EdgeElementTypes[3] = {MPI_INT, MPI_INT, MPI_DOUBLE};
	MPI_Datatype MPI_EDGE;

	MPI_Type_create_struct(3, blockcount, offsets, EdgeElementTypes, &MPI_EDGE);
	MPI_Type_commit(&MPI_EDGE);

	return MPI_EDGE;
}

Graph mst_sequential(const Graph &graph)
{
	Graph mst_result(graph.n_nodes, graph.n_nodes - 1);

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
		for (int i = 0; i < graph.n_edges; i++)
		{
			Edge edge = graph.edges[i];

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
				add_edge_to_graph(mst_result, it->second);

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

void add_edge_to_graph(Graph &graph, Edge edge)
{
	graph.edges[graph.n_edges] = edge;
	graph.n_edges++;
}

int find_root(unordered_map<int, int> &mst_i, int vertex)
{
	if (mst_i.at(vertex) != vertex)
		mst_i.find(vertex)->second = find_root(mst_i, mst_i.at(vertex));

	return mst_i.at(vertex);
}

int find_root(int *mst_forest, int vertex)
{
	if (mst_forest[vertex] != vertex)
		mst_forest[vertex] = find_root(mst_forest, mst_forest[vertex]);

	return mst_forest[vertex];
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

void contract_subtrees(int *mst_forest, int root1, int root2)
{
	mst_forest[root1] = root2;
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