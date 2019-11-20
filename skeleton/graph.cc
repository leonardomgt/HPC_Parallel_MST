#include "graph.h"
#include <stdio.h>
#include <cstring>

bool load_graph(const char *filename, Graph &graph)
{
	int n_srcs, n_dests;

	std::vector<Edge> vectorElems;

	if (!read_matrix_market(filename, vectorElems, n_srcs, n_dests))
		return false;

	if (n_srcs != n_dests)
		return false;

	graph.n_edges = vectorElems.size();
	graph.n_nodes = n_srcs;

	graph.edges = new Edge[graph.n_edges];

	memcpy(graph.edges, vectorElems.data(), graph.n_edges * sizeof(Edge));

	return true;
}

void print_graph(const Graph &graph)
{

	printf("\nNum nodes: %d\n", graph.n_nodes);
	printf("Num edges: %d\n", graph.n_edges);

	if (graph.n_edges > 100)
		return;

	bool big = graph.n_nodes > 26 ? true : false;

	for (int i = 0; i < graph.n_edges; i++)
	{
		Edge edge = graph.edges[i];
		if (!big)
			printf("Edge: %c --> %c, w = %f\n", 65 + edge.row, 65 + edge.col, edge.val);
		else
			printf("Edge: %d --> %d, w = %f\n", edge.row, edge.col, edge.val);
	}
}

double total_mst_weight(const Graph &graph)
{
	double sum = 0;

	for (int i = 0; i < graph.n_edges; i++)
	{
		sum += graph.edges[i].val;
	}

	return sum;
}

void print_graph_info(const Graph &graph)
{
	printf("\nNum nodes: %d\n", graph.n_nodes);
	printf("Num edges: %d\n", graph.n_edges);
	printf("Total MST weight: %f\n\n", total_mst_weight(graph));
}