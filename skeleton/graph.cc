#include "graph.h"
#include <stdio.h>

bool load_graph(const char *filename, Graph &graph)
{
	int n_rows, n_cols;

	if (!read_matrix_market(filename, graph.edges, n_rows, n_cols))
		return false;

	if (n_rows != n_cols)
		return false;

	graph.n_nodes = n_rows;

	return true;
}

void print_graph(const Graph &graph)
{

	printf("\n\nNum nodes: %d\n", graph.n_nodes);
	printf("Num edges: %lu\n", graph.edges.size());

	if (graph.edges.size() > 100)
		return;

	bool big = graph.n_nodes > 26 ? true : false;

	for (size_t i = 0; i < graph.edges.size(); i++)
	{
		Edge edge = graph.edges.at(i);
		if (!big)
			printf("Edge: %c --> %c, w = %f\n", 65 + edge.row, 65 + edge.col, edge.val);
		else
			printf("Edge: %d --> %d, w = %f\n", edge.row, edge.col, edge.val);
	}
}

double total_mst_weight(const Graph &graph)
{
	double sum = 0;

	for (size_t i = 0; i < graph.edges.size(); i++)
	{
		sum += graph.edges.at(i).val;
	}

	return sum;
}