#ifndef __GRAPH_H__
#define __GRAPH_H__

#include <vector>
#include "matrix.h"

typedef struct Element Edge;

struct Graph
{
	int n_nodes;
	int n_edges;

	Edge *edges;

	Graph() {}
	Graph(int n_nodes, int max_n_edges) : n_nodes(n_nodes), n_edges(0), edges(new Edge[max_n_edges]) {}
	Graph(int n_nodes, int n_edges, Edge *edges) : n_nodes(n_nodes), n_edges(n_edges), edges(edges) {}
};

bool load_graph(const char *filename, Graph &graph);

void print_graph(const Graph &graph);

double total_mst_weight(const Graph &graph);

void print_graph_info(const Graph &graph);

#endif /* __GRAPH_H__ */
