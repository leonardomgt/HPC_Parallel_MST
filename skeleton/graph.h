#ifndef __GRAPH_H__
#define __GRAPH_H__

#include <vector>
#include "matrix.h"

typedef Element Edge;

struct Graph
{
	int n_nodes;
	std::vector<Edge> edges;

	Graph() : n_nodes(0), edges(std::vector<Edge>()) {}
	Graph(int n_nodes, std::vector<Edge> edges) : n_nodes(n_nodes), edges(edges) {}
	Graph(int n_nodes) : n_nodes(n_nodes), edges(std::vector<Edge>()) {}
};

// struct Vertex
// {
// 	std::vector<Edge> edges;
// };

// struct GraphAlt
// {
// 	std::vector<Vertex> vertices;
// };

bool load_graph(const char *filename, Graph &graph);

void print_graph(const Graph &graph);

double total_mst_weight(const Graph &graph);

#endif /* __GRAPH_H__ */
