#include <iostream>
#include <unordered_set>
#include<map>

class UndirectedGraph {
public:
    void addEdge(int u, int v, float c) {
        vertices.insert(u);
        vertices.insert(v);
        edges[{u, v}] = c;
    }

    int getNumVertices() {
        return vertices.size();
    }

    int getNumEdges() {
        return edges.size();
    };

    float getEdgeCost(int u, int v) {
        // since this is an undirected graph we check for an edge in both {u, v} or {v, u} directions
        std::map<std::pair<int, int>, float>::iterator edge;

        // check for u,v
        edge = edges.find({u, v});
        if (edge != edges.end()) {
            return edge->second;
        }

        // else check for v,u
        edge = edges.find({v, u});
        if (edge != edges.end()) {
            return edge->second;
        }

        // else return +inf
        return std::numeric_limits<float>::max();;
    }

private:
    std::unordered_set<int> vertices;
    std::map<std::pair<int, int>, float> edges;
};

/** Returns the shortest path distance that visits every node once. */
float BellmanHeldKarpTsp(UndirectedGraph g) {
    return 123.45;
}

UndirectedGraph readGraph() {
    // TODO: read graph from input file
    UndirectedGraph g;

    g.addEdge(0, 1, 2.1);
    g.addEdge(1, 2, 3.4);

    return g;
}

int main() {
    UndirectedGraph g = readGraph();
    std::cout << "Num vertices: " << g.getNumVertices() << std::endl;
    std::cout << "Num edges: " << g.getNumEdges() << std::endl;

    float shortestPathLength = BellmanHeldKarpTsp(g);
    std::cout << "Shortest travelling salesman path length: " << shortestPathLength << " (" << (int) shortestPathLength
              << ")" << std::endl;

    return 0;
}