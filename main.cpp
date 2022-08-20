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

    float getEdge(int u, int v) {
        // TODO: implement
        // check for u,v
        //   return edges[{u,v}]
        // else check for v,u
        //   return edges[{v,u}]
        // else return +inf
        return 0;
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