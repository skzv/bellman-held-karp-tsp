#include <iostream>
#include <unordered_set>
#include <map>
#include <algorithm>
#include <fstream>

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

/** Binary representation of a subset S of vector space V = {1, 2, ..., n} of size n where inclusion of element i is
 * represented by the bit at the i'th position. */
bool isElementInSet(int element, int set) {
    return ((set >> element) & 1) == 1;
}

int removeElementFromSet(int element, int set) {
    if (isElementInSet(element, set)) {
        return set ^ (1 << element);
    }
    return set;
}

/** Binary representation of a subset S of vector space V = {1, 2, ..., n} of size n where inclusion of element i is
 * represented by the bit at the i'th position. */
std::unordered_set<int> generateAllSubsetsOfSize(int vectorSpaceSize, int subsetSize) {
    // TODO: memoize this
    std::unordered_set<int> sets;

    if (subsetSize == 0) {
        sets.insert(0);
        return sets;
    }

    if (subsetSize == 1) {
        for (int i = 0; i < vectorSpaceSize; i++) {
            sets.insert(1 << i);
        }
        return sets;
    }

    std::unordered_set<int> subsets = generateAllSubsetsOfSize(vectorSpaceSize, subsetSize - 1);

    for (int subset: subsets) {
        for (int i = 0; i < vectorSpaceSize; i++) {
            if (!isElementInSet(i, subset)) {
                sets.insert((1 << i) | subset);
            }
        }
    }

    return sets;
}

/** Binary representation of a subset S of vector space V = {1, 2, ..., n} of size n where inclusion of element i is
 * represented by the bit at the i'th position. */
std::unordered_set<int> generateAllSubsetsIncludingFirstElementOfSize(int vectorspaceSize, int subsetSize) {
    std::unordered_set<int> bitmasks = generateAllSubsetsOfSize(vectorspaceSize - 1,
                                                                subsetSize - 1);
    std::unordered_set<int> sets;
    for (int bitmask: bitmasks) {
        sets.insert(bitmask << 1 | 1);
    }

    return sets;
}

/** Returns the shortest path distance that visits every node once. */
float BellmanHeldKarpTsp(UndirectedGraph g) {
    // note: index vertices from 0 to (n - 1)

    int n = g.getNumVertices();
    int V = (int) pow(2, n) - 1; // vector space (set including every element)

    // sub-problems (1 e S, |S| >= 2, j e V - {0})
    // (only sub-problems with j e S are ever used)
    float A[V + 1][n - 1];

    // base cases (|S| = 2)
    for (int j = 1; j < n; j++) {
        int S = 1 | (1 << j);
        A[S][j] = g.getEdgeCost(0, j);
    }

    // systematically solve all sub-problems
    for (int s = 3; s <= n; s++) {
        // for S with |S| = s and 1 e S
        for (int S: generateAllSubsetsIncludingFirstElementOfSize(n, s)) {
            // for j e S - {1}
            for (int j = 1; j < n; j++) {
                // if j is in set S
                if (isElementInSet(j, S)) {
                    // search for minPath from k (where k is in S - {j}, k != 1, j) to j
                    float minPath = std::numeric_limits<float>::max();
                    for (int k = 1; k < n; k++) {
                        if (k == j) {
                            continue;
                        }
                        if (isElementInSet(k, S)) {
                            int setWithoutJ = removeElementFromSet(j, S);
                            minPath = std::min(minPath, A[setWithoutJ][k] + g.getEdgeCost(j, k));
                        }
                    }
                    A[S][j] = minPath;
                }
            }
        }
    }

    // compute optimal tour cost
    float minCost = std::numeric_limits<float>::max();
    for (int j = 1; j < n; j++) {
        minCost = std::min(minCost, A[V][j] + g.getEdgeCost(j, 0));
    }
    return minCost;
}

float calculateEuclideanDistance(float x0, float y0, float x1, float y1) {
    return std::sqrt(std::pow(x1 - x0, 2) + std::pow(y1 - y0, 2));
}

UndirectedGraph readGraph() {
    // TODO: read graph from input file
    UndirectedGraph g;

    // hardcoded values for testing
//    g.addEdge(0, 1, 1);
//    g.addEdge(1, 2, 5);
//    g.addEdge(2, 3, 2);
//    g.addEdge(3, 0, 7);
//    g.addEdge(1, 3, 3);
//    g.addEdge(0, 2, 4);

    std::ifstream input("data.txt");

    int numVertices;
    input >> numVertices;

    // positions of vertices
    std::pair<float, float> p[numVertices];
    float x, y;

    for (int i = 0; i < numVertices; i++) {
        input >> x >> y;
        printf("%f, %f \n", x, y);
        p[i] = {x, y};
    }

    for (int i = 0; i < numVertices; i++) {
        for (int j = 0; j < numVertices; j++) {
            g.addEdge(i, j, calculateEuclideanDistance(p[i].first, p[j].first, p[i].second, p[j].second));
        }
    }

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