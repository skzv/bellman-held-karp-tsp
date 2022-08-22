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

/**
 * Binary representation of a subset S of vector space V = {1, 2, ..., n} of size n where inclusion of element i is
 * represented by the bit at the i'th position.
 *
 * @param element
 * @param set
 * @return true if element is in set, false otherwise
 */
bool isElementInSet(int element, int set) {
    return ((set >> element) & 1) == 1;
}

int removeElementFromSet(int element, int set) {
    if (isElementInSet(element, set)) {
        return set ^ (1 << element);
    }
    return set;
}

/**
 * Binary representation of a subset S of vector space V = {1, 2, ..., n} of size n where inclusion of element i is
 * represented by the bit at the i'th position.
 *
 * @param vectorSpaceSize
 * @param subsetSize
 * @return binary represenations of all subsets of V of size |S|
 */
std::unordered_set<int> generateAllSubsetsOfSize(int vectorSpaceSize, int subsetSize) {
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

/**
 * Binary representation of a subset S of vector space V = {1, 2, ..., n} of size n where inclusion of element i is
 * represented by the bit at the i'th position.
 *
 * @param vectorSpaceSize
 * @param subsetSize
 * @return all subsets of vector space V of size |S| containing element {1}
 */
std::unordered_set<int> generateAllSubsetsIncludingFirstElementOfSize(int vectorSpaceSize, int subsetSize) {
    std::unordered_set<int> bitmasks = generateAllSubsetsOfSize(vectorSpaceSize - 1,
                                                                subsetSize - 1);
    std::unordered_set<int> sets;
    for (int bitmask: bitmasks) {
        sets.insert(bitmask << 1 | 1);
    }

    return sets;
}

/**
 * Returns the shortest path distance that visits every node once.
 *
 * @param g
 * @return shortest path distance visiting every node of graph g
 */
float BellmanHeldKarpTsp(UndirectedGraph g) {
    // note: index vertices from 0 to (n - 1)

    int n = g.getNumVertices();
    int V = (int) pow(2, n) - 1; // vector space (set including every element)

    // sub-problems (1 e S, |S| >= 2, j e V - {0})
    // (only sub-problems with j e S are ever used)
    std::map<std::pair<int, int>, float> A;

    // base cases (|S| = 2)
    for (int j = 1; j < n; j++) {
        int S = 1 | (1 << j);
        A[{S, j}] = g.getEdgeCost(0, j);
    }

    // systematically solve all sub-problems
    for (int s = 3; s <= n; s++) {
        std::cout << "(" << s << "/" << n << ")" << std::endl;
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
                            minPath = std::min(minPath, A[{setWithoutJ, k}] + g.getEdgeCost(j, k));
                        }
                    }
                    A[{S, j}] = minPath;
                }
            }
        }
    }

    // compute optimal tour cost
    float minCost = std::numeric_limits<float>::max();
    for (int j = 1; j < n; j++) {
        minCost = std::min(minCost, A[{V, j}] + g.getEdgeCost(j, 0));
    }
    return minCost;
}

/**
 * Calculates euclidean distance between two points.
 *
 * @param p0
 * @param p1
 * @return euclidean distance p0 and p1
 */
float calculateEuclideanDistance(std::pair<float, float> p0, std::pair<float, float> p1) {
    return (float) std::sqrt(std::pow(p0.first - p1.first, 2) + std::pow(p0.second - p1.second, 2));
}

UndirectedGraph readGraph() {
    UndirectedGraph g;

    const bool TEST = false;

    if (TEST) {
        // hardcoded values for testing
        g.addEdge(0, 1, 1);
        g.addEdge(1, 2, 5);
        g.addEdge(2, 3, 2);
        g.addEdge(3, 0, 7);
        g.addEdge(1, 3, 3);
        g.addEdge(0, 2, 4);
    } else {
        std::ifstream input("data.txt");

        int numVertices;
        input >> numVertices;

        // positions of vertices
        std::pair<float, float> p[numVertices];
        float x, y;

        for (int i = 0; i < numVertices; i++) {
            input >> x >> y;
            p[i] = {x, y};
        }

        for (int i = 0; i < numVertices; i++) {
            for (int j = 0; j < numVertices; j++) {
                g.addEdge(i, j, calculateEuclideanDistance(p[i], p[j]));
            }
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