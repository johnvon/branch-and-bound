/*
 * @file prim.h
 *      
 * @author: Fabio Cruz B. de Albuquerque
 * 
 * Algoritmo de prim adaptado de 
 * http://geeksforgeeks.org/greedy-algorithms-set-5-prims-minimum-spanning-tree-mst-2
 * 
 * @date: 04/10/14
 */

#include "../include/prim.h"

// A utility function to find the vertex with minimum key value, from
// the set of vertices not yet included in MST
template<typename type>
int minKey(type key[], bool mstSet[], const unsigned dim) {
    // Initialize min value
    int min_index = -1;
    type min = inf+1;

    for (unsigned v = 1; v < dim; v++)
        if (mstSet[v] == false && min > key[v])
            min = key[v], min_index = v;

    return min_index;
}

// Function to construct and print MST for a graph represented using adjacency
// matrix representation
template<typename type>
void prim1Tree(const unsigned dim, type ** graph, unsigned * degree, bool ** sol1Tree, type &cost) {
    int parent[dim]; // Array to store constructed MST
    type key[dim];   // Key values used to pick minimum weight edge in cut
    bool mstSet[dim];  // To represent set of vertices not yet included in MST

    unsigned count, i, v;
    int u;

    // Initialize all keys as INFINITE
    for (i = 0; i < dim; i++)
        key[i] = inf, mstSet[i] = false, degree[i] = 0, parent[i] = 0;

    // Always include first 1st vertex in MST.
    key[1] = 0;     // Make key 0 so that this vertex is picked as first vertex
    parent[1] = -1; // First node is always root of MST

    // The MST will have dim vertices
    for (count = 0; count < dim-2; count++) {
        // Pick thd minimum key vertex from the set of vertices
        // not yet included in MST

        u = minKey<type>(key, mstSet, dim);
        // solucao inviavel, descarta construcao de MST
        if (u == -1)
            return;

        // Add the picked vertex to the MST Set
        mstSet[u] = true;

        // Update key value and parent index of the adjacent vertices of
        // the picked vertex. Consider only those vertices which are not yet
        // included in MST
        for (v = 1; v < dim; v++)

            // graph[u][v] is non zero only for adjacent vertices of m
            // mstSet[v] is false for vertices not yet included in MST
            // Update the key only if graph[u][v] is smaller than key[v]
            if (graph[u][v] != inf && mstSet[v] == false && graph[u][v] <  key[v])
                parent[v]  = u, key[v] = graph[u][v];
    }

    for (i = 2; i < dim; i++) {
        degree[i]++, degree[parent[i]]++;
        sol1Tree[i][parent[i]] = sol1Tree[parent[i]][i] = true;
        cost += graph[i][parent[i]];
    }

    //    delete[] parent;
    //    delete[] key;
    //    delete[] mstSet;
}

template int minKey<int>(int key[], bool mstSet[], const unsigned dim);
template int minKey<double>(double key[], bool mstSet[], const unsigned dim);

template void prim1Tree<int>(const unsigned dim, int ** graph, unsigned * degree, bool ** sol1Tree, int &cost);
template void prim1Tree<double>(const unsigned dim, double ** graph, unsigned * degree, bool ** sol1Tree, double &cost);

