/*
 * @file headers.h
 *      
 * @author: Fabio Cruz B. de Albuquerque
 * 
 * @date: 09/29/2014
 */

#ifndef HEADERS_H
#define HEADERS_H

#include <iostream>
#include <string>
#include <unistd.h>
#include <utility>
#include <vector>

#include "../include/hungarian.h"

#include "Util.h"

struct Node {
    unsigned cost;
    int n;
    std::vector<int> route;
    std::vector<std::pair<int,int> > arrows;
    std::vector<std::pair<int,int> > prohibited;
};

const int inf = 123456;

Node dfs(std::vector<Node>& nodes);
Node bfs(std::vector<Node>& nodes);
Node bestb(std::vector<Node>& nodes);

double ** copyMatrix2Double(const int ** matrix, const unsigned dim);
inline bool isNewUB(Node& node, const unsigned dim, const unsigned ub) {
    return node.route.size() == dim + 1 && node.cost < ub;
}

inline bool isValidCH(std::vector<int>& cycle, const unsigned dim) {
    return cycle.size() == dim + 1;
}

inline bool isRootOptimal(Node& root, const unsigned dim, unsigned& lb, unsigned& ub) {
    // define lower bound
    lb = root.cost;
    // assignment equivalente a um tour completo = solucao viavel ao TSP
    if (isValidCH(root.route, dim)) {
        ub = root.cost;
        return true;
    }
    return false;
}

std::vector<int> getVectorSolution(int ** matrix, int dim);
void hungarian(Node& nodeCurr, double ** matrix, const int dim);
void bnb(std::vector<Node>& nodes, const int ** matrix, unsigned dim, unsigned& lb, unsigned& ub);
void bnb(std::vector<int>& bestRoute, std::vector<Node>& nodes, const int ** matrix, unsigned dim, unsigned& lb, unsigned& ub);
void initBranchAndBound(const int ** matrix, const unsigned dim);
void printNode(const Node& node);
void verifyCycle(std::vector<int> &sol, std::vector< std::pair<int,int> > &cycleArrows, unsigned dim);
template<typename type> void free(type ** matrix, const unsigned dim);

#endif 
