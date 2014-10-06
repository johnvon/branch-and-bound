/*
 * @file headers.h
 *      
 * @author: Fabio Cruz B. de Albuquerque
 * 
 * @date: 09/29/2014
 */

#ifndef HEADERS_H
#define HEADERS_H

#include <iomanip>
#include <iostream>
#include <list>
#include <string>
#include <unistd.h>
#include <utility>
#include <vector>

#include "../include/hungarian.h"
#include "../include/prim.h"

#include "Util.h"

struct Node {
    unsigned cost;
    int n;
    std::vector<int> route;
    std::vector<std::pair<int,int> > arrows;
    std::vector<std::pair<int,int> > prohibited;
};


Node dfs(std::vector<Node>& nodes);
Node bfs(std::vector<Node>& nodes);
Node bestb(std::vector<Node>& nodes);

bool ** newBoolMatrix(const unsigned dim);

template<typename from, typename to> to ** copyMatrixFromTo(const from ** matrix, const unsigned dim, const unsigned s = 0);

template<typename type> void free(type ** matrix, const unsigned dim);

inline double gap(const unsigned lb, const unsigned ub) {
    return (1.0 - lb*1.0/ub)*100;
}

inline void doLog(unsigned ub, unsigned lb, unsigned size, unsigned long count, std::string strat) {
    std::cout << "UB: " << ub << std::setw(4) << " LB: " << lb << std::setw(4) 
        << " GAP: " <<  gap(lb, ub) << "%" << std::setw(4) << " Numero de nos abertos: " << size 
        << std::setw(4) << " Iteracao: " << count << " " << strat << std::endl;

}

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

inline void updateLB(std::vector<Node>& nodes, unsigned& lb) {
    unsigned min = inf;
    for (unsigned j = 0; j < nodes.size(); j++) {
        if (nodes[j].cost < min) {
            min = nodes[j].cost;
        }
    }
    if (min != (unsigned) inf && min > lb) {
        lb = min;
    }
}

bool isFeasible(const unsigned dim, unsigned * degree, unsigned& k, unsigned& z);

template<typename type> std::vector<int> getVectorSolution(type ** matrix, int dim);

std::vector<int> get1TreeVectorSolution(bool ** matrix, unsigned dim);

void hungarian(Node& nodeCurr, double ** matrix, const unsigned dim);

void prim(const unsigned dim, unsigned& cost, int ** graph, bool ** sol1Tree, unsigned * degree);

void oneTree(Node& node, int ** matrix, const unsigned dim);

Node rootBBHung(const int ** matrix, const unsigned dim);

Node rootBB1Tree(const int ** matrix, const unsigned dim);

void initBranchAndBound(const int ** matrix, const unsigned dim, unsigned b = 0);

void bnb(std::vector<int>& bestRoute, std::vector<Node>& nodes, const int ** matrix, const unsigned dim, unsigned& lb, unsigned& ub, unsigned b, unsigned x);

void printNode(const Node& node);

void printDegrees(unsigned * degree, const unsigned dim);

void verifyCycle(std::vector<int> &sol, std::vector< std::pair<int,int> > &cycleArrows, unsigned dim);
#endif 
