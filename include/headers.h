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
#include <cmath>
#include <string>
#include <unistd.h>
#include <utility>
#include <vector>

#include "../include/hungarian.h"
#include "../include/prim.h"

#include "Construct.h"
#include "LocalSearch.h"
#include "Util.h"

struct Node {
    double cost;
    int n;
    std::vector<int> route;
    std::vector<int> pi;
    std::vector<double> u;
    std::vector<std::pair<int,int> > arrows;
    std::vector<std::pair<int,int> > prohibited;
};

inline bool isNewUB(Node& node, const unsigned dim, const unsigned ub) {
    return node.route.size() == dim + 1 && node.cost <= ub;
}

inline bool isValidCH(std::vector<int>& cycle, const unsigned dim) {
    return cycle.size() == dim + 1;
}

inline double gap(const unsigned lb, const unsigned ub) {
    return (1.0 - lb*1.0/ub)*100;
}

inline void updateLB(std::list<Node>& nodes, unsigned& lb) {
    unsigned min = inf;
    for (std::list<Node>::iterator it = nodes.begin(); it != nodes.end(); ++it) {
        if (it->cost < min) {
            min = it->cost;
        }
    }
    if (min != (unsigned) inf && min > lb) {
        lb = min;
    }
}

inline unsigned sumsqr(double * subgradient, const unsigned dim) {
    unsigned i, t = 0;
    for (i = 0; i < dim; i++) {
        t += pow((2 - subgradient[i]), 2);
    }
    return t;
}

inline double sum(std::vector<double>& u) {
    double t = 0;
    for (std::vector<double>::iterator it = u.begin();
            it != u.end(); ++it) {
        t += *it;
    }
    return t;
}

/**
 * Best-first
 */
inline std::list<Node>::iterator bestb(std::list<Node>& nodes) {
    unsigned bb = 0, i, indexBB = 0;
    std::list<Node>::iterator it;
    for (it = nodes.begin(); it != nodes.end(); ++it) {
        if (it->cost > bb) {
            bb = it->cost;
            indexBB = i;
        }
    }
    it = nodes.begin();
    std::advance(it, indexBB);
    return it;
}

/**
 * Breadth-first
 */
inline std::list<Node>::iterator bfs(std::list<Node>& nodes) {
    return nodes.begin();
}

/**
 * Depth-first
 */
inline std::list<Node>::iterator dfs(std::list<Node>& nodes) {
    return --nodes.end();
}

/**
 * "Random"-first
 */
inline std::list<Node>::iterator randb(std::list<Node>& nodes) {
    std::list<Node>::iterator it = nodes.begin();
    std::advance(it, rand() % nodes.size());
    return it;
}



Node rootBB1Tree(const int ** matrix, const unsigned dim);
Node rootBBHung(const int ** matrix, const unsigned dim);
Node rootBBLR(const int ** matrix, const unsigned dim, unsigned ub);

void bnbHung(std::vector<int>& bestRoute, std::list<Node>& nodes, const int ** matrix,
        const unsigned dim, unsigned& lb, unsigned& ub, unsigned b);

void bnb1Tree(std::vector<int>& bestRoute, std::list<Node>& nodes, const int ** matrix,
        const unsigned dim, unsigned& lb, unsigned& ub, unsigned b);

void bnbLR(std::vector<int>& bestRoute, std::list<Node>& nodes, const int ** matrix,
        const unsigned dim, unsigned& lb, unsigned& ub, unsigned b);

bool ** newBoolMatrix(const unsigned dim);
bool isFeasible(const unsigned dim, unsigned * degree, unsigned& k);
bool isRootOptimal(Node& root, const unsigned dim, unsigned& lb, unsigned& ub);
std::vector<int> get1TreeVectorSolution(bool ** matrix, unsigned dim);
template<typename from, typename to> to ** copyMatrixFromTo(const from ** matrix, 
        const unsigned dim, const unsigned s = 0);
template<typename type> std::vector<int> getVectorSolution(type ** matrix, int dim);
template<typename type> void free(type ** matrix, const unsigned dim);
void bnb(std::vector<int>& bestRoute, std::list<Node>& nodes, const int ** matrix, 
        const unsigned dim, unsigned& lb, unsigned& ub, unsigned b, unsigned x);
void doLog(unsigned ub, unsigned lb, unsigned size, unsigned long count, 
        std::string strat);
void hungarian(Node& nodeCurr, double ** matrix, const unsigned dim);
void lagrangean(Node& nodeCurr, double ** cMatrix, bool ** sol1Tree, const unsigned dim, 
        unsigned ub);
void initBranchAndBound(const int ** matrix, const unsigned dim, unsigned b = 0, 
        unsigned x = 1);
template<typename type>
void oneTree(Node& node, type ** matrix, const unsigned dim, bool ** sol1Tree, 
        unsigned * degree);
void printDegrees(unsigned * degree, const unsigned dim);
void printNode(const Node& node);
void resetBoolMatrix(bool ** matrix, const unsigned dim);
void strategy(std::list<Node>::iterator (*curr)(std::list<Node>&), std::string& strat, 
        unsigned b);
void verifyCycle(std::vector<int> &sol, 
        std::vector< std::pair<int,int> > &cycleArrows, unsigned dim);

#endif 
