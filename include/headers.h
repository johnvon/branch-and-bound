/*
 * @file headers.h
 *      
 * @author: Fabio Cruz B. de Albuquerque
 * 
 * @date: 09/29/2014
 */

#include <exception>
#include <iostream>
#include <list>
#include <sstream>
#include <string>
#include <unistd.h>
#include <utility>
#include <vector>

#include "../include/hungarian.h"

#include "Instance.h"
#include "Util.h"

struct Node {
    unsigned solution;
    int n;
    std::vector<int> route;
    std::vector<std::pair<int,int>> arrows;
    std::vector<std::pair<int,int>> prohibited;

    Node& operator=(Node l) // copy assignment
    {
        n = l.n;
        solution = l.solution;

        route.swap(l.route);
        arrows.swap(l.arrows);
        prohibited.swap(l.prohibited);

        return *this;
    }
};

double ** copyMatrix2Double(const int ** matrix, const unsigned dim);
inline bool isValidCH(std::vector<int>& vCycle, const unsigned dim) {
    return vCycle.size() == dim + 1;
}
std::vector<int> getVectorSolution(int ** matrix, int dim);
std::vector<int> hungarian(Node& nodeCurr, const double ** matrix, const int dim);
Node bnb(std::vector<Node>& nodes, const int ** matrix, unsigned dim, unsigned& lb, unsigned& ub);

template<typename type> 
void free(type ** matrix, const unsigned dim);

void initBranchAndBound(const int ** matrix, const unsigned dim);

bool isRootOptimal(Node& root, const unsigned dim, std::vector<int>& vCycle, unsigned& lb, unsigned& ub);
void printNode(const Node& node);
void verifyCycle(std::vector<int> &sol, std::vector<int> &cycle, std::vector< std::pair<int,int> > &cycleArrows, unsigned dim);
