/*
 * @file main.cpp
 *      
 * @author: Fabio Cruz B. de Albuquerque
 * 
 * @date: 27/09/14
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
    /* 
     *     Node& operator=(Node l) // copy assignment
     *     {
     *         n = l.n;
     *         parent = l.parent;
     *         solution = l.solution;
     * 
     *         route.swap(l.route);
     *         arrows.swap(l.arrows);
     *         prohibited.swap(l.prohibited);
     * 
     *         return *this;
     *     }
     */

};

void printNode(const Node& node) {

    std::cout << std::endl << "Node " << &node << " " << node.n << "solution: " << &(node.solution) << " " << node.solution << std::endl;
    std::cout << "Route " << &node.route << std::endl;
    tsp::printRoute(node.route);
    std::cout << "Arrows " << &node.arrows << std::endl;
    for (unsigned i = 0; i < node.arrows.size(); i++) {
        std::cout << node.arrows[i].first << " " << node.arrows[i].second << ", ";
    }
    std::cout << std::endl << "Prohibited " << std::endl;
    for (unsigned i = 0; i < node.prohibited.size(); i++) {
        std::cout << node.prohibited[i].first << " " << node.prohibited[i].second << ", ";
    }
    std::cout << std::endl;
}

void verifyCycle(std::vector<int> &sol, std::vector<int> &cycle, std::vector< std::pair<int,int> > &cycleArrows, unsigned dim) {
    bool cDetected = false;
    int cSize = -1;
    unsigned cBegin = 0, cBeginTemp = 0, i;

    // busca menor ciclo na rota
    if (sol.size() > dim + 1) {
        for (i = 1; i < sol.size(); i++) {
            if (cDetected) {
                cBeginTemp = i;
                cDetected = false;
            } else {
                if (sol[i] == sol[cBeginTemp]) {
                    cDetected = true;
                    // se encontrou ciclo menor ou primeiro ciclo
                    if ((int) (i - cBeginTemp) < (int) cSize || (cSize == -1)) {
                        cBegin = cBeginTemp;
                        cSize = i - cBeginTemp;
                    }
                }
            }
        }
    }

    cycle.clear();
    cycleArrows.clear();

    // clientes/cidades pertencentes ao intervalo
    // e respectivos arcos
    if (sol.size() > dim + 1) {
        // intervalo da rota com ciclo
        for (i = cBegin; i <= cBegin + cSize; i++) {
            cycle.push_back(sol[i]);
        }

        for (i = 0; i < cycle.size() - 1; i++) {
            cycleArrows.push_back(std::make_pair(cycle[i], cycle[i+1]));
        }
    }
}

std::vector<int> getVectorSolution(int ** matrix, int dim) {
    int i, j;
    std::vector<bool> visited;
    std::vector<int> route;

    for (i = 0; i < dim; i++) {
        visited.push_back(false);
    }

    // partindo da primeira cidade
    route.push_back(0);
    for (i = 0; i < dim; i++) {
        if (visited[i])
            continue;

        if (route.back()!=i) {
            route.push_back(i);
        }
        visited[i] = true;

        for (j = 0; j < dim; j++) {
            if (matrix[i][j] == 1) {
                route.push_back(j);
                visited[true];
                i = j - 1;
                break;
            } 
        }
    }
    return route;
}

std::vector<int> best;
/**
 * Algoritmo Branch-and-bound
 * Usando busca em largura
 */
void bnb(std::vector<Node *>& nodes, const int ** matrix, unsigned dim, unsigned& lb, unsigned& ub) {
    bool firstTime = false;
    double ** cMatrix = new double * [dim];
    int mode;
    std::vector<int> vCycle;
    unsigned i, j, nChild = 0;

    Node *nodeAux;
    Node *nodeCurr = NULL;

    // primeira execucao
    if (nodes.empty()) {
        nodeCurr = new Node;
        nodeCurr->n = 0;
        firstTime = true;
        std::cout << "first time" << std::endl;
    } else {
        nodeCurr = nodes.front();
    }

    // copia matrix de custos para double, usado pelo hungaro
    for (i = 0; i < dim; i++) {
        cMatrix[i] = new double[dim];
        for (j = 0; j < dim; j++) {
            if (i==j) {
                // inviabiliza a -> a
                cMatrix[i][j] = std::numeric_limits<int>::max();
            } else {
                cMatrix[i][j] = matrix[i][j];
            }
        }
    }

    // inviabiliza o uso dos arcos proibidos
    for (i = 0; i < nodeCurr->prohibited.size(); i++) {
        cMatrix[nodeCurr->prohibited[i].first]
            [nodeCurr->prohibited[i].second] = std::numeric_limits<int>::max();
    }

    // chamada ao hungaro para resolucao do problema de assignment
    hungarian_problem_t p;
    mode = HUNGARIAN_MODE_MINIMIZE_COST;
    hungarian_init(&p, cMatrix, dim, dim, mode);
    nodeCurr->solution = hungarian_solve(&p);
    // converte matrix binaria da solucao para rota (possivelmente c/ subciclo)
    nodeCurr->route = getVectorSolution(p.assignment, dim);
    // verifica ciclo e armazena como vetor de pares (arcos a serem proibidos)
    verifyCycle(nodeCurr->route, vCycle, nodeCurr->arrows, dim);

    // finaliza hungaro
    hungarian_free(&p);

    // primeira solucao do assignment 
    if (firstTime) {
        // assignment equivalente a um tour completo = solucao viavel ao TSP
        if (nodeCurr->route.size() == dim + 1 && vCycle.size() == dim + 1) {
            ub = nodeCurr->solution;
            std::cout << "Solucao otima encontrada! " << std::endl;
            tsp::printVector<int>(nodeCurr->route);
        }
        // define lower bound
        lb = nodeCurr->solution;
    }

    // branch-and-bound
    if (nodeCurr->route.size() > dim + 1 && nodeCurr->solution < ub) {
        // cada novo no eh uma copia do no atual adicionado de um dos arcos como proibidos
        for (i = 0; i < nodeCurr->arrows.size(); i++) {
            nodeAux = new Node;
            nodeAux->n = (int) nodes.size() + i;
            nodeAux->solution = nodeCurr->solution;
            nodeAux->route = nodeCurr->route;
            nodeAux->arrows = nodeCurr->arrows;
            nodeAux->prohibited = nodeCurr->prohibited;
            nodeAux->prohibited.push_back(nodeCurr->arrows[i]);
            nodes.push_back(nodeAux);
            nChild++;
        }
        //std::cout << "branching de fator " << nChild << std::endl;
    } else if (nodeCurr->route.size() == dim + 1 && nodeCurr->solution >= lb && nodeCurr->solution < ub) {
        std::cout << "novo best" << std::endl;
        best.swap(nodeCurr->route);
        ub = nodeCurr->solution;
        int s = 0;
        std::vector<Node *>::iterator prev, it = nodes.begin();
        it++;
        while (it != nodes.end()) {
            if ((*it)->solution >= ub) {
                delete *it;
                nodes.erase(it);
                it--;
                s++;
            }
            it++;
        }
        std::cout << s << " nos eliminados por bound, |nodes| = " << nodes.size() << std::endl;
        tsp::printVector<int>(nodeCurr->route);
    }

    if(!firstTime) {
        delete *(nodes.begin());
        nodes.erase(nodes.begin());
    }

    // desaloca memoria
    for (i = 0; i < dim; i++) {
        delete[] cMatrix[i];
    }
    delete[] cMatrix;

}


int main(int argc, char *argv[]) {
    std::string file = "";
    std::stringstream ss;
    std::vector<Node *> nodes;
    unsigned lb, ub;

    //Node *root = new Node;

    if (argc >= 2) {
        if (argc >= 3) {
            ss.str(argv[2]);
            ss >> lb;
        } else {
            lb = 0;
        }
        file.assign(argv[1]);
    } else {
        std::cout << "Uso " << argv[0] << " arquivo_entrada [lb] " << std::endl;
        return 1;
    }

    srand(time(NULL));
    Instance instance(file);

    try {
        instance.readInfo(tsp::TSP);
        instance.printInfo();
        instance.initMatrix();
        instance.readMatrixData();
    } catch (std::exception &e) {
        std::cout << e.what() << std::endl;
        std::cout << "Falha na inicializacao... encerrando." << std::endl;
        return 1;
    }

    Node aux;
    ub = 118361;// std::numeric_limits<int>::max();
    int i = 0;
    do {
        i++;
        bnb(nodes, instance.getMatrix(), instance.getDim(), lb, ub);
        std::cout << "Lista de nos a serem expandidos..." << std::endl;
        for(int k=0; k<nodes.size(); k++) {
            printNode(*nodes[k]);
        }
        if (i%500==0) {
            std::cout << "(" << i << ") LB= " << lb << ", UB=" << ub << std::endl;
            std::cout << "|nodes| = " << nodes.size() << std::endl;
        }
        usleep(1000000);
    } while (!nodes.empty());
    std::cout << 299;
}
