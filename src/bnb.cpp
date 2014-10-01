/*
 * @file bnb.cpp
 *      
 * @author: Fabio Cruz B. de Albuquerque
 * 
 * @date: 01/10/14
 */

#include "../include/headers.h"

/**
 * Avalia solucao inicial de assignment (raiz do b&b)
 * se nao for otimo para o TSP, inicia B&B
 */
void initBranchAndBound(const int ** matrix, const unsigned dim) {
    Node nodeCurr;
    double ** dMatrix = copyMatrix2Double(matrix, dim);
    std::vector<Node> nodes;
    std::vector<int> bestRoute;
    unsigned lb, ub = inf;

    // raiz
    nodeCurr.n = 0;
    hungarian(nodeCurr, dMatrix, dim);
    if (isRootOptimal(nodeCurr, dim, lb, ub)) {
        std::cout << "Solucao otima encontrada! " << std::endl;
        tsp::printVector<int>(nodeCurr.route);
        return;
    } else {
        nodes.push_back(nodeCurr);
    }
    free<double>(dMatrix, dim);

    // comeca branch-and-bound
    bnb(bestRoute, nodes, matrix, dim, lb, ub);
    std::cout << "Rota (custo = " << tsp::cost(bestRoute, matrix) << ")" << std::endl;
    tsp::printVector<int>(bestRoute);
}

/**
 * Algoritmo hungaro para problema do assignment
 */
void hungarian(Node& nodeCurr, double ** matrix, const int dim) {
    int mode;

    // chamada ao hungaro para resolucao do problema de assignment
    hungarian_problem_t p;
    mode = HUNGARIAN_MODE_MINIMIZE_COST;
    hungarian_init(&p, matrix, dim, dim, mode);
    nodeCurr.cost = hungarian_solve(&p);
    // converte matrix binaria da solucao para rota (possivelmente c/ subciclo)
    nodeCurr.route = getVectorSolution(p.assignment, dim);
    // verifica ciclo e armazena como vetor de pares (arcos a serem proibidos)
    verifyCycle(nodeCurr.route, nodeCurr.arrows, dim);
    // finaliza hungaro
    hungarian_free(&p);

}

/**
 * Converte matrix binaria do algoritmo hungaro em uma rota
 */
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

/**
 * Verifica existencia de subciclos em uma rota
 */
void verifyCycle(std::vector<int> &sol, std::vector< std::pair<int,int> > &cycleArrows, unsigned dim) {
    bool cDetected = false;
    long cSize = -1, cBegin = 0, cBeginTemp = 0, i;
    std::vector<int> cycle;

    // busca menor ciclo na rota
    if (sol.size() > dim + 1) {
        for (i = 1; i < (int) sol.size(); i++) {
            if (cDetected) {
                cBeginTemp = i;
                cDetected = false;
            } else {
                if (sol[i] == sol[cBeginTemp]) {
                    cDetected = true;
                    // se encontrou ciclo menor ou primeiro ciclo
                    if ((i - cBeginTemp) < cSize || (cSize == -1)) {
                        cBegin = cBeginTemp;
                        cSize = i - cBeginTemp;
                    }
                }
            }
        }
    }

    cycle.clear();
    cycleArrows.clear();

    // intervalo da rota com ciclo
    for (i = cBegin; i <= cBegin + cSize; i++) {
        cycle.push_back(sol[i]);
    }

    // clientes/cidades pertencentes ao intervalo
    // e respectivos arcos
    if (sol.size() - dim > 1) {
        for (i = 0; i < (int) cycle.size() - 1; i++) {
            cycleArrows.push_back(std::make_pair(cycle[i], cycle[i+1]));
        }
    }
}

/**
 * Branch-and-bound
 */
void bnb(std::vector<int>& bestRoute, std::vector<Node>& nodes, const int ** matrix, unsigned dim, unsigned& lb, unsigned& ub) {
    Node nodeAux;
    double ** dMatrix, t;
    unsigned i, nBound;
    unsigned long c = 0;
    std::string strat = "";

    // referencia para funcao/estrategia usada
    Node (* curr)(std::vector<Node>&) = NULL;
    int s = rand() % 2; 
    switch(s) { 
        case 0:
            curr = &dfs;
            strat = "DFS";
            break;
        case 1:
            curr = &bfs;
            strat = "BFS";
            break;
        case 2:
            curr = &bestb;
            strat = "BEST-BOUND";
            break;
    }
    std::cout << std::endl;

    while (!nodes.empty()) {

        if (c % 2000 == 0)
            std::cout << "UB: " << ub << std::setw(4) << " LB: " << lb << std::setw(4) 
                << " GAP: " <<  gap(lb, ub) << "%" << std::setw(4) << " Numero de nos abertos: " << nodes.size() 
                << std::setw(4) << " Iteracao: " << c << std::setw(4) << strat << std::endl;

        dMatrix = copyMatrix2Double(matrix, dim);

        // seleciona (e retira da lista) no de acordo com estrategia
        Node nodeCurr = curr(nodes);

        // inviabiliza o uso dos arcos proibidos
        for (i = 0; i < nodeCurr.prohibited.size(); i++) {
            dMatrix[nodeCurr.prohibited[i].first]
                [nodeCurr.prohibited[i].second] = inf; //std::numeric_limits<int>::max();
        }

        // cada novo no eh uma copia do no atual adicionado de um dos arcos como proibidos
        for (i = 0; i < nodeCurr.arrows.size(); i++) {
            c++;
            // branch-and-bound
            nodeAux.n = (int) nodes.size() + 1;
            nodeAux.prohibited = nodeCurr.prohibited;
            nodeAux.prohibited.push_back(nodeCurr.arrows[i]);

            t = dMatrix[nodeCurr.arrows[i].first][nodeCurr.arrows[i].second];
            dMatrix[nodeCurr.arrows[i].first][nodeCurr.arrows[i].second] = inf;// std::numeric_limits<int>::max();

            hungarian(nodeAux, dMatrix, dim);

            // se nao eh uma solucao viavel TSP (ciclo hamiltoniano) mas custo esta abaixo do UB
            if (!isValidCH(nodeAux.route, dim) && nodeAux.cost < ub) {
                nodes.push_back(nodeAux);
            } else { 
                // solucao viavel e de menor custo, novo UB
                if (isNewUB(nodeAux, dim, ub)) {
                    ub = nodeAux.cost;
                    bestRoute = nodeAux.route;

                    std::cout << "Novo valor Upper Bound: " << ub; 
                    nBound = 0;
                    std::vector<Node>::iterator prev, it = nodes.begin();
                    for (it != nodes.begin(); it != nodes.end(); ++it) {
                        if (it->cost >= ub) {
                            nodes.erase(it);
                            nBound++;
                            it--;
                        }
                    }
                    std::cout << std::setw(4) << nBound << " nos eliminados por bound" << std::endl;
                    //tsp::printVector<int>(nodeAux.route);
                }
            }
            dMatrix[nodeCurr.arrows[i].first][nodeCurr.arrows[i].second] = t;
        }

        // atualiza lower bound
        lb = ub;
        for (i = 0; i < nodes.size(); i++) {
            if (nodes[i].cost < lb)
                lb = nodes[i].cost;
        }
        free<double>(dMatrix, dim);
        c++;
    }
    std::cout << "FIM Branch-and-bound - LB=" << lb << " UB=" << ub << std::endl;
}

// ##################
// Funcoes auxiliares
// ##################

double ** copyMatrix2Double(const int ** matrix, const unsigned dim) {
    double ** dMatrix = new double * [dim];
    unsigned i, j;
    for (i = 0; i < dim; i++) {
        dMatrix[i] = new double[dim];
        for (j = 0; j < dim; j++) {
            if (i==j)
                dMatrix[i][j] = inf;//std::numeric_limits<int>::max();
            else
                dMatrix[i][j] = matrix[i][j];
        }
    }
    return dMatrix;
}

template<typename type> void free(type ** matrix, const unsigned dim) {
    unsigned i;
    for (i = 0; i < dim; i++) {
        delete[] matrix[i];
    }
    delete[] matrix;
}


void printNode(const Node& node) {
    std::cout << std::endl << "Node " << &node << " " << node.n << ", cost: " << node.cost << std::endl;
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

/**
 * Depth-first
 */
Node dfs(std::vector<Node>& nodes) {
    Node node = nodes.back();
    nodes.pop_back();
    return node;
}

/**
 * Breadth-first
 */
Node bfs(std::vector<Node>& nodes) {
    Node node = nodes.front();
    nodes.erase(nodes.begin());
    return node;
}

/**
 * Best-first
 */
Node bestb(std::vector<Node>& nodes) {
    unsigned bb = inf, i, indexBB = 0;
    for (i = 0; i < nodes.size(); i++) {
        if (nodes[i].cost < bb) {
            bb = nodes[i].cost;
            indexBB = i;
        }
    }
    Node node = nodes[indexBB];
    nodes.erase(nodes.begin() + indexBB);
    return node;
}
