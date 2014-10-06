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
void initBranchAndBound(const int ** matrix, const unsigned dim, unsigned b) {
    std::vector<Node> nodes;
    std::vector<int> bestRoute;
    unsigned lb, ub = inf;

    unsigned x = 0;
    Node nodeCurr = (x) ? rootBBHung(matrix,dim) : rootBB1Tree(matrix,dim);

    if (isRootOptimal(nodeCurr, dim, lb, ub)) {
        std::cout << "Solucao otima encontrada! " << std::endl;
        tsp::printVector<int>(nodeCurr.route);
        return;
    } else {
        nodes.push_back(nodeCurr);
    }

    // comeca branch-and-bound
    bnb(bestRoute, nodes, matrix, dim, lb, ub, b, x);
    std::cout << "Rota (custo = " << tsp::cost(bestRoute, matrix) << ")" << std::endl;
    tsp::printVector<int>(bestRoute);
}

/**
 * Relaxacao pelo algoritmo 1-tree
 */
void oneTree(Node& node, int ** matrix, const unsigned dim) {
    bool ** sol1Tree  = newBoolMatrix(dim); // Total de arestas em cada vertice
    unsigned * degree = new unsigned[dim]; // Total de arestas em cada vertice
    unsigned fn, sn, k, z, i, cost = 0;

    prim1Tree(dim, matrix, degree, sol1Tree, cost);

    // selecao de duas arestas mais baratas
    fn = 1, sn = 2;
    if (matrix[0][fn] > matrix[0][sn]) {
        fn = 2;
        sn = 1;
    }

    for (i = 3; i < dim; i++) {
        if (matrix[0][i] < matrix[0][fn]) {
            sn = fn;
            fn = i;
        } else {
            if(matrix[0][i] < matrix[0][sn]) {
                sn = i;
            }
        }
    }

    degree[0] = 2;
    degree[fn]++;
    degree[sn]++;
    sol1Tree[0][fn] = sol1Tree[fn][0] = true;
    sol1Tree[0][sn] = sol1Tree[sn][0] = true; 
    cost += matrix[0][fn] + matrix[0][sn];

    node.cost = cost;

    if (isFeasible(dim, degree, k, z)) {
        node.route = get1TreeVectorSolution(sol1Tree, dim);
    } else {
        node.arrows.clear();
        for (i = 0; i < dim; i++) {
            if (sol1Tree[i][k]) {
                node.arrows.push_back(std::make_pair<int,int>(i,k));
            }
        }
    }

    delete[] degree;
    free<bool>(sol1Tree, dim);
}

/**
 * Relaxacao pelo algoritmo hungaro para problema do assignment
 */
void hungarian(Node& nodeCurr, double ** matrix, const unsigned dim) {
    // chamada ao hungaro para resolucao do problema de assignment
    hungarian_problem_t p;
    int mode = HUNGARIAN_MODE_MINIMIZE_COST;
    hungarian_init(&p, matrix, dim, dim, mode);
    nodeCurr.cost = hungarian_solve(&p);
    // converte matrix binaria da solucao para rota (possivelmente c/ subciclo)
    nodeCurr.route = getVectorSolution<int>(p.assignment, dim);
    // verifica ciclo e armazena como vetor de pares (arcos a serem proibidos)
    verifyCycle(nodeCurr.route, nodeCurr.arrows, dim);
    // finaliza hungaro
    hungarian_free(&p);
}

/**
 * Branch-and-bound
 */
void bnb(std::vector<int>& bestRoute, std::vector<Node>& nodes, const int ** matrix, const unsigned dim, unsigned& lb, unsigned& ub, unsigned b, unsigned x) {
    Node nodeAux;
    double ** dMatrix = NULL, d = 0;
    int    ** cMatrix = NULL, c = 0;
    unsigned i, nBound;
    unsigned long count = 0;
    std::string strat = "";

    // referencia para funcao/estrategia usada
    Node (* curr)(std::vector<Node>&) = NULL;

    switch(b) { 
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
        if (count % 7500 == 0) {
            doLog(ub, lb, nodes.size(), count, strat);
        }

        if (x) {
            // hungaro
            dMatrix = copyMatrixFromTo<int,double>(const_cast<const int **>(matrix), dim);
        } else {
            // 1-tree
            cMatrix = copyMatrixFromTo<int, int>(const_cast<const int **>(matrix), dim);
        }

        // seleciona (e retira da lista) no de acordo com estrategia
        Node nodeCurr = curr(nodes);

        // atualiza lower bound
        if (nodeCurr.cost > lb) {
            updateLB(nodes, lb);
        }

        // inviabiliza o uso dos arcos proibidos
        for (i = 0; i < nodeCurr.prohibited.size(); i++) {
            if (x) { // hungaro
                dMatrix[nodeCurr.prohibited[i].first][nodeCurr.prohibited[i].second] = inf;
            } else { // 1-tree
                cMatrix[nodeCurr.prohibited[i].first][nodeCurr.prohibited[i].second] = inf;
                cMatrix[nodeCurr.prohibited[i].second][nodeCurr.prohibited[i].first] = inf;
            }
        }

        // cada novo no eh uma copia do no atual adicionado de um dos arcos como proibidos
        for (i = 0; i < nodeCurr.arrows.size(); i++) {
            // branch-and-bound
            nodeAux.n = (int) nodes.size() + 1;
            nodeAux.prohibited = nodeCurr.prohibited;
            nodeAux.prohibited.push_back(nodeCurr.arrows[i]);
            nodeAux.route.clear();

            if (x) { // hungaro
                d = dMatrix[nodeCurr.arrows[i].first][nodeCurr.arrows[i].second];
                dMatrix[nodeCurr.arrows[i].first][nodeCurr.arrows[i].second] = inf;
                hungarian(nodeAux, dMatrix, dim);
            } else { // 1-tree
                c = cMatrix[nodeCurr.arrows[i].first][nodeCurr.arrows[i].second];
                cMatrix[nodeCurr.arrows[i].first][nodeCurr.arrows[i].second] = inf;
                cMatrix[nodeCurr.arrows[i].second][nodeCurr.arrows[i].first] = inf;
                oneTree(nodeAux, cMatrix, dim);
            }

            // se nao eh uma solucao viavel TSP (ciclo hamiltoniano) mas custo esta abaixo do UB
            if (!isValidCH(nodeAux.route, dim) && nodeAux.cost < ub) {
                nodes.push_back(nodeAux);
                count++;
            } else { // solucao viavel e de menor custo, novo UB
                if (isNewUB(nodeAux, dim, ub)) {
                    ub = nodeAux.cost;
                    bestRoute = nodeAux.route;

                    nBound = 0;
                    std::vector<Node>::iterator prev, it;
                    for (it = nodes.begin(); it != nodes.end(); ++it) {
                        if (it->cost >= ub) {
                            nodes.erase(it);
                            nBound++;
                            it--;
                        }
                    }
                    doLog(ub, lb, nodes.size(), count, strat);
                    std::cout << nBound << " nos eliminados por bound" << std::endl;
                }
            }

            if (x) { // hungaro 
                dMatrix[nodeCurr.arrows[i].first][nodeCurr.arrows[i].second] = d;
            } else { // 1-tree
                cMatrix[nodeCurr.arrows[i].first][nodeCurr.arrows[i].second] = c;
                cMatrix[nodeCurr.arrows[i].second][nodeCurr.arrows[i].first] = c;
            }
        }

        if (x)
            free<double>(dMatrix, dim);
        else
            free<int>(cMatrix, dim);
    }
    std::cout << "FIM Branch-and-bound - LB=" << lb << " UB=" << ub << std::endl;
}

Node rootBB1Tree(const int ** matrix, const unsigned dim) {
    Node nodeCurr;

    nodeCurr.n = 0;
    // one tree
    int ** cMatrix = copyMatrixFromTo<int, int>(matrix, dim);
    oneTree(nodeCurr, cMatrix, dim);
    free<int>(cMatrix, dim);

    return nodeCurr;
}

Node rootBBHung(const int ** matrix, const unsigned dim) {
    Node nodeCurr;

    nodeCurr.n = 0;
    //hungaro
    double ** cMatrix = copyMatrixFromTo<int, double>(matrix, dim);
    hungarian(nodeCurr, cMatrix, dim);
    free<double>(cMatrix, dim);

    return nodeCurr;
}

bool isFeasible(const unsigned dim, unsigned * degree, unsigned& k, unsigned& z) {
    unsigned i, max = 0, min = inf;
    bool f = true;
    for (i = 0; i < dim; i++) {
        f = (degree[i]==2) ? f : false;

        if (degree[i] > max) {
            max = degree[i];
            k = i;
        }

        if (degree[i] < min) {
            min = degree[i];
            z = i;
        }
    }
    return f;
}

/**
 * Converte matrix binaria do algoritmo hungaro em uma rota
 */
template<typename type> std::vector<int> getVectorSolution(type ** matrix, int dim) {
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
            if (matrix[i][j]) {
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
 * Converte matrix binaria do algoritmo hungaro em uma rota
 */
std::vector<int> get1TreeVectorSolution(bool ** matrix, unsigned dim) {
    unsigned i, j;
    std::vector<int> route;
    bool visited[dim];

    //    tsp::printMatrix(const_cast<const bool **>(matrix), dim, 3);
    // partindo da primeira cidade

    for(i=0; i<dim; i++)
        visited[i] = false;

    visited[0] = true;
    i = 0;
    route.push_back(0);
    while (route.size() < dim) {
        for (j = 0; j < dim; j++) {
            if (!visited[j] && matrix[i][j]) {
                route.push_back(j);
                visited[j] = true;
                i = j;
                break;
            } 
        }
    }
    route.push_back(0);
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

// ##################
// Funcoes auxiliares
// ##################

bool ** newBoolMatrix(const unsigned dim) {
    bool ** bMatrix = new bool * [dim];
    unsigned i, j;
    for (i = 0; i < dim; i++) {
        bMatrix[i] = new bool[dim];
        for (j = 0; j < dim; j++) {
            bMatrix[i][j] = false; 
        }
    }
    return bMatrix;
}

template<typename from, typename to> to ** copyMatrixFromTo(const from ** matrix, const unsigned dim, const unsigned s = 0) {
    to ** dMatrix = new to * [dim-s];
    unsigned i, j;
    for (i = 0; i < dim - s; i++) {
        dMatrix[i] = new to[dim-s];
        for (j = 0; j < dim - s; j++) {
            if (i==j)
                dMatrix[i][j] = (to) inf;
            else
                dMatrix[i][j] = (to) matrix[i+s][j+s];
        }
    }
    return dMatrix;
}

template<typename type> void free(type ** matrix, const unsigned dim) {
    unsigned i;
    if (matrix != NULL) {
        for (i = 0; i < dim; i++) {
            if (matrix[i] != NULL)
                delete[] matrix[i];
        }
        delete[] matrix;
    }
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

void printDegrees(unsigned * degree, const unsigned dim) {
    for (unsigned i = 0; i < dim; i++) {
        std::cout << std::setw(3) << i; 
    }
    std::cout << std::endl;
    for (unsigned i = 0; i < dim; i++) {
        std::cout << std::setw(3) << degree[i]; 
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
    unsigned bb = 0, i, indexBB = 0;
    for (i = 0; i < nodes.size(); i++) {
        if (nodes[i].cost > bb) {
            bb = nodes[i].cost;
            indexBB = i;
        }
    }
    Node node = nodes[indexBB];
    nodes.erase(nodes.begin() + indexBB);
    return node;
}
