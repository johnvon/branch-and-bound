/*
 * @file bnb.cpp
 *      
 * @author: Fabio Cruz B. de Albuquerque
 * 
 * @date: 01/10/14
 */

#include "../include/bnb.h"

/**
 * Avalia solucao inicial de assignment (raiz do b&b)
 * se nao for otimo para o TSP, inicia B&B
 */
void initBranchAndBound(const int ** matrix, const unsigned dim, unsigned b, unsigned x) {
    Node nodeCurr;
    std::list<Node> nodes;
    std::vector<int> bestRoute;
    unsigned lb, ub;

    // heuristica ILS-RVND
    Construct c(matrix, dim);
    LocalSearch l(matrix, dim);

    // upper-bound inicia
    ub = l.ilsRvnd(c, 0, 10) + 1;

    switch(x) {
        case 0: // Hungaro
            nodeCurr = rootBBHung(matrix,dim);
            if (isRootOptimal(nodeCurr, dim, lb, ub)) {
                bestRoute = nodeCurr.route;
            } else {
                nodes.push_back(nodeCurr);
            }
            bnbHung(bestRoute, nodes, matrix, dim, lb, ub, b);
            break;
        case 1: // 1-Tree
            nodeCurr = rootBB1Tree(matrix,dim);
            if (isRootOptimal(nodeCurr, dim, lb, ub)) {
                bestRoute = nodeCurr.route;
            } else {
                nodes.push_back(nodeCurr);
            }
            bnb1Tree(bestRoute, nodes, matrix, dim, lb, ub, b);
            break;
        case 2: // Lagrangeana
            nodeCurr = rootBBLR(matrix, dim, ub);
            lb = nodeCurr.cost;
            std::cout << nodeCurr.cost << std::endl;
            if (!nodeCurr.route.empty()) {
                ub = nodeCurr.cost = cost(nodeCurr.route, matrix);
                bestRoute = nodeCurr.route;
                std::cout << "Solucao otima encontrada na raiz com custo: " << ub << std::endl;
            } else if (std::ceil(nodeCurr.cost) + 1 == ub) {
                std::cout << "Solucao fracionaria ~= UB" << std::endl;
                bestRoute = l.getRoute();
                bestRoute.push_back(bestRoute.front());
            } else {
                nodes.push_back(nodeCurr);
            }
            bnbLR(bestRoute, nodes, matrix, dim, lb, ub, b);
            break;
    }

    std::cout << "Rota (custo = " << tsp::cost(bestRoute, matrix) << ")" << std::endl;
    tsp::printVector<int>(bestRoute);
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

Node rootBB1Tree(const int ** matrix, const unsigned dim) {
    Node nodeCurr;
    nodeCurr.n = 0;
    int ** cMatrix   = copyMatrixFromTo<int, int>(matrix, dim);
    bool ** sol1Tree = newBoolMatrix(dim); // Total de arestas em cada vertice
    unsigned * degree = new unsigned[dim];

    oneTree<int>(nodeCurr, cMatrix, dim, sol1Tree, degree);

    delete[] degree;
    free<int>(cMatrix, dim);
    free<bool>(sol1Tree, dim);

    return nodeCurr;
}

Node rootBBLR(const int ** matrix, const unsigned dim, unsigned ub) {
    double ** dMatrix = copyMatrixFromTo<int,double>(const_cast<const int **>(matrix), dim);
    bool   ** sol1Tree = newBoolMatrix(dim);

    Node nodeCurr;
    nodeCurr.n = 0;
    for (unsigned i = 0; i < dim; i++) {
        nodeCurr.u.push_back(0);
    }

    lagrangean(nodeCurr, matrix, dMatrix, sol1Tree, dim, ub, true);

    free<bool>(sol1Tree, dim);
    free<double>(dMatrix, dim);

    return nodeCurr;
}

bool isRootOptimal(Node& root, const unsigned dim, unsigned& lb, unsigned& ub) {
    // define lower bound
    lb = root.cost;
    // assignment equivalente a um tour completo = solucao viavel ao TSP
    if (isValidCH(root.route, dim)) {
        ub = root.cost;
        return true;
    }
    return false;
}

void bnbHung(std::vector<int>& bestRoute, std::list<Node>& nodes, const int ** matrix, const unsigned dim, unsigned& lb, unsigned& ub, unsigned b) {
    Node nodeAux;
    bool   ** sol1Tree = newBoolMatrix(dim);
    double ** dMatrix  = NULL, d = 0;
    time_t t0, t1, t2;
    unsigned i, nBound;
    unsigned long count = 0;
    std::list<Node>::iterator prev, it, itCurr;

    std::string strat = "";
    std::list<Node>::iterator (* curr)(std::list<Node>&) = NULL;

    switch(b) {
        case 0:
            curr = &dfs;
            strat = "Hung-DFS";
            break;
        case 1:
            curr = &bfs;
            strat = "Hung-BFS";
            break;
        case 2:
            curr = &bestb;
            strat = "Hung-BEST-BOUND";
            break;
        case 3:
            curr = &randb;
            strat = "Hung-RAND-BOUND";
            break;
    }

    t0 = t1 = t2 = time(NULL); 
    doLog(ub, lb, nodes.size(), count, strat);

    while (!nodes.empty()) {
        if (difftime(time(NULL),t1) > 1) {
            time(&t1);
            doLog(ub, lb, nodes.size(), count, strat);
        }

        dMatrix = copyMatrixFromTo<int,double>(const_cast<const int **>(matrix), dim);
        // seleciona (e retira da lista) no de acordo com estrategia
        itCurr = curr(nodes);

        // atualiza lower bound
        if (itCurr->cost > lb) {
            updateLB(nodes, lb);
        }
        // inviabiliza o uso dos arcos proibidos
        for (i = 0; i < itCurr->prohibited.size(); i++) 
            dMatrix[itCurr->prohibited[i].first][itCurr->prohibited[i].second] = inf;

        // cada novo no eh uma copia do no atual adicionado de um dos arcos como proibidos
        for (i = 0; i < itCurr->arrows.size(); i++) {
            // branch-and-bound
            nodeAux.n = (int) nodes.size() + 1;
            nodeAux.prohibited = itCurr->prohibited;
            nodeAux.prohibited.push_back(itCurr->arrows[i]);
            nodeAux.route.clear();

            d = dMatrix[itCurr->arrows[i].first][itCurr->arrows[i].second];
            dMatrix[itCurr->arrows[i].first][itCurr->arrows[i].second] = inf;
            hungarian(nodeAux, dMatrix, dim);

            // se nao eh uma solucao viavel TSP (ciclo hamiltoniano) mas custo esta abaixo do UB
            if (!isValidCH(nodeAux.route, dim) && nodeAux.cost < ub) {
                nodes.push_back(nodeAux);
                count++;
            } else { // solucao viavel e de menor custo, novo UB
                if (isNewUB(nodeAux, dim, ub)) {
                    ub = nodeAux.cost;
                    bestRoute = nodeAux.route;

                    nBound = 0;
                    for (it = nodes.begin(); it != nodes.end();) {
                        if (it->cost > ub && it != itCurr) {
                            it = nodes.erase(it);
                            nBound++;
                        } else {
                            ++it;
                        }
                    }
                }
            }
            dMatrix[itCurr->arrows[i].first][itCurr->arrows[i].second] = d;
        }

        free<double>(dMatrix, dim);
        if (difftime(time(NULL), t2) > 30) { 
            time(&t2);
            std::cout << "---------------------------------------" << std::endl;
            std::cout << "# Tempo gasto: " << difftime(t2,t0) << "s" << std::endl;
            std::cout << "---------------------------------------" << std::endl;
        }
        nodes.erase(itCurr);
    }
    std::cout << "FIM Branch-and-bound - LB=" << lb << " UB=" << ub << std::endl;
    free<bool>(sol1Tree, dim);
}

void bnb1Tree(std::vector<int>& bestRoute, std::list<Node>& nodes, const int ** matrix,
        const unsigned dim, unsigned& lb, unsigned& ub, unsigned b) {
    Node nodeAux;
    bool   ** sol1Tree = newBoolMatrix(dim);
    int    ** cMatrix  = NULL, c = 0;
    time_t t0, t1, t2;
    unsigned i, j, nBound, * degree = new unsigned[dim];
    unsigned long count = 0;

    std::list<Node>::iterator prev, it, itCurr;

    std::string strat = "";
    std::list<Node>::iterator (* curr)(std::list<Node>&) = NULL;

    switch(b) {
        case 0:
            curr = &dfs;
            strat = "1-Tree-DFS";
            break;
        case 1:
            curr = &bfs;
            strat = "1-Tree-BFS";
            break;
        case 2:
            curr = &bestb;
            strat = "1-Tree-BEST-BOUND";
            break;
        case 3:
            curr = &randb;
            strat = "1-Tree-RAND-BOUND";
            break;
    }

    t0 = t1 = t2 = time(NULL); 
    doLog(ub, lb, nodes.size(), count, strat);

    while (!nodes.empty()) {
        if (difftime(time(NULL),t1) > 1) {
            time(&t1);
            doLog(ub, lb, nodes.size(), count, strat);
        }

        cMatrix = copyMatrixFromTo<int, int>(const_cast<const int **>(matrix), dim);

        // seleciona (e retira da lista) no de acordo com estrategia
        itCurr = curr(nodes);

        // atualiza lower bound
        if (itCurr->cost > lb) {
            updateLB(nodes, lb);
        }

        // fortalece os bounds
        for (i = 0; i < dim; i++) {
            for (j = 0; j < dim; j++) {
                if (cMatrix[i][j] != inf)
                    cMatrix[i][j] += itCurr->pi[i] + itCurr->pi[j];
            }
        }

        // inviabiliza o uso dos arcos proibidos
        for (i = 0; i < itCurr->prohibited.size(); i++) {
            cMatrix[itCurr->prohibited[i].first][itCurr->prohibited[i].second] = inf;
            cMatrix[itCurr->prohibited[i].second][itCurr->prohibited[i].first] = inf;
        }


        // cada novo no eh uma copia do no atual adicionado de um dos arcos como proibidos
        for (i = 0; i < itCurr->arrows.size(); i++) {
            // branch-and-bound
            nodeAux.n = (int) nodes.size() + 1;
            nodeAux.prohibited = itCurr->prohibited;
            nodeAux.prohibited.push_back(itCurr->arrows[i]);
            nodeAux.route.clear();

            c = cMatrix[itCurr->arrows[i].first][itCurr->arrows[i].second];
            cMatrix[itCurr->arrows[i].first][itCurr->arrows[i].second] = inf;
            cMatrix[itCurr->arrows[i].second][itCurr->arrows[i].first] = inf;
            oneTree<int>(nodeAux, cMatrix, dim, sol1Tree, degree);

            // se nao eh uma solucao viavel TSP (ciclo hamiltoniano) mas custo esta abaixo do UB
            if (!isValidCH(nodeAux.route, dim) && nodeAux.cost < ub) {
                nodes.push_back(nodeAux);
                count++;
            } else { // solucao viavel e de menor custo, novo UB
                if (isNewUB(nodeAux, dim, ub)) {
                    ub = nodeAux.cost;
                    bestRoute = nodeAux.route;

                    nBound = 0;
                    for (it = nodes.begin(); it != nodes.end();) {
                        if (it->cost > ub && it != itCurr) {
                            it = nodes.erase(it);
                            nBound++;
                        } else {
                            ++it;
                        }
                    }
                }
            }
            cMatrix[itCurr->arrows[i].first][itCurr->arrows[i].second] = c;
            cMatrix[itCurr->arrows[i].second][itCurr->arrows[i].first] = c;
        }
        free<int>(cMatrix, dim);

        if (difftime(time(NULL), t2) > 30) { 
            time(&t2);
            std::cout << "---------------------------------------" << std::endl;
            std::cout << "# Tempo gasto: " << difftime(t2,t0) << "s" << std::endl;
            std::cout << "---------------------------------------" << std::endl;
        }
        nodes.erase(itCurr);
    }
    std::cout << "FIM Branch-and-bound - LB=" << lb << " UB=" << ub << std::endl;
    free<bool>(sol1Tree, dim);
}

void bnbLR(std::vector<int>& bestRoute, std::list<Node>& nodes, const int ** matrix,
        const unsigned dim, unsigned& lb, unsigned& ub, unsigned b) {
    Node nodeAux;
    bool     ** sol1Tree = newBoolMatrix(dim);
    int      ** cMatrix  = NULL, c = 0;
    double   ** lMatrix  = NULL;
    time_t t0, t1, t2;
    unsigned i, nBound;
    unsigned long count = 0;
    std::list<Node>::iterator prev, it, itCurr;

    std::string strat = "";
    std::list<Node>::iterator (* curr)(std::list<Node>&) = NULL;

    switch(b) {
        case 0:
            curr = &dfs;
            strat = "LR-DFS";
            break;
        case 1:
            curr = &bfs;
            strat = "LR-BFS";
            break;
        case 2:
            curr = &bestb;
            strat = "LR-BEST-BOUND";
            break;
        case 3:
            curr = &randb;
            strat = "LR-RAND-BOUND";
            break;
    }

    t0 = t1 = t2 = time(NULL); 
    doLog(ub, lb, nodes.size(), count, strat);

    while (!nodes.empty()) {
        if (difftime(time(NULL),t1) > 1) {
            time(&t1);
            doLog(ub, lb, nodes.size(), count, strat);
        }

        cMatrix = copyMatrixFromTo<int,int>(const_cast<const int **>(matrix), dim);
        // seleciona (e retira da lista) no de acordo com estrategia
        itCurr = curr(nodes);

        // atualiza lower bound
        if (itCurr->cost > lb) {
            updateLB(nodes, lb);
        }

        // inviabiliza o uso dos arcos proibidos
        for (i = 0; i < itCurr->prohibited.size(); i++) {
            cMatrix[itCurr->prohibited[i].first][itCurr->prohibited[i].second] = inf;
            cMatrix[itCurr->prohibited[i].second][itCurr->prohibited[i].first] = inf;
        }

        // cada novo no eh uma copia do no atual adicionado de um dos arcos como proibidos
        for (i = 0; i < itCurr->arrows.size(); i++) {
            // branch-and-bound
            nodeAux.n = (int) nodes.size() + 1;
            nodeAux.prohibited = itCurr->prohibited;
            nodeAux.prohibited.push_back(itCurr->arrows[i]);
            nodeAux.route.clear();
            nodeAux.u = itCurr->u;

            c = cMatrix[itCurr->arrows[i].first][itCurr->arrows[i].second];
            cMatrix[itCurr->arrows[i].first][itCurr->arrows[i].second] = inf;
            cMatrix[itCurr->arrows[i].second][itCurr->arrows[i].first] = inf;

            lMatrix = copyMatrixFromTo<int,double>(const_cast<const int **>(cMatrix), dim);
            lagrangean(nodeAux, const_cast<const int **>(cMatrix), lMatrix, sol1Tree, dim, ub);
            free<double>(lMatrix, dim);

            // se nao eh uma solucao viavel TSP (ciclo hamiltoniano) mas custo esta abaixo do UB
            if (!isValidCH(nodeAux.route, dim) && nodeAux.cost < ub) {
                nodes.push_back(nodeAux);
                count++;
            } else { // solucao viavel e de menor custo, novo UB
                if (isValidCH(nodeAux.route, dim)) {
                    nodeAux.cost = cost(nodeAux.route, matrix);
                    if (nodeAux.cost <= ub) {
                        ub = nodeAux.cost;
                        bestRoute = nodeAux.route;

                        nBound = 0;
                        for (it = nodes.begin(); it != nodes.end();) {
                            if (it->cost > ub && it != itCurr) {
                                it = nodes.erase(it);
                                nBound++;
                            } else {
                                ++it;
                            }
                        }
                    }
                }       
            }

            cMatrix[itCurr->arrows[i].first][itCurr->arrows[i].second] = c;
            cMatrix[itCurr->arrows[i].second][itCurr->arrows[i].first] = c;
        }

        free<int>(cMatrix, dim);

        if (difftime(time(NULL), t2) > 30) { 
            time(&t2);
            std::cout << "---------------------------------------" << std::endl;
            std::cout << "# Tempo gasto: " << difftime(t2,t0) << "s" << std::endl;
            std::cout << "---------------------------------------" << std::endl;
        }
        nodes.erase(itCurr);

    }
    std::cout << "FIM Branch-and-bound - LB=" << lb << " UB=" << ub << std::endl;
    free<bool>(sol1Tree, dim);

}

/**
 * Relaxacao pelo algoritmo hungaro para o problema do assignment
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
 * Relaxacao 1-tree
 */
template <typename type>
void oneTree(Node& node, type ** matrix, const unsigned dim, bool ** sol1Tree, unsigned * degree) {
    unsigned fn, sn, k, i;
    type cost = 0;

    resetBoolMatrix(sol1Tree, dim); // Total de arestas em cada vertice
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

    node.pi.clear();
    for (i = 0; i < dim; i++) {
        node.pi.push_back(degree[i] - 2);
    }

    if (isFeasible(dim, degree, k)) {
        node.route = get1TreeVectorSolution(sol1Tree, dim);
    } else {
        node.arrows.clear();
        for (i = 0; i < dim; i++) {
            if (sol1Tree[i][k]) {
                node.arrows.push_back(std::make_pair<int,int>(i,k));
            }
        }
    }
}

bool isFeasible(const unsigned dim, unsigned * degree, unsigned& k) {
    unsigned i, max = 0;
    bool f = true;
    for (i = 0; i < dim; i++) {
        f = (degree[i]==2) ? f : false;

        if (degree[i] > max) {
            max = degree[i];
            k = i;
        }
    }
    return f;
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

void lagrangean(Node& nodeCurr, const int ** originalMatrix, double ** cMatrix, bool ** sol1Tree, const unsigned dim, unsigned ub, bool root) {
    Node nodeTemp;
    double * subgradient = new double[dim], e = 1.0;
    unsigned * degree = new unsigned[dim], i, j, counter = 0;

    nodeTemp.cost = 0;
    while (e > 0.001) {
        // atualiza matrix com multiplicadores
        for (i = 0; i < dim; i++) {
            for (j = 0; j < dim; j++) {
                if (cMatrix[i][j] != inf)
                    cMatrix[i][j] = originalMatrix[i][j] - nodeCurr.u[i] - nodeCurr.u[j];
            }
        }

        // calcula 1-Tree
        oneTree<double>(nodeCurr, cMatrix, dim, sol1Tree, degree);

        // atualiza custo da funcao objetivo com multiplicadores
        nodeCurr.cost += 2 * sum(nodeCurr.u);

        if (root && std::ceil(nodeCurr.cost) == ub) {
            std::cout << "(break raiz)custo fracionario ≃ UB" << nodeCurr.cost << std::endl;
            break;
        }

        // solucao otima lagrange
        if (!nodeCurr.route.empty()) {
            break;
        }

        // calcula subgradient
        for (i = 0; i < dim; i++) {
            subgradient[i] = 2 - degree[i] * 1.0;
        }

        // atualiza multiplicadores
        for (i = 0; i < dim; i++)
            nodeCurr.u[i] += e * ((ub  - nodeCurr.cost ) / sumsqr(subgradient, dim)) * subgradient[i];

        // dual lagrange, max z(RLu)
        if (nodeCurr.cost > nodeTemp.cost)
            nodeTemp = nodeCurr, counter = 0;
        else
            counter++;

        // atualiza passo
        if (counter > 10)
            e *= 0.5, counter = 0;
    }

    // se nao for raiz e nao encontrou solucao viavel, usa o maior custo
    if (!root && nodeCurr.route.empty())
        nodeCurr = nodeTemp;

    delete[] subgradient;
    delete[] degree;
}

// ##################
// Funcoes auxiliares
// ##################

bool ** newBoolMatrix(const unsigned dim) {
    bool ** bMatrix = new bool * [dim];
    for (unsigned i = 0; i < dim; i++) {
        bMatrix[i] = new bool[dim];
        for (unsigned j = 0; j < dim; j++) {
            bMatrix[i][j] = false; 
        }
    }
    return bMatrix;
}

void resetBoolMatrix(bool ** bMatrix, const unsigned dim) {
    for (unsigned i = 0; i < dim; i++) {
        for (unsigned j = i; j < dim; j++) {
            bMatrix[i][j] = false; 
            bMatrix[j][i] = false;
        }
    }
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
    std::cout << std::endl << "Node " << node.n << ", cost: " << node.cost << std::endl;
    std::cout << "Route " << std::endl;
    tsp::printRoute(node.route);
    std::cout << "Arrows " << std::endl;
    for (unsigned i = 0; i < node.arrows.size(); i++) {
        std::cout << node.arrows[i].first << " " << node.arrows[i].second << ", ";
    }
    std::cout << std::endl << "Prohibited " << std::endl;
    for (unsigned i = 0; i < node.prohibited.size(); i++) {
        std::cout << node.prohibited[i].first << " " << node.prohibited[i].second << ", ";
    }
    std::cout << std::endl << "U " << std::endl;
    for (unsigned i = 0; i < node.u.size(); i++) {
        std::cout << node.u[i] << " ";
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


void doLog(unsigned ub, unsigned lb, unsigned size, unsigned long count, std::string strat) {
    std::cout.precision(3);
    std::cout << "UB: " << ub << std::setw(4) << " LB: " << lb << std::setw(4) 
        << " GAP: " << std::fixed << gap(lb, ub) << "%" << std::setw(4) << " Numero de nos abertos: " << size 
        << std::setw(4) << " Iteracao: " << count << " " << strat << std::endl;
}
