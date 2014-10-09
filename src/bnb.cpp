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
    std::list<Node> nodes;
    std::vector<int> bestRoute;
    unsigned lb, ub;

    // heuristica construtiva para definir UB
    Construct c(matrix, dim);
    c.nearestNeighbor();
    ub = tsp::cost(c.getRoute(), matrix);

    unsigned x = 1;
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
    int ** cMatrix = copyMatrixFromTo<int, int>(matrix, dim);
    bool ** sol1Tree  = newBoolMatrix(dim); // Total de arestas em cada vertice
    oneTree(nodeCurr, cMatrix, dim, sol1Tree);
    free<int>(cMatrix, dim);
    free<bool>(sol1Tree, dim);

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

/**
 * Branch-and-bound
 */
void bnb(std::vector<int>& bestRoute, std::list<Node>& nodes, const int ** matrix, const unsigned dim, unsigned& lb, unsigned& ub, unsigned b, unsigned x) {
    Node nodeAux;
    bool   ** sol1Tree = newBoolMatrix(dim);
    double ** dMatrix  = NULL, d = 0;
    int    ** cMatrix  = NULL, c = 0;
    time_t t0, t1, t2;
    unsigned i, j, nBound;
    unsigned long count = 0;

    std::string strat = "";
    std::list<Node>::iterator prev, it, itCurr;
    // referencia para funcao/estrategia usada
    std::list<Node>::iterator (* curr)(std::list<Node>&) = NULL;

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
        case 3:
            curr = &randb;
            strat = "RAND-BOUND";
            break;
    }

    t0 = t1 = t2 = time(NULL); 
    doLog(ub, lb, nodes.size(), count, strat);

    while (!nodes.empty()) {
        if (difftime(time(NULL),t1) > 1) {
            time(&t1);
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
        itCurr = curr(nodes);

        // atualiza lower bound
        if (itCurr->cost > lb) {
            updateLB(nodes, lb);
        }

        // fortalece os bounds
        if (x==0) {
            for (i = 0; i < dim; i++) {
                for (j = 0; j < dim; j++) {
                    if (cMatrix[i][j] != inf)
                        cMatrix[i][j] += itCurr->pi[i] + itCurr->pi[j];
                }
            }
        }

        // inviabiliza o uso dos arcos proibidos
        for (i = 0; i < itCurr->prohibited.size(); i++) {
            if (x) { // hungaro
                dMatrix[itCurr->prohibited[i].first][itCurr->prohibited[i].second] = inf;
            } else { // 1-tree
                cMatrix[itCurr->prohibited[i].first][itCurr->prohibited[i].second] = inf;
                cMatrix[itCurr->prohibited[i].second][itCurr->prohibited[i].first] = inf;
            }
        }

        // cada novo no eh uma copia do no atual adicionado de um dos arcos como proibidos
        for (i = 0; i < itCurr->arrows.size(); i++) {
            // branch-and-bound
            nodeAux.n = (int) nodes.size() + 1;
            nodeAux.prohibited = itCurr->prohibited;
            nodeAux.prohibited.push_back(itCurr->arrows[i]);
            nodeAux.route.clear();

            if (x) { // hungaro
                d = dMatrix[itCurr->arrows[i].first][itCurr->arrows[i].second];
                dMatrix[itCurr->arrows[i].first][itCurr->arrows[i].second] = inf;
                hungarian(nodeAux, dMatrix, dim);
            } else { // 1-tree
                c = cMatrix[itCurr->arrows[i].first][itCurr->arrows[i].second];
                cMatrix[itCurr->arrows[i].first][itCurr->arrows[i].second] = inf;
                cMatrix[itCurr->arrows[i].second][itCurr->arrows[i].first] = inf;
                oneTree(nodeAux, cMatrix, dim, sol1Tree);
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
                    for (it = nodes.begin(); it != nodes.end();) {
                        if (it->cost >= ub) {
                            it = nodes.erase(it);
                            nBound++;
                        } else {
                            ++it;
                        }
                    }
                }
            }

            if (x) { // hungaro 
                dMatrix[itCurr->arrows[i].first][itCurr->arrows[i].second] = d;
            } else { // 1-tree
                cMatrix[itCurr->arrows[i].first][itCurr->arrows[i].second] = c;
                cMatrix[itCurr->arrows[i].second][itCurr->arrows[i].first] = c;
            }
        }

        if (x)
            free<double>(dMatrix, dim);
        else
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
void oneTree(Node& node, int ** matrix, const unsigned dim, bool ** sol1Tree) {
    unsigned * degree = new unsigned[dim]; // Total de arestas em cada vertice
    unsigned fn, sn, k, i, cost = 0;

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

    delete[] degree;
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


void doLog(unsigned ub, unsigned lb, unsigned size, unsigned long count, std::string strat) {
    std::cout << "UB: " << ub << std::setw(4) << " LB: " << lb << std::setw(4) 
        << " GAP: " <<  gap(lb, ub) << "%" << std::setw(4) << " Numero de nos abertos: " << size 
        << std::setw(4) << " Iteracao: " << count << " " << strat << std::endl;
}
