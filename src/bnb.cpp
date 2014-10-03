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
    std::vector<Node> nodes;
    std::vector<int> bestRoute;
    unsigned lb, ub = inf;

    // raiz
    nodeCurr.n = 0;

    // one tree
    oneTree(nodeCurr, matrix, dim);
    //    printNode(nodeCurr);

    //hungaro
    //double ** cMatrix = copyMatrixFromTo<int, double>(matrix, dim);
    //hungarian(nodeCurr, cMatrix, dim);
    //free<double>(cMatrix, dim);

    if (isRootOptimal(nodeCurr, dim, lb, ub)) {
        std::cout << "Solucao otima encontrada! " << std::endl;
        tsp::printVector<int>(nodeCurr.route);
        return;
    } else {
        nodes.push_back(nodeCurr);
    }

    // comeca branch-and-bound
    bnb(bestRoute, nodes, matrix, dim, lb, ub);
    std::cout << "Rota (custo = " << tsp::cost(bestRoute, matrix) << ")" << std::endl;
    tsp::printVector<int>(bestRoute);
}

void oneTree(Node& node, const int ** matrix, const unsigned dim) {
    int a, fn,sn; // primeiro e segundo clientes
    unsigned * degree = new unsigned[dim];
    bool ** sol1Tree = newBoolMatrix(dim);
    unsigned i, k;

    //    tsp::printMatrix(matrix, dim, 10);
    node.cost = 0;
    for (i = 0; i < dim; i++) {
        degree[i] = 0; // graus
    }

    // MST de vetor[1..dim], i.e., ignorando primeiro elemento
    prim(dim - 1, node.cost, matrix, sol1Tree, degree);

    fn = 1, sn = 2;
    a = inf;
    for (i = 1; i < dim; i++) {
        if (matrix[0][i] < a) {
            sn = fn;
            fn = i;
            a = matrix[0][i];
        }
    }

    // inclui arcos ao custo e a matriz 
    node.cost += matrix[0][fn] + matrix[0][sn];
    sol1Tree[0][fn] = sol1Tree[fn][0] = true;
    sol1Tree[0][sn] = sol1Tree[sn][0] = true;

    // atualiza os graus
    degree[0] = 2;
    degree[fn]++;
    degree[sn]++;
/*     for(i=0;i<dim;i++)
 *         std::cout<<degree[i]<<" ";
 *     std::cout<<std::endl;
 */

    if (isFeasible(dim, degree, k)) {
        std::cout << "Feasible" << std::endl;
        node.route = get1TreeVectorSolution(sol1Tree, dim);
    } else {
        node.arrows.clear();
        for (i = 0; i < dim; i++) {
            if (sol1Tree[k][i]) {
                node.arrows.push_back(std::make_pair(k, i));
            }
        }
    }

    free<bool>(sol1Tree, dim);
    delete[] degree;
}

void prim(const unsigned dim, unsigned& cost, const int ** graph, bool ** sol1Tree, unsigned * degree) {
    bool * selected = new bool[dim];
    int min, x = 0, y = 0; // numero de arestas/edges
    unsigned i, j, ne;

    for (i = 1; i < dim; i++) {
        selected[i] = false;
    }
    selected[0] = true;
    ne = 0;

    while (ne < dim - 1) {
        min = inf;
        for (i = 0; i < dim; i++) {
            if (selected[i]) {
                for (j = 0; j < dim; j++) {
                    if (!selected[j]) {
                        if (graph[i][j] < min) {
                            min = graph[i][j];
                            x = i;
                            y = j;
                        }
                    }
                }
            }
        }

        selected[y] = true;
        // deslocado em 1 para desconsiderar primeiro vertice
        degree[x+1]++;
        degree[y+1]++;

        cost += graph[x][y];
        sol1Tree[y+1][x+1] = sol1Tree[x+1][y+1] = true;
        ne++;
    }
}

bool isFeasible(const unsigned dim, unsigned * degree, unsigned& k) {
    unsigned i, max = 0;
    for (i = 0; i < dim; i++) {
        if (degree[i] > max) {
            max = degree[i];
            k = i;
        }
    }
    return (max==2);
}

std::vector<int> get1TreeVectorSolution(bool ** sol1Tree, const unsigned dim) {
    unsigned i, j, aux;
    std::vector<int> route;
    std::list<int> clients;
    std::list<int>::iterator cit;

    for (i = 1; i < dim; i++) {
        clients.push_back(i);
    }
    route.push_back(0);
    aux = 0;

    while (!clients.empty()) {
        for (cit = clients.begin(); cit != clients.end(); ++cit) {
            j = *cit;
            if (sol1Tree[aux][j] || sol1Tree[j][aux]) {
                if (sol1Tree[aux][j])
                    sol1Tree[aux][j] = false;
                else
                    sol1Tree[j][aux] = false;
                route.push_back(j);
                clients.erase(cit);
                aux = j;
                --cit;
                break;
            }
        }
    }
    route.push_back(0);
    return route;
}

/**
 * Algoritmo hungaro para problema do assignment
 */
void hungarian(Node& nodeCurr, double ** matrix, const int dim) {
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
 * Branch-and-bound
 */
void bnb(std::vector<int>& bestRoute, std::vector<Node>& nodes, const int ** matrix, const unsigned dim, unsigned& lb, unsigned& ub) {
    Node nodeAux;

    //hungaro
    // double ** dMatrix, t;
    //one-tree
    int ** dMatrix, t;

    unsigned i, nBound;
    unsigned long c = 0;
    std::string strat = "";

    // referencia para funcao/estrategia usada
    Node (* curr)(std::vector<Node>&) = NULL;

    int s = 0; //rand() % 3; 
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
        if (c % 2000 == 0) { // log
            std::cout << "UB: " << ub << std::setw(4) << " LB: " << lb << std::setw(4) 
                << " GAP: " <<  gap(lb, ub) << "%" << std::setw(4) << " Numero de nos abertos: " << nodes.size() 
                << std::setw(4) << " Iteracao: " << c << " " << strat << std::endl;
/*             std::cout << "Lista de nos: " << std::endl;
 *             for(std::vector<Node>::iterator it = nodes.begin(); it != nodes.end(); ++it) {
 *                 printNode(*it);
 *             }
 *             std::cout << std::endl;
 * 
 *             std::cin >> t;
 */
        }

        // dMatrix = copyMatrixFromTo<int, double>(matrix, dim);
        dMatrix = copyMatrixFromTo<int, int>(matrix, dim);

        // seleciona (e retira da lista) no de acordo com estrategia
        Node nodeCurr = curr(nodes);

        // inviabiliza o uso dos arcos proibidos
        for (i = 0; i < nodeCurr.prohibited.size(); i++) {
            dMatrix[nodeCurr.prohibited[i].first]
                [nodeCurr.prohibited[i].second] = inf;
            dMatrix[nodeCurr.prohibited[i].second]
                [nodeCurr.prohibited[i].first] = inf;
        }

        // cada novo no eh uma copia do no atual adicionado de um dos arcos como proibidos
        for (i = 0; i < nodeCurr.arrows.size(); i++) {
            //            std::cout << "proibindo " << nodeCurr.arrows[i].first << " " << nodeCurr.arrows[i].second << std::endl;
            c++;
            // branch-and-bound
            nodeAux.n = (int) nodes.size() + 1;
            nodeAux.prohibited = nodeCurr.prohibited;
            nodeAux.prohibited.push_back(nodeCurr.arrows[i]);

            t = dMatrix[nodeCurr.arrows[i].first][nodeCurr.arrows[i].second];
            dMatrix[nodeCurr.arrows[i].first][nodeCurr.arrows[i].second] = inf;
            dMatrix[nodeCurr.arrows[i].second][nodeCurr.arrows[i].first] = inf;

            //hungarian(nodeAux, dMatrix, dim);
            oneTree(nodeAux, const_cast<const int **>(dMatrix), dim);

            // se nao eh uma solucao viavel TSP (ciclo hamiltoniano) mas custo esta abaixo do UB
            if (!nodeAux.arrows.empty() && nodeAux.cost < ub) {//!isValidCH(nodeAux.route, dim) && nodeAux.cost < ub) {
                nodes.push_back(nodeAux);
            } else {
                // solucao viavel e de menor custo, novo UB
                if (nodeAux.arrows.empty() && nodeAux.cost < ub) {// isNewUB(nodeAux, dim, ub)) {
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
                    std::cout << std::setw(4) << " " << nBound << " nos eliminados por bound" << std::endl;
                }
            }
            dMatrix[nodeCurr.arrows[i].first][nodeCurr.arrows[i].second] = t;
            dMatrix[nodeCurr.arrows[i].second][nodeCurr.arrows[i].first] = t;
        }

        // atualiza lower bound
        lb = ub;
        for (i = 0; i < nodes.size(); i++) {
            if (nodes[i].cost < lb)
                lb = nodes[i].cost;
        }
        //        free<double>(dMatrix, dim);
        free<int>(dMatrix, dim);
    }
    std::cout << "FIM Branch-and-bound - LB=" << lb << " UB=" << ub << std::endl;
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

template<typename from, typename to> to ** copyMatrixFromTo(const from ** matrix, const unsigned dim) {
    to ** dMatrix = new to * [dim];
    unsigned i, j;
    for (i = 0; i < dim; i++) {
        dMatrix[i] = new to[dim];
        for (j = 0; j < dim; j++) {
            if (i==j)
                dMatrix[i][j] = (to) inf;
            else
                dMatrix[i][j] = (to) matrix[i][j];
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
