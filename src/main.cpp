/*
 * @file main.cpp
 *      
 * @author: Fabio Cruz B. de Albuquerque
 * 
 * @date: 27/09/14
 */

#include "../include/headers.h"

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

std::vector<int> hungarian(Node& nodeCurr, double ** matrix, const int dim) {
    int mode;
    std::vector<int> vCycle;

    // chamada ao hungaro para resolucao do problema de assignment
    hungarian_problem_t p;
    mode = HUNGARIAN_MODE_MINIMIZE_COST;
    hungarian_init(&p, matrix, dim, dim, mode);
    nodeCurr.solution = hungarian_solve(&p);
    // converte matrix binaria da solucao para rota (possivelmente c/ subciclo)
    nodeCurr.route = getVectorSolution(p.assignment, dim);
    // verifica ciclo e armazena como vetor de pares (arcos a serem proibidos)
    verifyCycle(nodeCurr.route, vCycle, nodeCurr.arrows, dim);
    // finaliza hungaro
    hungarian_free(&p);

    return vCycle;
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

bool isRootOptimal(Node& root, const unsigned dim, std::vector<int>& vCycle, unsigned& lb, unsigned& ub) {
    // define lower bound
    lb = root.solution;
    // assignment equivalente a um tour completo = solucao viavel ao TSP
    if (isValidCH(vCycle, dim)) {
        ub = root.solution;
        return true;
    }
    return false;
}

/**
 * Branch-and-bound
 */
Node bnb(std::vector<Node>& nodes, const int ** matrix, unsigned dim, unsigned& lb, unsigned& ub) {
    double ** dMatrix = copyMatrix2Double(matrix, dim);
    unsigned i, nChild = 0;

    Node nodeAux;
    Node nodeCurr = nodes.back(); // DFS
    nodes.pop_back();
    std::cout << nodeCurr.n << std::endl;

    // inviabiliza o uso dos arcos proibidos
    for (i = 0; i < nodeCurr.prohibited.size(); i++) {
        dMatrix[nodeCurr.prohibited[i].first]
            [nodeCurr.prohibited[i].second] = std::numeric_limits<int>::max();
    }

    std::vector<int> vCycle = hungarian(nodeCurr, dMatrix, dim);

    // branch-and-bound
    if (nodeCurr.route.size() > dim + 1 && nodeCurr.solution < ub) {
        // cada novo no eh uma copia do no atual adicionado de um dos arcos como proibidos
        for (i = 0; i < nodeCurr.arrows.size(); i++) {
            nodeAux.n = (int) nodes.back().n + 1 + i;
            nodeAux.solution = nodeCurr.solution;
            nodeAux.route = nodeCurr.route;
            nodeAux.arrows = nodeCurr.arrows;
            nodeAux.prohibited = nodeCurr.prohibited;
            nodeAux.prohibited.push_back(nodeCurr.arrows[i]);
            nodes.push_back(nodeAux);
            nChild++;
        }
        //          std::cout << "branching de " << nChild << std::endl;
    } else { 
        if (nodeCurr.route.size() == dim + 1 && nodeCurr.solution >= lb && nodeCurr.solution < ub) {
            std::cout << "Novo valor Upper Bound" << std::endl;
            //                best.swap(nodeCurr.route);
            ub = nodeCurr.solution;
            int s = 0;
            std::vector<Node>::iterator prev, it = nodes.begin();
            while (it != nodes.end()) {
                if (it->solution >= ub) {
                    nodes.erase(it);
                    it--;
                    s++;
                }
                it++;
            }
            std::cout << s << " nos eliminados por bound, |nodes| = " << nodes.size() << std::endl;
            tsp::printVector<int>(nodeCurr.route);
        }
    }

    free<double>(dMatrix, dim);
    return nodeCurr;
}

void initBranchAndBound(const int ** matrix, const unsigned dim) {
    std::vector<Node> nodes;
    Node nodeCurr;
    unsigned i;
    unsigned lb, ub = 1234567;

    double ** dMatrix = copyMatrix2Double(matrix, dim);

    // raiz
    nodeCurr.n = 0;
    std::vector<int> vCycle = hungarian(nodeCurr, dMatrix, dim);
    if (isRootOptimal(nodeCurr, dim, vCycle, lb, ub)) {
        std::cout << "Solucao otima encontrada! " << std::endl;
        tsp::printVector<int>(nodeCurr.route);
        return;
    } else {
        nodes.push_back(nodeCurr);
    }
    free<double>(dMatrix, dim);

    // comeca branch-and-bound
    i = 0;
    while (!nodes.empty()) {
        i++;
        nodeCurr = bnb(nodes, matrix, dim, lb, ub);
        std::cout << "(" << i << ") LB= " << lb << ", UB=" << ub << std::endl;
        std::cout << "|nodes| = " << nodes.size() << std::endl;
    }
    std::cout << tsp::cost(nodeCurr.route, matrix) << std::endl;
}

int main(int argc, char *argv[]) {
    std::string file = "";
    std::stringstream ss;

    if (argc >= 2) {
        file.assign(argv[1]);
    } else {
        std::cout << "Uso " << argv[0] << " arquivo_entrada" << std::endl;
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

    initBranchAndBound(instance.getMatrix(), instance.getDim());
    std::cout << "OK" << std::endl;

    return 0;
}

double ** copyMatrix2Double(const int ** matrix, const unsigned dim) {
    double ** dMatrix = new double * [dim];
    unsigned i, j;
    for (i = 0; i < dim; i++) {
        dMatrix[i] = new double[dim];
        for (j = 0; j < dim; j++) {
            if (i==j)
                dMatrix[i][j] = std::numeric_limits<int>::max();
            else
                dMatrix[i][j] = matrix[i][j];
        }
    }
    return dMatrix;
}

template<typename type> 
void free(type ** matrix, const unsigned dim) {
    unsigned i;
    for (i = 0; i < dim; i++) {
        delete[] matrix[i];
    }
    delete[] matrix;
}

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
