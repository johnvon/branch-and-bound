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
    nodeCurr.cost = hungarian_solve(&p);
    // converte matrix binaria da solucao para rota (possivelmente c/ subciclo)
    nodeCurr.route = getVectorSolution(p.assignment, dim);
    // verifica ciclo e armazena como vetor de pares (arcos a serem proibidos)
    verifyCycle(nodeCurr.route, vCycle, nodeCurr.arrows, dim);
    // finaliza hungaro
    hungarian_free(&p);

    return vCycle;
}

void verifyCycle2() {
}

void verifyCycle(std::vector<int> &sol, std::vector<int> &cycle, std::vector< std::pair<int,int> > &cycleArrows, unsigned dim) {
    bool cDetected = false;
    int cSize = -1, cBegin = 0, cBeginTemp = 0, i;

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
        for (i = 0; i < cycle.size() - 1; i++) {
            cycleArrows.push_back(std::make_pair(cycle[i], cycle[i+1]));
        }
    }
}

bool isRootOptimal(Node& root, const unsigned dim, std::vector<int>& vCycle, unsigned& lb, unsigned& ub) {
    // define lower bound
    lb = root.cost;
    // assignment equivalente a um tour completo = solucao viavel ao TSP
    if (isValidCH(vCycle, dim)) {
        ub = root.cost;
        return true;
    }
    return false;
}

std::vector<int> bestRoute;
int bestCost;

/**
 * Branch-and-bound
 */
void bnb(std::vector<Node>& nodes, const int ** matrix, unsigned dim, unsigned& lb, unsigned& ub) {
    while (!nodes.empty()) {
        double ** dMatrix = copyMatrix2Double(matrix, dim);
        unsigned i, j, nChild = 0, s;

        if (nodes.empty()) {
            std::cerr << "ERRO" << std::endl;
        }
        Node * nodeAux = NULL;
        Node nodeCurr = nodes.back(); // DFS
        nodes.pop_back();

        for (j = 0; j < nodeCurr.arrows.size(); j++) {


            // inviabiliza o uso dos arcos proibidos
            for (i = 0; i < nodeCurr.prohibited.size(); i++) {
                dMatrix[nodeCurr.prohibited[i].first]
                    [nodeCurr.prohibited[i].second] = std::numeric_limits<int>::max();
            }

            std::vector<int> vCycle = hungarian(nodeCurr, dMatrix, dim);

            // branch-and-bound

            // se nao eh uma solucao viavel TSP, ramifica
            if (nodeCurr.route.size() > dim + 1 && nodeCurr.cost < ub) {
                // cada novo no eh uma copia do no atual adicionado de um dos arcos como proibidos
                for (i = 0; i < nodeCurr.arrows.size(); i++) {
                    nodeAux = new Node;
                    nodeAux->n = (int) nodes.size() + 1;
                    nodeAux->cost = nodeCurr.cost;
                    nodeAux->route = nodeCurr.route;
                    //nodeAux->arrows = nodeCurr.arrows;
                    nodeAux->prohibited = nodeCurr.prohibited;
                    nodeAux->prohibited.push_back(nodeCurr.arrows[i]);
                    nodes.push_back(*nodeAux);
                    delete nodeAux;

                    nChild++;
                }
            } else { // caso contrario
                if (nodeCurr.route.size() == dim + 1 && nodeCurr.cost < ub) {
                    if (nodeCurr.cost < bestCost) {
                        ub = nodeCurr.cost;
                        bestRoute = nodeCurr.route;
                        bestCost = nodeCurr.cost;


                        std::cout << "Novo valor Upper Bound: " << ub << std::endl;
                        /*s = 0;
                          std::vector<Node>::iterator prev, it = nodes.begin();
                          for (it != nodes.begin(); it != nodes.end(); ++it) {
                          if (it->cost >= ub) {
                          nodes.erase(it);
                          s++;
                          it--;
                          }
                          }
                          std::cout << s << " nos eliminados por bound, |nodes| = " << nodes.size() << std::endl;
                          */
                        tsp::printVector<int>(nodeCurr.route);
                        std::cout << tsp::cost(nodeCurr.route, matrix) << std::endl;

                    } else {
                        std::cout << "183" << nodes.size() << std::endl;
                    }
                }
                else {
                    std::cout << "187" << nodes.size() << std::endl;
                }
            }

        }
        free<double>(dMatrix, dim);
        //return nodeCurr;
    }
}

void initBranchAndBound(const int ** matrix, const unsigned dim) {
    std::vector<Node> nodes;
    Node nodeCurr;
    unsigned long i;
    unsigned lb = 123456, ub = 1234567;

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
        std::cout << "Raiz: " << std::endl;
        //printNode(nodes.back());
    }
    free<double>(dMatrix, dim);

    // comeca branch-and-bound
    i = 0;
    int x;
    std::cout << "starting bnb " << ub << " " << lb << " " << nodes.size() << " " << (!nodes.empty()) << std::endl;
    bnb(nodes, matrix, dim, lb, ub);

    std::cout << "end bnb" << std::endl;
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
        instance.closeFile();
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
