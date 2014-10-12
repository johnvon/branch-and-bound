/*
 * @file main.cpp
 *      
 * @author: Fabio Cruz B. de Albuquerque
 * 
 * @date: 27/09/14
 */

#include "../include/headers.h"

#include <exception>
#include <iostream>
#include <sstream>
#include <string>

#include "Instance.h"
#include "Util.h"

int main(int argc, char *argv[]) {
    std::string file = "";
    std::stringstream ss;
    unsigned s = 0, x = 1;
    srand(time(NULL));

    if (argc >= 2) {
        file.assign(argv[1]);
        if (argc >= 3) {
            ss.str(argv[2]);
            ss >> s;
            if (argc >= 4) {
                ss.clear();
                ss.str(argv[3]);
                ss >> x;
            }
        }
    } else {
        std::cout << "Uso " << argv[0] << " arquivo_entrada [relaxacao(0=hungaro,1=1tree,2=lagrangeana)] [branching(0=dfs,1=bfs,2=bestb)]" << std::endl;
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

    double bef = tsp::cpuTime();
    initBranchAndBound(instance.getMatrix(), instance.getDim(), x, s);
    double aft = tsp::cpuTime();
    std::cout << (aft-bef) << "ms" << std::endl;

    return 0;
}


