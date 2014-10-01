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

    double bef = tsp::cpuTime();
    initBranchAndBound(instance.getMatrix(), instance.getDim());
    double aft = tsp::cpuTime();

    std::cout << (aft-bef) << "ms" << std::endl;

    return 0;
}


