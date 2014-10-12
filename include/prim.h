/*
 * @file prim.h
 *      
 * @author: Fabio Cruz B. de Albuquerque
 * 
 * @date: 04/10/14
 */

#ifndef PRIM_H_
#define PRIM_H_

#include <iostream>
#include <iomanip>

const int inf = 99999;

template<typename type>
int minKey(type key[], bool mstSet[], const unsigned dim);

template<typename type>
void prim1Tree(const unsigned dim, type ** graph, unsigned * degree, bool ** sol1Tree, type &cost);

extern template int minKey<int>(int key[], bool mstSet[], const unsigned dim);
extern template int minKey<double>(double key[], bool mstSet[], const unsigned dim);

extern template void prim1Tree<int>(const unsigned dim, int ** graph, unsigned * degree, bool ** sol1Tree, int &cost);
extern template void prim1Tree<double>(const unsigned dim, double ** graph, unsigned * degree, bool ** sol1Tree, double &cost);


#endif
