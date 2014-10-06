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

int minKey(int key[], bool mstSet[], const unsigned dim);
void printMST(int parent[], int ** graph, const unsigned dim);
void prim1Tree(const unsigned dim, int ** graph, unsigned * degree, bool ** sol1Tree, unsigned &cost);
#endif
