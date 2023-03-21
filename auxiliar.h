/**********************************************************************************
* (h) Tomás Fonseca - 66325
*     Afonso Alemão - 96135
*     Rui Daniel    - 96317
*
* Last Modified: 23/03/2023
*
* Name: auxiliar.h
*
* Description: auxiliar.c header
*
* Comments:
*
**********************************************************************************/

#ifndef _AUXILIAR_H
#define _AUXILIAR_H

#include "structs.h"
#include "queue.h"

/* Methods */

double **aloc2d(int, int);
char compare (void *, void *);
double distance(int, int, Inputs *);
void error();
void free2d(double **, int);
void *free_safe (void *);
void InitializeIsInTour(int *, int);
double InitialLowerBound(Inputs *);
double newLowerBound(Inputs *, Path *, int);
FILE *OpenFile(char *, char *);
Inputs *parse_inputs(int, char **);
void print_matrix(int, Inputs *);
void print_result(Solution *, int);
Inputs *ReadFileIn(FILE *, double);
void ValidInputFileName(char *);

#endif
