/**********************************************************************************
* (h) Tomás Fonseca - 66325
*     Afonso Alemão - 96135
*     Rui Daniel    - 96317
*
* Last Modified: 12/04/2023
*
* Name: tsp_mpi.h
*
* Description: tsp_mpi.c header
*
* Comments:
*
**********************************************************************************/

#ifndef _TSP_OMP_H
#define _TSP_OMP_H

#include <omp.h>

#include "auxiliar.h"
#include "queue.h"

Solution *tsp_mpi(Inputs *);
void work(priority_queue_t *, int , double *, Inputs* , Solution *, Path*, int *);

#endif


