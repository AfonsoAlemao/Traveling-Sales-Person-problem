/**********************************************************************************
* (h) Tomás Fonseca - 66325
*     Afonso Alemão - 96135
*     Rui Daniel    - 96317
*
* Last Modified: 23/03/2023
*
* Name: tsp_omp.h
*
* Description: tsp_omp.c header
*
* Comments:
*
**********************************************************************************/

#ifndef _TSP_OMP_H
#define _TSP_OMP_H

#include <omp.h>
#include <mpi.h>

#include "auxiliar.h"
#include "queue.h"

Solution *tsp_omp(Inputs *);
void work(priority_queue_t *, int , double *, Inputs* , Solution *, Path*, int *);

#endif


