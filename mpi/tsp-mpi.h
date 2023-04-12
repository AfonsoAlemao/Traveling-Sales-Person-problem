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

Solution *tsp_mpi(Inputs *input, int argc, char *argv[]);
void work(priority_queue_t *queue, int n_cities, double *BestTourCost, Inputs* input, Solution *sol, Path* current_path, int *flag, int rank, int numprocs, double *BestTourCostAuxx, bool comm, int *procBestCost, bool *updateCost, int auxProcBestCost, int receiver, int MYTAG);
#endif


