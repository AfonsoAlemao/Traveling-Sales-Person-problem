/**********************************************************************************
* (c) Tomás Fonseca - 66325
*     Afonso Alemão - 96135
*     Rui Daniel    - 96317
*
* Last Modified: 27/02/2023
*
* Name: tsp_omp.c
*
* Description:	Serial implementation of TSP Branch and Bound
*
* Comments:
*
**********************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "tsp-omp.h"
#include "queue.h"

/**********************************************************************************
* tsp_omp()
*
* Arguments:    input - pointer to the input
*
* Returns:      Solutions pointer
*
* Side-Effects: Allocates memory for solution
*
* Description:  Serial implementation of the algorithm that gets
*   the shortest path through all cities
*
**********************************************************************************/

Solution *tsp_omp(Inputs *input) {
    if (input == NULL) return NULL;
    int i = 0, n_cities = get_n_cities(input), *isInTour, *tour;
    double BestTourCost = get_max_value(input), aux_distance = 0, newBound = 0;
    priority_queue_t *queue;
    Path *current_path, *initial_path, *new_path;
    Solution *sol;

    queue = queue_create(compare);
    if (queue == NULL) {
        /* All needed frees and exits in error */
        free_inputs(input);
        error();
    }

    initial_path = create_path(n_cities);
    if (initial_path == NULL) {
        /* All needed frees and exits in error */
        queue_delete(queue);
        free_safe(queue);
        free_inputs(input);
        error();
    }

    /* isInTour will mark the elements in the current tour */
    isInTour = (int *) malloc(n_cities * sizeof(int));
    if (isInTour == NULL) {
        /* All needed frees and exits in error */
        queue_delete(queue);
        free_safe(queue);
        free_inputs(input);
        free_path(initial_path);
        error();
    }

    sol = create_solution(n_cities);
    if (sol == NULL) {
        /* All needed frees and exits in error */
        queue_delete(queue);
        free_safe(queue);
        free_inputs(input);
        free_path(initial_path);
        free_safe(isInTour);
        error();
    }

    set_bound(initial_path, InitialLowerBound(input));

    queue_push(queue, initial_path);

    
    while (queue->size != 0) {
        current_path = queue_pop(queue);

        /* All remaining nodes worse than best */
        if (get_bound(current_path) >= BestTourCost) {
            /* All needed frees */
            free_safe(isInTour);
            while (queue->size != 0) {
                free_path(queue_pop(queue));
            }
            queue_delete(queue);
            free_safe(queue);
            free_path(current_path);

            if (get_BestTourCost(sol) == -1) {
                free_solution(sol);
                return NULL;
            }
            else {
                return sol;
            }
        }

        /* Tour complete, check if it is best */
        if (get_length(current_path) == n_cities) { 
            aux_distance = distance(get_node(current_path), 0, input);
            if (get_cost(current_path) + aux_distance < BestTourCost && aux_distance >= 0) {
                set_BestTour(sol, get_Tour(current_path), n_cities);
                set_BestTour_item(n_cities, sol, 0, n_cities);
                BestTourCost = get_cost(current_path) + aux_distance;
                set_BestTourCost(sol, BestTourCost);
            }
        }
        else {
            tour = get_Tour(current_path);

            /* isInTour mark the elements in the current tour */
            InitializeIsInTour(isInTour, n_cities);

            
            for (i = 0; i < n_cities; i++) {
                if (tour[i] != -1) {
                   isInTour[tour[i]] = 0;
                }
                else {
                    break;
                }
            }
        
            omp_set_num_threads(8);
            #pragma omp parallel private(i, newBound, new_path) 
            {
            #pragma omp for schedule(dynamic, 1)
            for (i = 0; i < n_cities; i++) {
                /* Connection does not exist */
                if (distance(get_node(current_path), i, input) < 0) {
                    continue;
                }

                /* Check if city in already in the current tour, except for the start city (0) */
                if (isInTour[i] == 0 && !(i == 0 && get_length(current_path) == n_cities)) {
                    continue;
                }
                newBound = newLowerBound(input, current_path, i);              

                if (newBound > BestTourCost) { 
                    continue;
                }

                new_path = create_path(n_cities);
                if (new_path == NULL) {
                    /* All needed frees and exits in error */
                    free_safe(isInTour);
                    while (queue->size != 0) {
                        free_path(queue_pop(queue));
                    }
                    queue_delete(queue);
                    free_safe(queue);
                    free_path(current_path);
                    free_solution(sol);
                    free_inputs(input);
                    error();                  
                }
                
                /* Create new path element */
                /* newTour ← Tour ∪ {i} */
                set_Tour(new_path, get_Tour(current_path), n_cities);
                set_Tour_item(get_length(current_path), new_path, i, n_cities);
                /* newCost ← cost + Distances(Node, i) */
                set_cost(new_path, get_cost(current_path) + distance(get_node(current_path), i, input));
                /* newLength ← length + 1 */ 
                set_length(new_path, get_length(current_path) + 1);
                /* newBound ← newBound */ 
                set_bound(new_path, newBound);
                /* newNode ← i*/ 
                set_node(new_path, i);

                #pragma omp critical 
                {
                /* Insert it in queue */
                queue_push(queue, new_path);
                }
            }
            }
        }
        /* Frees already processed element */
        free_path(current_path);
    }

    /* Frees auxiliar structures */
    free_safe(isInTour);
    queue_delete(queue);
    free_safe(queue);

    /* Check if a valid solution was found */
    if (valid_BestTour(sol, n_cities)) {
        set_BestTourCost(sol, BestTourCost);
        return sol;
    }
    free_solution(sol);
    return NULL;
}
