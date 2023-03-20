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
    int n_cities = get_n_cities(input), whistle = -1, global_tid = 0;
    double BestTourCost = get_max_value(input);
    Path *initial_path;
    Solution *sol;

    initial_path = create_path(n_cities);
    if (initial_path == NULL) {
        /* All needed frees and exits in error */
        free_inputs(input);
        error();
    }

    sol = create_solution(n_cities);
    if (sol == NULL) {
        /* All needed frees and exits in error */
        free_inputs(input);
        free_path(initial_path);
        error();
    }

    set_bound(initial_path, InitialLowerBound(input));

    
    Path *new_path[32];
    int exit_global = 0;
    priority_queue_t *queue[32];

    omp_set_num_threads(8);
    #pragma omp parallel shared(BestTourCost)
    {
        // printf("%d\n",omp_get_thread_num());

        int tid = omp_get_thread_num();

        queue[tid] = queue_create(compare);
        Path *current_path;
        int flag = 0, pop_counter = 0;
        
        if (queue[tid] == NULL) {
            /* All needed frees and exits in error */
            error();
        }

        #pragma omp single
        {
            queue_push(queue[tid], initial_path);
        }


        while ((queue[tid]->size != 0) && (flag != 1) && (queue[tid]->size < omp_get_num_threads())) {

            // printf("%d\n",omp_get_thread_num());

            // printf("Thread %d\n", omp_get_thread_num());
            // printf("queuesize %ld\n", queue->size);
            current_path = queue_pop(queue[tid]);
            work(queue[tid], n_cities, &BestTourCost, input, sol, current_path, &flag);
            free_path(current_path);    
        }

        if (flag == 1) {
            #pragma omp atomic
                exit_global += 1;
        }
        //printf("1: Thread: %d, flag = %d, queue[tid]-size = %ld\n",omp_get_thread_num(), flag, queue[tid]->size);
        //printf("exit = %d\n", exit_global);

        #pragma omp barrier

        if ((exit_global != 1)) {
                // printf("%d\n",omp_get_thread_num());

            if (queue[tid]->size != 0) {
                //printf("Olá, sou a thread %d\n", omp_get_thread_num());
                for (int i = 0 ; i < omp_get_num_threads(); i++) {
                    new_path[i] = queue_pop(queue[tid]);
                }
            }
            #pragma omp barrier
            
            //printf("2: Thread: %d, flag = %d, queue[tid]-size = %ld\n",omp_get_thread_num(), flag, queue[tid]->size);

            #pragma omp for
            for (int i = 0 ; i < omp_get_num_threads(); i++) {
                queue_push(queue[tid], new_path[i]);
                // printf("%d\n",omp_get_thread_num());

            }

            //printf("3: Thread: %d, flag = %d, queue[tid]-size = %ld\n",omp_get_thread_num(), flag, queue[tid]->size);

        }

        while ((queue[tid]->size != 0) && (flag != 1)) {

            // printf("%d\n",omp_get_thread_num());

            // printf("%d\n", omp_get_thread_num());
            current_path = queue_pop(queue[tid]);
            pop_counter++;
            work(queue[tid], n_cities, &BestTourCost, input, sol, current_path, &flag);
            free_path(current_path);
            if (whistle != -1) { //flag tem o numero da thread que precisa de uma thread
                if (whistle == global_tid) {
                    if (global_tid < omp_get_num_threads() - 1) {
                        global_tid += 1;
                    }
                    else {
                        global_tid = 0;
                    }
                }
                
                if (tid == global_tid) {
                    #pragma omp critical 
                    {
                        //printf("Thread %d vai mandar um push para o whistle %d\n", tid, whistle);
                        //printf("Thread %d queue size %d\n", tid, queue[tid]->size);
                        //printf("Thread %d queue size %d\n", whistle, queue[whistle]->size);
                        if (queue[tid]->size != 0 && whistle != -1) {
                                
                            // int distribute = 1 + queue[tid]->size / (omp_get_num_threads() * 8);
                            // int distribute = 2;
                            // for (int ii = 0; ii < distribute && queue[tid]->size != 0; ii++) {
                            queue_push(queue[whistle], queue_pop(queue[tid]));
                            // }
                            // global_tid = 0;
                            whistle = -1;  
                        }
                        else {
                            if (global_tid < omp_get_num_threads() - 1) {
                                global_tid += 1;
                            }
                            else {
                                global_tid = 0;
                            }
                            
                        } 
                    }
                }
            }
            if (!((queue[tid]->size != 0) && (flag != 1))) {
                
                // printf("%d\n",omp_get_thread_num());
                
                while (queue[tid]->size != 0) {
                    free_path(queue_pop(queue[tid]));
                }

                //thread vazia
                #pragma omp critical 
                {
                    if (whistle == -1) {
                        whistle = tid;
                    }
                }

                int count = 0;
                while(whistle == tid) {
                    
                    // printf("1: %d\n",omp_get_thread_num());
                    // printf("whistle = %d\n", whistle);
                    
                    count = 0;
                    for (int q = 0; q < omp_get_num_threads(); q++){
                        // printf("size = %d\n", queue[q]->size);
                        if (queue[q]->size != 0) {
                            count = 1;
                            break;
                        }
                    }
                    if (count == 0) {
                        break;
                    }
                }

                // if (whistle != -1) {
                // printf("Thread %d \t whistle = %d \t count = %d, queue->size = %ld\n", tid, whistle, count, queue[tid]->size);
                // }

                if (queue[tid]->size != 0 || count != 0) {
                   flag = 0; 
                }
                
                // #pragma omp critical 
                // {
                //     whistle = -1;
                // }
                
            }    
        }
        

        

        while (queue[tid]->size != 0) {
            free_path(queue_pop(queue[tid]));
        }

        printf("Thread %d, pop_counter = %d\n", omp_get_thread_num(), pop_counter);
        
        #pragma omp barrier
        /* Frees auxiliar structures */
        queue_delete(queue[tid]);
        free_safe(queue[tid]);
    }

    

    /* Check if a valid solution was found */
    if (valid_BestTour(sol, n_cities)) {
        return sol;
    }
    free_solution(sol);
    return NULL;
}

void work(priority_queue_t *queue, int n_cities, double *BestTourCost, Inputs* input, Solution *sol, Path* current_path, int *flag) {
    int i = 0;
    int *isInTour, *tour;
    double newBound = 0, aux_distance = 0;
    Path *new_path;
    
    /* All remaining nodes worse than best */
    if (get_bound(current_path) >= *BestTourCost) {
        *flag = 1;
    }

    /* isInTour will mark the elements in the current tour */
    isInTour = (int *) malloc(n_cities * sizeof(int));
    if (isInTour == NULL) {
        *flag = 2;
        error();
    }
    //#pragma omp barrier

    /* Tour complete, check if it is best */
    if (get_length(current_path) == n_cities && *flag == 0) { 
        aux_distance = distance(get_node(current_path), 0, input);
        #pragma omp critical
        {
            if (get_cost(current_path) + aux_distance < *BestTourCost && aux_distance >= 0) {
                set_BestTour(sol, get_Tour(current_path), n_cities);
                set_BestTour_item(n_cities, sol, 0, n_cities);
                *BestTourCost = get_cost(current_path) + aux_distance;
                set_BestTourCost(sol, *BestTourCost);
            }

        }
        

    }
    
    else if (*flag == 0) {
        tour = get_Tour(current_path);

        /* isInTour mark the elements in the current tour */
        InitializeIsInTour(isInTour, n_cities);

        for (i = 0; i < n_cities; i++) {
            if (tour[i] != -1) {
                isInTour[tour[i]] = 0;
            }
        }
    
        
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

            if (newBound > *BestTourCost) { 
                continue;
            }

            new_path = create_path(n_cities);
            if (new_path == NULL) {
                *flag = 2;             
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

            /* Insert it in queue */
            queue_push(queue, new_path);
            
        }
        
    }
    // free_path(current_path);
    free_safe(isInTour);

    return;
}