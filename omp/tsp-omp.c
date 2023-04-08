/**********************************************************************************
* (c) Tomás Fonseca - 66325
*     Afonso Alemão - 96135
*     Rui Daniel    - 96317
*
* Last Modified: 23/03/2023
*
* Name: tsp_omp.c
*
* Description:	OpenMP implementation of TSP Branch and Bound
*
* Comments:
*
**********************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "tsp-omp.h"

#define MAX_N_THREADS 64

/**********************************************************************************
* tsp_omp()
*
* Arguments:    input - pointer to the input
*
* Returns:      Solutions pointer
*
* Side-Effects: Allocates memory for solution
*
* Description:  OpenMP implementation of the algorithm that gets
*   the shortest path through all cities
*
**********************************************************************************/

Solution *tsp_omp(Inputs *input) {
    if (input == NULL) return NULL;
    int n_cities = get_n_cities(input), whistle = -1, global_tid = 0, exit_global = 0, twice = 0, dealer = 0;
    double BestTourCost = get_max_value(input);
    Path *initial_path, *new_path[MAX_N_THREADS * 2];
    Solution *sol;
    priority_queue_t *queue[MAX_N_THREADS];
        
    /* density = edges / (2 * nodes) */
    int twice_density = get_n_edges(input) / get_n_cities(input);

    initial_path = create_path(n_cities);
    if (initial_path == NULL) {
        /* All needed frees and exits in error. */
        free_inputs(input);
        error();
    }

    sol = create_solution(n_cities);
    if (sol == NULL) {
        /* All needed frees and exits in error. */
        free_inputs(input);
        free_path(initial_path);
        error();
    }

    set_bound(initial_path, InitialLowerBound(input));

    /* Creates N parallel threads. All threads execute the subsequent block.
    All threads wait for each other at the end of this executing block: implicit barrier synchronization */
    /* BestTourCost is a shared variable that exists in a single location 
    and all threads can read and write it. */
    #pragma omp parallel shared(BestTourCost, dealer)
    {
        /* Each thread gets its thread id. */
        int tid = omp_get_thread_num();
        
        Path *current_path;
        int flag = 0, count = 0;
        
        /* Creates a queue for each thread. */
        queue[tid] = queue_create(compare);
        if (queue[tid] == NULL) {
            /* Exits in error */
            error();
        }
 
        /* Single thread execute this region: the master thread. The 1st element goes to a single thread's queue. */
        #pragma omp master
        {
            queue_push(queue[tid], initial_path);
        
            /* Relevant computation to secure that further ahead we
            can evenly distribute elements by the queue of each thread. */
            if (((int) (twice_density * 0.8)) < 1) {
                dealer = omp_get_num_threads();
            }
            else {
                dealer = (int) (omp_get_num_threads() * twice_density * 0.8);
            }
        }

        /* Process (Multiple: Pop + Work) the queue that has received first element until:
        - the program ends: flag = 1 (only irrelevant elements in queue) or queue[tid]->size = 0 (no more elements)
        or
        - the queue size is large enough to evenly distribute elements (1 or 2) among the queues of the 
        other threads. This queue size and the number of distributed elements were optimized through a 
        calculation that uses the density of the graph in question.  */
        while (((int) queue[tid]->size != 0) && (flag != 1) && ((int) queue[tid]->size < dealer)) {
            current_path = queue_pop(queue[tid]);
            work(queue[tid], n_cities, &BestTourCost, input, sol, current_path, &flag);
            free_path(current_path);  
        }

        #pragma omp master
        {
            /* Checks if all the processing of the graph was previously done. */
            if (flag == 1 || (int) queue[tid]->size == 0) {
                /* Guarantees that reading and writing of one memory location is atomic. */
                #pragma omp atomic
                    exit_global += 1;
            }
        }

        /* All threads wait for each other at the end of the executing block. */
        #pragma omp barrier

        /* Evenly distribute elements (1 or 2) among the queues of the 
        other threads. */
        if (exit_global != 1) {
            /* Check which thread have the elements to distribute. */
            if ((int) queue[tid]->size != 0) {
                /* Check if the number of elements to distribute allows to distribute 
                1 or 2 elements per queue thread. Then pop that elements. */
                if ((int) queue[tid]->size > omp_get_num_threads() * 2 + 2) {
                    #pragma omp atomic
                        twice += 1;
                    for (int i = 0 ; i < omp_get_num_threads() * 2; i++) {  
                        new_path[i] = queue_pop(queue[tid]);
                    }
                }
                else {
                    for (int i = 0 ; i < omp_get_num_threads(); i++) {
                        new_path[i] = queue_pop(queue[tid]);
                    }
                }
            }
            #pragma omp barrier


            /* Distribute the elements: twice > 0 means that 2 elements will be distributed for each thread queue, 
            otherwise distribute 1.
            The distribution of the elements is made taking into consideration a balance relative to 
            the priority of each element: threads will receive the elements with highest priorities. */
            queue_push(queue[tid], new_path[tid]);
            if (twice) {
                queue_push(queue[tid], new_path[omp_get_num_threads() + tid]);
            }
        }


        /* Process the queue of each thread until the program ends: 
        flag = 1 (only irrelevant elements in queues) or queue size = 0 for all threads. */
        while (((int) queue[tid]->size != 0) && (flag != 1)) {
            current_path = queue_pop(queue[tid]);
            work(queue[tid], n_cities, &BestTourCost, input, sol, current_path, &flag);
            free_path(current_path);

            /* Load balancing: when a thread finish the execution of its tasks it asks for elements 
            from the other threads in order to process them. The process is the following:
            We have a shared variable whistle. When a thread ask elements it will change the value of whistle 
            to its thread id. Then the other threads will check if there is a whistling thread (whistle != -1). 
            In that case one thread, whose id is equal to global_tid (shared variable), will give an element 
            (with high priority) to the whistling thread.
            The global tid is initialized to 0 because it is the master's id. Then the global_tid will change when necessary:
            when the thread whose id is the global tid can't give more elements (because its size is 0) or when the
            whistling thread id is equal to global id. 
            We use critical regions in order to synchronize the process of communication between threads. */

            /* Section dedicated to give elements to the whistling thread. */
            if (whistle != -1) {
                if (whistle == global_tid) {
                    if (global_tid < omp_get_num_threads() - 1) {
                        global_tid += 1;
                    }
                    else {
                        global_tid = 0;
                    }
                }
                
                if (tid == global_tid) {

                    /* A thread waits at the the beginning of a critical region until no
                    other thread is executing a critical region with the same name. */
                    #pragma omp critical 
                    {
                        if ((int) queue[tid]->size != 0 && whistle != -1) {
                            queue_push(queue[whistle], queue_pop(queue[tid]));
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

            /* Section dedicated to the asking / receiving of elements by the whistling thread. */
            if (!(((int) queue[tid]->size != 0) && (flag != 1))) {
                
                /* If the whistle thread have only irrelevant elements in queue, free them. */
                while ((int) queue[tid]->size != 0) {
                    free_path(queue_pop(queue[tid]));
                }

                /* Become the whistling thread. */
                #pragma omp critical 
                {
                    if (whistle == -1) {
                        whistle = tid;
                    }
                }

                /* Ask for elements while the other threads still have elements in their queues. */
                while (whistle == tid && queue[tid]->size != 0) {
                    count = 0;

                    for (int k = 0; k < omp_get_num_threads(); k++){
                        if ((int) queue[k]->size != 0) {
                            count = 1;
                            break;
                        }
                    }
                    if (count != 1) {
                        break;
                    }
                }

                /* Check if the whistling thread received any element. */
                if ((int) queue[tid]->size != 0 || count != 0) {
                   flag = 0; 
                }             
            }    
        }

        /* Free the irrelevant elements in queue. */
        while ((int) queue[tid]->size != 0) {
            free_path(queue_pop(queue[tid]));
        }

        #pragma omp barrier

        /* Frees auxiliar structures. */
        queue_delete(queue[tid]);
        free_safe(queue[tid]);
    }

    /* Check if a valid solution was found. */
    if (valid_BestTour(sol, n_cities)) {
        return sol;
    }
    free_solution(sol);
    return NULL;
}


/**********************************************************************************
* work()
*
* Arguments:    queue - pointer to the queue
*               n_cities - number of nodes in the graph
*               BestTourCost - pointer to the current best solution cost for the problem
*               input - pointer to the input
*               sol - pointer to the solution
*               current_path - path element to be process
*               flag - check if the queue has any useful work
*
* Returns:      (void)
*
* Side-Effects: (none)
*
* Description:  Process a popped path element from the queue and executes 
*    the problem algorithm, pushing new elements to the queue
*
**********************************************************************************/

void work(priority_queue_t *queue, int n_cities, double *BestTourCost, Inputs* input, Solution *sol, Path* current_path, int *flag) {
    int i = 0, *isInTour, *tour;
    double newBound = 0, aux_distance = 0;
    Path *new_path;
    
    /* Checks if all remaining nodes in queue are worse than BestTourCost. */
    if (get_bound(current_path) >= *BestTourCost) {
        *flag = 1;
    }

    /* isInTour will mark the elements in the current tour. */
    isInTour = (int *) malloc(n_cities * sizeof(int));
    if (isInTour == NULL) {
        *flag = 2;
        error();
    }

    /* Tour complete, check if it is best. */
    if (get_length(current_path) == n_cities && *flag == 0) { 
        aux_distance = distance(get_node(current_path), 0, input);
        /* Update the best solution tour and cost. We are updating the solution that is a
        shared structure, so we need a critical region to secure that a thread waits at the 
        the beginning of that block until no other thread is executing that section. */
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

        /* isInTour mark the elements in the current tour. */
        InitializeIsInTour(isInTour, n_cities);

        for (i = 0; i < n_cities; i++) {
            if (tour[i] != -1) {
                isInTour[tour[i]] = 0;
            }
        }

        for (i = 0; i < n_cities; i++) {
            /* Connection does not exist. */
            if (distance(get_node(current_path), i, input) < 0) {
                continue;
            }

            /* Check if city in already in the current tour, except for the start city (0). */
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
            
            /* Create new path element. */
            /* newTour ← Tour ∪ {i} */
            set_Tour(new_path, get_Tour(current_path), n_cities);
            set_Tour_item(get_length(current_path), new_path, i, n_cities);
            /* newCost ← cost + Distances(Node, i) */
            set_cost(new_path, get_cost(current_path) + distance(get_node(current_path), i, input));
            /* newLength ← length + 1 */ 
            set_length(new_path, get_length(current_path) + 1);
            /* newBound ← newBound */ 
            set_bound(new_path, newBound);
            /* newNode ← i */ 
            set_node(new_path, i);

            /* Insert it in queue. */
            queue_push(queue, new_path);
            
        }  
    }

    free_safe(isInTour);

    return;
}