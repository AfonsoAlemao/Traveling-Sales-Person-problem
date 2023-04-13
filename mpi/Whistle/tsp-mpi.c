/**********************************************************************************
* (c) Tomás Fonseca - 66325
*     Afonso Alemão - 96135
*     Rui Daniel    - 96317
*
* Last Modified: 12/04/2023
*
* Name: tsp_mpi.c
*
* Description:	MPI + OpenMP implementation of TSP Branch and Bound
*
* Comments:
*
**********************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "tsp-mpi.h"


#define MAX_N_THREADS 64
#define MAX_N_PROCS 8
#define MYTAG 1010
#define MYTAG2 2211
#define TAG_SOS 2233
#define TAG_WORK_SOS 3434
#define TAG_SOSFINAL 4545
#define INFINTY_DOUBLE 1.7E+308


typedef struct _path {
    double cost;
    double bound;
	int length;
	int node;
    long isInTour;
    int Tour[];
}Path;

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

Solution *tsp_mpi(Inputs *input, int argc, char *argv[]) {
    int numprocs = 1, rank = 0;
    int sender = 0, receiver = 0;
    int solved = 0;
    int workStatus = 1;
    int workStatusAux = 1;
    bool mega_global_exit = false, listenSOS = false;
    int whistlerSOS = -1;
    int flag_recSOS = 0;
    bool final = false, listenSOS2 = false, giveWork = false;
    int flag_recSOSmain = 0, sos_received_aux = 0, flag_recSOS2 = 0, sos_received_count = 0;
    int count8 = 0, count10=0;


    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Request req, reqsend, req2, reqsend2, reqSOSmain, reqSOS, reqSOSsender, reqSOS2, reqsend8, reqsend9, reqsend10, reqsend11, reqsend12;
    MPI_Status stat, stat2, statusSOSmain, statusSOS, statusSOS2;

    


    if (input == NULL) return NULL;
    double BestTourCost = get_max_value(input), BestTourCostAuxx[2];
    int twice_density = get_n_edges(input) / get_n_cities(input);
    int auxProcBestCost = 0;
    bool received = false, received2 = false;
    bool updateCost = false;

    Path *initial_path, *initial_path_aux, *new_path[MAX_N_THREADS * 2], *path_to_send;
    Solution *sol;
    priority_queue_t *queue[MAX_N_THREADS];
    int n_cities = get_n_cities(input), whistle = -1, global_tid = 0, exit_global = 0, twice = 0, dealer = 0, flag_rec = 0, flag_rec2 = 0;

    int procBestCost = rank;

    // mpi structure name
    MPI_Datatype path_mpi;

    // number of structure members
    const int nitems = 6;
 
    // array of structure member sizes
    int blocklengths[6];
    blocklengths[0] = 1;
    blocklengths[1] = 1;
    blocklengths[2] = 1;
    blocklengths[3] = 1;
    blocklengths[4] = 1;
    blocklengths[5] = n_cities + 1;

    // structure member types
    MPI_Datatype types[6] = { MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT, MPI_LONG, MPI_INT};

    // status
    //MPI_Status status;

    // offset of structure members
    MPI_Aint offsets[6];
    offsets[0] = offsetof( Path,cost);
    offsets[1] = offsetof( Path,bound);
    offsets[2] = offsetof( Path,length);
    offsets[3] = offsetof( Path,node);
    offsets[4] = offsetof( Path,isInTour);
    offsets[5] = offsetof( Path,Tour);

    // create mpi struct
    MPI_Type_create_struct( nitems, blocklengths, offsets, types, &path_mpi);
    MPI_Type_commit( &path_mpi);


    initial_path_aux = create_path(n_cities);


    for(int i = 1; i < numprocs; i++) {
        // printf("needWork %d \n", needWork);
        workStatus = workStatus + (1 << i);
    }

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

    Path *current_path_p0;
    int flag_p0 = 0, dealer_p0 = 0;
    priority_queue_t *queue_p0;
    Path *initial_path_p0, *new_path_p0[MAX_N_PROCS];

    initial_path_p0 = create_path(n_cities);
    if (initial_path_p0 == NULL) {
        /* All needed frees and exits in error. */
        free_inputs(input);
        error();
    }


    /* Creates a queue for each thread. */
    queue_p0 = queue_create(compare);
    if (queue_p0 == NULL) {
        /* Exits in error */
        error();
    }

    set_bound(initial_path_p0, InitialLowerBound(input));

    queue_push(queue_p0, initial_path_p0);

    /* Relevant computation to secure that further ahead we
    can evenly distribute elements by the queue of each thread. */
    if (((int) (twice_density * 0.8)) < 1) {
        dealer_p0 = numprocs;
    }
    else {
        dealer_p0 = (int) (numprocs * twice_density * 0.8);
    }

    // printf("dealer_p0=%d\n",dealer_p0);
    

    /* Process (Multiple: Pop + Work) the queue that has received first element until:
    - the program ends: flag = 1 (only irrelevant elements in queue) or queue[tid]->size = 0 (no more elements)
    or
    - the queue size is large enough to evenly distribute elements (1 or 2) among the queues of the 
    other threads. This queue size and the number of distributed elements were optimized through a 
    calculation that uses the density of the graph in question.  */
    while (((int) queue_p0->size != 0) && (flag_p0 != 1) && ((int) queue_p0->size < dealer_p0)) {
        current_path_p0 = queue_pop(queue_p0);
        work(queue_p0, n_cities, &BestTourCost, input, sol, current_path_p0, &flag_p0, rank, numprocs, BestTourCostAuxx, false, NULL, &procBestCost,&updateCost, auxProcBestCost, receiver, MYTAG);
        free_path(current_path_p0);  
    }

    // printf("queue_p0->size=%ld\n", queue_p0->size);

    /* Checks if all the processing of the graph was previously done. */
    if (flag_p0 == 1 || (int) queue_p0->size == 0) {
        if (rank == 0) {
            /* Guarantees that reading and writing of one memory location is atomic. */
            exit_global += 1; 
        }
        else {
            exit_global += 1; 
            solved = 1;
        }
    }
    else {
        /* Check which thread have the elements to distribute. */
        /* Check if the number of elements to distribute allows to distribute 
        1 or 2 elements per queue thread. Then pop that elements. */
        if ((int) queue_p0->size >= numprocs) {
            for (int i = 0 ; i < numprocs; i++) {  
                new_path_p0[i] = queue_pop(queue_p0);
            }
        }
        initial_path = new_path_p0[rank];
    }

    MPI_Barrier(MPI_COMM_WORLD);
        
    /* density = edges / (2 * nodes) */

    if (rank != 0) {
        sender = rank - 1; 
    }
    else {
        sender = numprocs - 1;
    }
    if (rank == (numprocs - 1)){
        receiver = 0;
    }
    else {
        receiver = rank + 1;
    }


    /* Creates N parallel threads. All threads execute the subsequent block.
    All threads wait for each other at the end of this executing block: implicit barrier synchronization */
    /* BestTourCost is a shared variable that exists in a single location 
    and all threads can read and write it. */

    

        if (exit_global == 0) {
        #pragma omp parallel shared(BestTourCost, dealer) private(updateCost)
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

            while (!mega_global_exit){

                /* Single thread execute this region: the master thread. The 1st element goes to a single thread's queue. */
                #pragma omp master
                {
                    queue_push(queue[tid], copy_path(initial_path, n_cities));

                    //printf("%d\n", (int) queue_p0->size);
                    while ((int) queue_p0->size != 0) {
                        if((count10 % numprocs) == rank) {
                            queue_push(queue[tid], (queue_pop(queue_p0)));
                        }
                        else {
                            free_path(queue_pop(queue_p0));
                        }
                        count10++;
                    }

                    /* Relevant computation to secure that further ahead we
                    can evenly distribute elements by the queue of each thread. */
                    if (((int) (twice_density * 0.8)) < 1) {
                        dealer = omp_get_num_threads();
                    }
                    else {
                        dealer = (int) (omp_get_num_threads() * twice_density * 0.8);
                    }

                    /* Process (Multiple: Pop + Work) the queue that has received first element until:
                    - the program ends: flag = 1 (only irrelevant elements in queue) or queue[tid]->size = 0 (no more elements)
                    or
                    - the queue size is large enough to evenly distribute elements (1 or 2) among the queues of the 
                    other threads. This queue size and the number of distributed elements were optimized through a 
                    calculation that uses the density of the graph in question.  */

                    

                    BestTourCostAuxx[0] = BestTourCost;
                    BestTourCostAuxx[1] = rank;
                
                }
                while (((int) queue[tid]->size != 0) && (flag != 1) && ((int) queue[tid]->size < dealer)) {
                    // printf("Sou o proc %d e o meu receiver é o %d \n", rank, receiver);
                    if (numprocs > 1)
                    {
                        if (!received2) {
                            MPI_Irecv(BestTourCostAuxx, 2, MPI_DOUBLE, sender, MYTAG2, MPI_COMM_WORLD, &req2); 
                            // printf("Sou o proc %d e quero receber de sender %d\n", rank, sender);
                            received2 = true;
                        }

                        // printf("%lf\n", BestTourCostAuxx[0]);
                    
                        MPI_Test(&req2, &flag_rec2, &stat2);
                        // printf("flag_rec=%d\n", flag_rec);
                        if (flag_rec2) {
                            // printf("Sou o processo %d e recebi do proc %d: %lf < %lf\n", rank, stat.MPI_SOURCE, BestTourCostAuxx[0], BestTourCost);
                            if ( (BestTourCostAuxx[0] < BestTourCost) || ((BestTourCostAuxx[0] == BestTourCost) && (procBestCost > (int) BestTourCostAuxx[1]))) {
                                printf("Sou o proc %d e mandei o cost %lf\n", rank, BestTourCostAuxx[0]);
                                //atenção que temos de verificar se ele já fez o send para poder executar outro
                                MPI_Isend(BestTourCostAuxx, 2, MPI_DOUBLE, receiver, MYTAG2, MPI_COMM_WORLD, &reqsend2);
                                updateCost = true;
                                auxProcBestCost = (int) BestTourCostAuxx[1];
                                
                            }
                            received2 = false;
                        }
                    }

                    current_path = queue_pop(queue[tid]);
                    work(queue[tid], n_cities, &BestTourCost, input, sol, current_path, &flag, rank, numprocs, BestTourCostAuxx, true, &reqsend2, &procBestCost, &updateCost, auxProcBestCost, receiver, MYTAG2);
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

                    #pragma omp master
                    {
                        if (numprocs > 1) {
                            if (!received2) {
                                MPI_Irecv(BestTourCostAuxx, 2, MPI_DOUBLE, sender, MYTAG2, MPI_COMM_WORLD, &req2);
                                // printf("Sou o proc %d e quero receber de todos\n", rank);
                                received2 = true;
                            }

                            // printf("%lf\n", BestTourCostAuxx[0]);
                            
                            MPI_Test(&req2, &flag_rec2, &stat2);
                            // printf("flag_rec=%d\n", flag_rec);
                            if (flag_rec2) {
                                // printf("Sou o processo %d e recebi do proc %d: %lf < %lf\n", rank, stat.MPI_SOURCE, BestTourCostAuxx[0], BestTourCost);
                                if ((BestTourCostAuxx[0] < BestTourCost) || (BestTourCostAuxx[0] == BestTourCost && (procBestCost > (int) BestTourCostAuxx[1]))) {
                                    // printf("Sou o proc %d e dei update do cost de %lf para %lf\n", rank, BestTourCost, BestTourCostAuxx[0]);
                                    //atenção que temos de verificar se ele já fez o send para poder executar outro
                                    MPI_Isend(BestTourCostAuxx, 2, MPI_DOUBLE, receiver, MYTAG2, MPI_COMM_WORLD, &reqsend2);
                                    updateCost = true;
                                    auxProcBestCost = (int) BestTourCostAuxx[1];
                                    
                                }
                                received2 = false;
                            }

                            if (!listenSOS) {
                                MPI_Irecv(&whistlerSOS, 1, MPI_INT, sender, TAG_SOS, MPI_COMM_WORLD, &reqSOSmain);
                                // printf("Dentro: Sou o proc %d e quero receber SOS de %d\n", rank, sender);
                                listenSOS = true;
                            }

                            if (giveWork) {
                                if (queue[tid]->size > 2) {
                                    
                                    path_to_send = queue_pop(queue[tid]);
                                    if ((path_to_send->length < (n_cities* 0.75)) && (path_to_send->length >= 1)) { //alterei para >=
                                        // printf("Dentro: Sou o proc %d e vou mandar um path ao %d\n", rank, sender);
                                        // print_path(path_to_send, n_cities);
                                        MPI_Isend(path_to_send, 1, path_mpi, sender, TAG_WORK_SOS, MPI_COMM_WORLD, &reqsend11);
                                        if ((workStatus & ( 1 << sender)) == 0){
                                            workStatus = workStatus ^ (1 << sender);
                                        }
                                    }
                                    else {
                                        workStatus = workStatus & whistlerSOS;
                                        queue_push(queue[tid], path_to_send);
                                        giveWork = false;
                                    }
                                    listenSOS = false;
                                }
                            }

                            MPI_Test(&reqSOSmain, &flag_recSOSmain, &statusSOSmain);
                            if (flag_recSOSmain) {
                                workStatus = workStatus & whistlerSOS;
                                if (queue[tid]->size > 2) {
                                    
                                    path_to_send = queue_pop(queue[tid]);
                                    if ((path_to_send->length < n_cities * (0.75) ) && (path_to_send->length >= 1)) { //alterei para >=
                                        // printf("Dentro: Sou o proc %d e vou mandar um path ao %d\n", rank, sender);
                                        // print_path(path_to_send, n_cities);
                                        MPI_Isend(path_to_send, 1, path_mpi, sender, TAG_WORK_SOS, MPI_COMM_WORLD, &reqsend12);
                                        if ((workStatus & ( 1 << rank)) == 0){
                                            workStatus = workStatus ^ (1 << sender);
                                        }
                                    }
                                    else {
                                        workStatus = workStatus & whistlerSOS;
                                        queue_push(queue[tid], path_to_send);
                                        giveWork = false;
                                    }
                                    listenSOS = false;
                                }
                                else {
                                    workStatus = workStatus & whistlerSOS;
                                    giveWork = true;
                                    // printf("Dentro: Sou o proc %d e não tenho suficientes paths para mandar ao %d\n", rank, whistlerSOS);
                                }
                                
                            }
                            
                        }
                        
                        
                    }

                    current_path = queue_pop(queue[tid]);
                    work(queue[tid], n_cities, &BestTourCost, input, sol, current_path, &flag, rank, numprocs, BestTourCostAuxx, true, &reqsend, &procBestCost, &updateCost, auxProcBestCost, receiver, MYTAG);

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

                // printf("Proc %d pede trabalho ao %d\n", rank, receiver);

                #pragma omp master
                {

                    workStatus = workStatus ^ ( 1 << rank);

                    if (numprocs > 1) {
                        // printf("Sou o processo %d e fiquei sem trabalho no count = %d \n", rank, count8);
                        MPI_Isend(&workStatus, 1, MPI_INT, receiver, TAG_SOS, MPI_COMM_WORLD, &reqsend8);
                        MPI_Isend(&workStatus, 1, MPI_INT, receiver, TAG_SOSFINAL, MPI_COMM_WORLD, &reqsend9);

                        MPI_Irecv(initial_path, 1, path_mpi, receiver, TAG_WORK_SOS, MPI_COMM_WORLD, &reqSOS);
                        

                        while (!flag_recSOS && !final) {

                            MPI_Test(&reqSOS, &flag_recSOS, &statusSOS);
                            // printf("flag_recSOS = %d \n", flag_recSOS);
                            if (!flag_recSOS) {

                                if (!listenSOS2) {
                                    // printf("Proc %d preparado para receber SOSFINAL\n", rank);
                                    MPI_Irecv(&sos_received_aux, 1, MPI_INT, sender, TAG_SOSFINAL, MPI_COMM_WORLD, &reqSOS2);
                                    listenSOS2 = true;
                                }

                                MPI_Test(&reqSOS2, &flag_recSOS2, &statusSOS2);
                                if (flag_recSOS2) {
                                    // printf("WORKSTATUS = %d \n", workStatus);
                                    workStatus = workStatus & sos_received_aux;
                                    // printf("Proc %d recebe SOSFINAL de %d: WORKSTATUS = %d \n", rank, statusSOS2.MPI_SOURCE, workStatus);
                                    MPI_Isend(&workStatus, 1, MPI_INT, receiver, TAG_SOSFINAL, MPI_COMM_WORLD, &reqsend10);
                                    listenSOS2 = false;
                                }

                                if (workStatus == 0) {
                                    final = true;
                                }
                            }
                        }

                        /* Not successful receiving. */
                        if (final) { 
                            free_path(initial_path);
                            mega_global_exit = true;

                        }
                        else {
                            if ((workStatus & ( 1 << rank)) == 0){
                                workStatus = workStatus ^ ( 1 << rank);
                            }
                            // printf("SOU O %d E RECEBI O PATH \n", rank);
                            // print_path(initial_path, n_cities);
                            final = false;
                            sos_received_count = 0;
                            exit_global = 0;
                            flag_recSOS = 0;
                        }
                    }
                    else {
                        mega_global_exit = true;
                    }
            }
            
            #pragma omp barrier
            
            }


            /* Frees auxiliar structures. */
            queue_delete(queue[tid]);
            free_safe(queue[tid]);
        }
    } 

    /* Frees auxiliar structures. */
    queue_delete(queue_p0);
    free_safe(queue_p0);

    // printf("%d: cheguei com BestTourCost=%lf\n", rank, BestTourCost);

    double BestTourCostAux = 0;
    procBestCost = 0;

    if (!valid_BestTour(sol, n_cities) || get_BestTourCost(sol) != BestTourCost) {
        BestTourCost = INFINTY_DOUBLE;
    }

    BestTourCostAux = BestTourCost; 

    MPI_Barrier(MPI_COMM_WORLD);
    bool valid = true;
    for (int iii = 0; iii < numprocs; iii++) {
        //if (valid) {
            MPI_Bcast (&BestTourCostAux, 1, MPI_DOUBLE, iii, MPI_COMM_WORLD); 
            MPI_Barrier(MPI_COMM_WORLD);
            if (BestTourCostAux < BestTourCost && (int) BestTourCostAux > 0) {
                valid = false;
            }
            else if (BestTourCostAux == BestTourCost && rank > iii) {
                valid = false;
            }
            BestTourCostAux = BestTourCost; 
        //}
    }

    MPI_Finalize();

    // printf("BestCost: %lf Rank: %d Elemento: %d\n", get_BestTourCost(sol), rank, (get_BestTour(sol)[2]));
    if ((valid) && (solved != 1)) {
        
        /* Check if a valid solution was found. */
        if (valid_BestTour(sol, n_cities)) {
            return sol;
        }
        free_solution(sol);
        return NULL;
    }
    else {
        exit(0);
    }
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
 
void work(priority_queue_t *queue, int n_cities, double *BestTourCost, Inputs* input, Solution *sol, Path* current_path, int *flag, int rank, int numprocs, double *BestTourCostAuxx, bool comm, MPI_Request *reqsend, int *procBestCost, bool *updateCost, int auxProcBestCost, int receiver, int myTAG) {
    int i = 0, min_city = -1, max_city = -1;
    double newBound = 0, aux_distance = 0;
    Path *new_path;

    //printf("293: Me: %d\n", get_node(current_path));
    /* Checks if all remaining nodes in queue are worse than BestTourCost. */
    if (get_bound(current_path) >= *BestTourCost) {
        *flag = 1;
    }

    /* Tour complete, check if it is best. */
    if ((get_length(current_path) == n_cities && *flag == 0) || *updateCost) { 
        aux_distance = distance(get_node(current_path), 0, input);
        // printf("Cost current path: %lf\n", get_cost(current_path) + aux_distance);
        /* Update the best solution tour and cost. We are updating the solution that is a
        shared structure, so we need a critical region to secure that a thread waits at the 
        the beginning of that block until no other thread is executing that section. */
        #pragma omp critical
        {
            if ((BestTourCostAuxx[0] < *BestTourCost && (*updateCost)) || ((BestTourCostAuxx[0] == *BestTourCost) && (*procBestCost > (int) BestTourCostAuxx[1]) && (*updateCost))){
                // printf("Sou o proc %d e dei update do cost de %lf para %lf\n", rank, *BestTourCost, BestTourCostAuxx[0]);
                *BestTourCost = BestTourCostAuxx[0];
                *procBestCost = auxProcBestCost;
                *updateCost = false;
            }
            if(get_length(current_path) == n_cities && *flag == 0){
                if (get_cost(current_path) + aux_distance < *BestTourCost && aux_distance >= 0) {
                    set_BestTour(sol, get_Tour(current_path), n_cities);
                    set_BestTour_item(n_cities, sol, 0, n_cities);
                    *BestTourCost = get_cost(current_path) + aux_distance;
                    BestTourCostAuxx[0] = *BestTourCost;
                    BestTourCostAuxx[1] = (double) rank;
                    set_BestTourCost(sol, *BestTourCost);
                    *procBestCost = rank;
                    if (comm && numprocs > 1) {

                        MPI_Isend(BestTourCostAuxx, 2, MPI_DOUBLE, receiver, myTAG, MPI_COMM_WORLD, reqsend);
                        // printf("Mandei %lf do processo %d para o %d\n", BestTourCostAuxx[0], rank, receiver);
                    }
                }
            }
            
            

            
        }
    }
    
    else if (*flag == 0) {
        min_city = get_mincity(get_node(current_path), input);
        max_city = get_maxcity(get_node(current_path), input) + 1;
        
        for (i = min_city; i < max_city; i++) {
            /* Connection does not exist. */
            if (distance(get_node(current_path), i, input) < 0) {
                continue;
            }

            /* Check if node is already in tour */
            if ((1 << i) & get_isInTour(current_path) && !(i == 0 && get_length(current_path) == n_cities)) {
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

            set_isInTour(new_path, i, get_isInTour(current_path));

            /* Insert it in queue. */
            queue_push(queue, new_path);
        }  
    }

    return;
}