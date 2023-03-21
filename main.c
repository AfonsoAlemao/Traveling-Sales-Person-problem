/**********************************************************************************
* (c) Tomás Fonseca - 66325
*     Afonso Alemão - 96135
*     Rui Daniel    - 96317
*
* Last Modified: 23/03/2023
*
* Name: main.c
*
* Description:	Main program to the submition of the CPD Project 22/23.
*				input file: <name>.in
*
* Comments:
*
**********************************************************************************/
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "queue.h"
#include "tsp-omp.h"
#include "auxiliar.h"

/**********************************************************************************
* main()
*
* Arguments:    argc - number of arguments.
*				argv - pointer to the arguments table.
*
* Returns:      int
*
* Side-Effects: Closes the program at the end.
*
* Description:  Main program.
*
********************************************************************************/

int main(int argc, char *argv[]) {
    double exec_time;
    Inputs *input;
    Solution *solution;

    input = parse_inputs(argc, argv);

    exec_time = -omp_get_wtime();

    solution = tsp_omp(input); 

    exec_time += omp_get_wtime();
    fprintf(stderr, "%.1fs\n", exec_time);

    print_result(solution, get_n_cities(input));

    free_solution(solution);
    free_inputs(input);
    return 0;
}



