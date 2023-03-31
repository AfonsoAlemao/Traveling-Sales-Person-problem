/********************************************************************************
* (c) Tomás Fonseca - 66325
*     Afonso Alemão - 96135
*     Rui Daniel    - 96317
*
* Last Modified: 23/03/2023
*
* Name: auxiliar.c
*
* Description: Manages the auxiliar methods of our program
*
* Comments:
*
**********************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "auxiliar.h"

#define SWAP(x, y) void* tmp = x; x = y; y = tmp;

/**********************************************************************************
* parse_inputs()
*
* Arguments:    argc - number of arguments.
*				argv - pointer to the arguments array.
*
* Returns:      Inputs pointer
*
* Side-Effects: Allocates memory
*
* Description:  Creates an Inputs element based on the program input arguments
*
**********************************************************************************/

Inputs *parse_inputs(int argc, char *argv[]) {
	FILE *fpIn;
	Inputs *inputs;
	double max_value;

    if (argc != 3) error();

    max_value = atof(argv[2]);

    /* Validate and open the input file */
	ValidInputFileName(argv[1]);

	fpIn = OpenFile(argv[1], "r");
    inputs = ReadFileIn(fpIn, max_value);
	fclose(fpIn);

	return inputs;
}


/******************************************************************************
 * OpenFile()
 *
 * Arguments: name - string with the name of the file to be opened
 * 			  mode - string with the opening mode of the file
 *
 * Returns:  Pointer to the open file
 * Side-Effects: (none)
 *
 * Description:
 *   Opens received file, and returns pointer to it
 *****************************************************************************/

FILE *OpenFile(char *name, char *mode) {
	if (name == NULL || mode == NULL) return NULL;

	FILE *fp;

	fp = fopen(name, mode);
	if (fp == NULL) {
		error();
	}
	return fp;
}


/******************************************************************************
 * ValidInputFileName()
 *
 * Arguments: name - string with the name of the file to be validated
 *
 * Returns: (void)
 *
 * Side-Effects: (none)
 *
 * Description:
 *   Checks if the input file name is valid
 *****************************************************************************/

void ValidInputFileName(char *name) {
	if (name == NULL) return;

    char *extension = ".in";
    int i = 0;
	int file_size = 0;

    file_size = strlen(name);

    /* File without extension */
    if (file_size < (int) strlen(extension)) {
		error();
	}
	/* Check the extension */
	for (i = 0; i < (int) strlen(extension); i++) {
		if (extension[i] != name[file_size - strlen(extension) + i]) {
			error();
		}
	}
	return;
}


/*************************************************************************************************************
 * ReadFileIn()
 *
 * Arguments: fp - pointer to file that will be read
 *			  max_value - double with the maximum value that we accept for the solution
 *
 * Returns:  Inputs pointer
 *
 * Side-Effects: (none)
 *
 * Description:
 * 	 Collects problem data in the input file
 ************************************************************************************************************/

Inputs *ReadFileIn(FILE *fp, double max_value) {
	if (fp == NULL) return NULL;

	int v1 = 0, v2 = 0;
	int counter = 0;
	double distance = 0;
	Inputs *inputs;
	int n_cities = 0, n_edges = 0;

	/* Read problem header */
	if (fscanf(fp, "%d %d", &n_cities, &n_edges) == EOF)
		return NULL;

	/* Creation of the Inputs structure, memory allocation and filling */
	inputs = create_inputs(n_cities);

	set_n_edges(inputs, n_edges);
	set_max_value(inputs, max_value);

	while (fscanf(fp, "%d %d %lf", &v1, &v2, &distance) != EOF) {
		counter++;
		set_ajd_matrix_item(v1, v2, inputs, distance);
		set_ajd_matrix_item(v2, v1, inputs, distance);
	}

	if (counter != n_edges) {
		free_inputs(inputs);
		error();
	}

	/* In this case it is a disconnected graph */
	if (n_edges < n_cities) {
		return NULL;
	}

	return inputs;
}

/******************************************************************************
 * print_matrix()
 *
 * Arguments: n_cities - number of cities of the graph
 * 			  inputs - input data
 *
 * Returns: (void)
 *
 * Side-Effects: (none)
 *
 * Description: Prints to stdout the adjacent matrix that represents 
 * 				the cities graph
 *****************************************************************************/

void print_matrix(int n_cities, Inputs *inputs) {
	if (inputs == NULL) return;
	printf("\t");
	for (int i = 0; i <  n_cities; i++) {
		printf("%d\t", i);
	}
	printf("\n");
	for (int i = 0; i <  n_cities; i++) {
		printf("%d\t", i);
		for (int j = 0; j < n_cities; j++){
			printf("%.1lf\t", get_ajd_matrix_item(i, j, inputs));
		}
		printf("\n");
		
	}
}

/******************************************************************************
 * error()
 *
 * Arguments: (none)
 *
 * Returns: (void)
 *
 * Side-Effects: Exits program
 *
 * Description:
 *   If an error occurs, it silently ends the program
 *****************************************************************************/

void error() {
    exit(0);
}


/******************************************************************************
 * aloc2d()
 *
 * Arguments: rows, cols - dimensions of the array to allocate
 *
 * Returns: Pointer to double*
 *
 * Side-Effects: Allocates memory
 *
 * Description:
 *   Allocation of memory for two-dimensional array of doubles
 *****************************************************************************/

double **aloc2d(int rows, int cols) {
    int i, j;
    double **array;
    if ((array = (double **) calloc(rows, sizeof(double *))) == NULL) {
        return NULL;
    }

    for (i = 0; i < rows; i++) {
        if ((array[i] = (double *) calloc(cols, sizeof(double))) == NULL) {
			for (j = 0; j < i; j++) {
				free_safe(array[j]);
			}
			free_safe(array);
            return NULL;
        }
    }
    return array;
}


/******************************************************************************
 * free2d()
 *
 * Arguments: array - two-dimensional array to be freed
 * 			   rows - size of the array to be freed
 *
 * Returns: (void)
 *
 * Side-Effects: frees memory for 2d array of double
 *
 * Description:
 *   Free memory for two-dimensional array of double
 *****************************************************************************/

void free2d(double **array, int rows) {
    int i;
    for (i = 0; i < rows; i++) {
		free_safe(array[i]);
    }
	free_safe(array);
    return;
}


/**********************************************************************************
* print_result()
*
* Arguments:    n_cities - number of cities of the graph
				solution - pointer to the solution
*
* Returns:      (void)
*
* Side-Effects: (none)
*
* Description:
* 	Prints solution to stdout
**********************************************************************************/

void print_result(Solution *solution, int n_cities) {
	if (solution == NULL || get_BestTourCost(solution) < 0) {
		fprintf(stdout, "NO SOLUTION\n");
	}
	else {
		fprintf(stdout, "%.1lf\n", get_BestTourCost(solution));

		for (int i = 0; i < n_cities + 1; i++) {
			fprintf(stdout, "%d", get_BestTour_item(i, solution, n_cities));
			if (i != n_cities) {
				fprintf(stdout, " ");
			}
			else {
				fprintf(stdout, "\n");
			}
		}
	}
}


/**********************************************************************************
* distance()
*
* Arguments:    row, col - index of the cities
*				got_input - input data
*
* Returns:      Double with the distance between two cities
*
* Side-Effects: (none)
*
* Description:
*	 Gets the distance between two cities
**********************************************************************************/

double distance(int row, int col, Inputs *got_input) {
    return get_ajd_matrix_item(row, col, got_input);
}


/**********************************************************************************
* InitialLowerBound()
*
* Arguments:    inputs - pointer to the input
*
* Returns:      double with the initial lower bound
*
* Side-Effects: (none)
*
* Description:
*	 Computes the initial lower bound
**********************************************************************************/

double InitialLowerBound(Inputs *inputs) {
	if (inputs == NULL) return -1;
	double result = 0;

	for (int i = 0; i < get_n_cities(inputs); i++) {
		result += get_min1(i, inputs) + get_min2(i, inputs);
	}

	return result / 2;
}


/**********************************************************************************
* newLowerBound()
*
* Arguments:    inputs - pointer to the input
*				path - pointer to the current path
*				t - index of destination city
*
* Returns:      double with the new lower bound
*
* Side-Effects: (none)
*
* Description:
*	 Computes the new lower bound
**********************************************************************************/

double newLowerBound(Inputs *inputs, Path *path, int t) {
	if (inputs == NULL || path == NULL) return -1;
	double dist = 0, cf = 0, ct = 0, LB = get_bound(path);
	int f = 0;

	f = get_node(path);

	dist = distance(f, t, inputs);

	if (dist >= get_min2(f, inputs)) {
		cf = get_min2(f, inputs);
	}
	else {
		cf = get_min1(f, inputs);
	}

	if (dist >= get_min2(t, inputs)) {
		ct = get_min2(t, inputs);
	}
	else {
		ct = get_min1(t, inputs);
	}

	return LB + dist - (cf + ct) / 2;
}

/**********************************************************************************
* compare()
*
* Arguments:    a, b - pointers to the elements to compare
*
* Returns:      char with the compare result (e1 < e2 ? 1 : 0)
*
* Side-Effects: (none)
*
* Description:
*	 Compares two elements based on their bound. If there is a tie, compares
*	the elements bases on their index.
**********************************************************************************/

char compare (void *a, void *b) {
	if (a == NULL || b == NULL) return -1;
	Path *e1 = (Path *) a;
	Path *e2 = (Path *) b;
	double b1 = get_bound(e1), b2 = get_bound(e2);
	if (b1 < b2) {
		return 0;
	}
	else if (b1 == b2) {
		if (get_node(e1) <= get_node(e2)) {
			return 0;
		}
		else {
			return 1;
		}
	}
	return 1;
}

/**********************************************************************************
 * InitializeIsInTour()
 *
 * Arguments: isInTour - tour to be initialized 
 *			  n_cities - number of cities of the graph
 * 
 * Return: (void) 
 *
 * Side effects: (none)
 *
 * Description: Initializes isInTour that marks the elements in the current tour
 *
 ********************************************************************************/

void InitializeIsInTour(bool *isInTour, int n_cities) {
	if (isInTour == NULL) return;

	for (int i = 0; i < n_cities; i++) {
		isInTour[i] = true;
	}
	return;
}

/**********************************************************************************
 * free_safe()
 *
 * Arguments: aux - pointer to the element to be freed (void *)
 *
 * Return: (void *) NULL
 *
 * Side effects: frees element
 *
 * Description: Free element
 *
 ********************************************************************************/

void *free_safe (void *aux) {
    if (aux != NULL) {
        free(aux);
    }
    return NULL;
}