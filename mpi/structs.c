/********************************************************************************
* (c) Tomás Fonseca - 66325
*     Afonso Alemão - 96135
*     Rui Daniel    - 96317
*
* Last Modified: 12/04/2023
*
* Name: structs.c
*
* Description: Manages the auxiliary structures of information.
*
* Comments:
*
**********************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "structs.h"

struct _solution {
    int *BestTour;
	double BestTourCost;
};

struct _inputs {
    double **adj_matrix;
    double max_value;
	int n_cities;
	int n_edges;
    double *min1;
    double *min2;
    int *max_city;
    int *min_city;
};

struct _path {
    double cost;
    double bound;
	int length;
	int node;
    long isInTour;
    int Tour[];
};


/**********************************************************************************
* create_inputs()
*
* Arguments:    n_cities - number of cities of the graph
*
* Returns:      Inputs pointer
*
* Side-Effects: Allocates memory
*
* Description:  Creates a inputs element
*
**********************************************************************************/

Inputs *create_inputs(int n_cities) {
    Inputs *new_input;
	new_input = NULL;
    int j = 0;

    new_input = (Inputs *) malloc(sizeof(Inputs));
    if (new_input == NULL) error();

    new_input->adj_matrix = aloc2d(n_cities, n_cities);
    if (new_input->adj_matrix == NULL) {
        /* All needed frees and exits in error */
        free_inputs(new_input);
        error();
    }

    new_input->max_value = 0;
	new_input->n_cities = n_cities;

    new_input->min1 = (double *) malloc((n_cities) * sizeof(double));
    if (new_input->min1 == NULL)  {
        /* All needed frees and exits in error */
        free_inputs(new_input);
        error();
    }

    new_input->min2 = (double *) malloc((n_cities) * sizeof(double));
    if (new_input->min2 == NULL)  {
        /* All needed frees and exits in error */
        free_inputs(new_input);
        error();
    }

    new_input->max_city = (int *) malloc((n_cities) * sizeof(int));
    if (new_input->max_city == NULL)  {
        /* All needed frees and exits in error */
        free_inputs(new_input);
        error();
    }

    new_input->min_city = (int *) malloc((n_cities) * sizeof(int));
    if (new_input->min_city== NULL)  {
        /* All needed frees and exits in error */
        free_inputs(new_input);
        error();
    }

    for (int i = 0; i < n_cities; i++) {
        new_input->min1[i] = -1;
        new_input->min2[i] = -1;
        new_input->min_city[i] = -1;
        new_input->max_city[i] = -1;
		for (j = 0; j < n_cities; j++) {
            new_input->adj_matrix[i][j] = -1;
        }
	}

    return new_input;
}

/********************************************************************************
* Getters of Inputs
********************************************************************************/

double get_ajd_matrix_item(int row, int col, Inputs *got_input) {
    if (got_input == NULL) return -1;
    if (row < 0 || row >= got_input->n_cities || col < 0 || col >= got_input->n_cities) {
        free_inputs(got_input);
        error();
    }

    return got_input->adj_matrix[row][col];
}

int get_mincity(int index, Inputs *got_input) {
    if (got_input == NULL) return -1;
    if (index < 0 || index >= got_input->n_cities) {
        free_inputs(got_input);
        error();
    }

    return got_input->min_city[index];
}

int get_maxcity(int index, Inputs *got_input) {
    if (got_input == NULL) return -1;
    if (index < 0 || index >= got_input->n_cities) {
        free_inputs(got_input);
        error();
    }

    return got_input->max_city[index];
}

void print_minmaxcity(Inputs *got_input, int n_cities) {
    if (got_input == NULL) return;
    for (int i = 0; i < n_cities; i++) {
        printf("Node = %d,\t min_city = %d, max_city = %d\n", i, got_input->max_city[i], got_input->min_city[i]);
    }
    printf("\n");
    return;
}

int get_n_cities(Inputs *got_input) {
    if (got_input == NULL) return -1;
    return got_input->n_cities;
}

double get_max_value(Inputs *got_input) {
    if (got_input == NULL) return -1;
    return got_input->max_value;
}

int get_n_edges(Inputs *got_input) {
    if (got_input == NULL) return -1;
    return got_input->n_edges;
}

double get_min1(int index, Inputs *got_input) {
    if (got_input == NULL) return -1;
    if (index < 0 || index >= got_input->n_cities) {
        free_inputs(got_input);
        error();
    }
    return got_input->min1[index];
}

double get_min2(int index, Inputs *got_input) {
    if (got_input == NULL) return -1;
    if (index < 0 || index >= got_input->n_cities) {
        free_inputs(got_input);
        error();
    }
    return got_input->min2[index];
}

/********************************************************************************
* Setters of Inputs
********************************************************************************/

void set_ajd_matrix_item(int row, int col, Inputs *got_input, double cost) {
    if (got_input == NULL) return;
    if (row < 0 || row >= got_input->n_cities || col < 0 || col >= got_input->n_cities) {
        free_inputs(got_input);
        error();
    }
    got_input->adj_matrix[row][col] = cost;

    /* Updates min1 and min2 of both the nodes in the new connection. 
    We have bidirectional graph, so we do the update only once. */
    if (row < col) {
        if (got_input->min1[row] == -1) {
            got_input->min1[row] = cost;
        }
        else if (cost <= got_input->min1[row]) {
            got_input->min2[row] = got_input->min1[row];
            got_input->min1[row] = cost;
        }
        else if (cost < got_input->min2[row] || got_input->min2[row] == -1) {
            got_input->min2[row] = cost;
        }

        if (got_input->min1[col] == -1) {
            got_input->min1[col] = cost;
        }
        else if (cost <= got_input->min1[col]) {
            got_input->min2[col] = got_input->min1[col];
            got_input->min1[col] = cost;
        }
        else if (cost < got_input->min2[col] || got_input->min2[col] == -1) {
            got_input->min2[col] = cost;
        }
    }

    if (got_input->max_city[row] == -1) {
        got_input->max_city[row] = col;
        got_input->min_city[row] = col;
    }
    else if (col > got_input->max_city[row]) {
        got_input->max_city[row] = col;
    }
    else if (col < got_input->min_city[row]) {
        got_input->min_city[row] = col;
    }


    if (got_input->max_city[col] == -1) {
        got_input->max_city[col] = row;
        got_input->min_city[col] = row;
    }
    else if (row > got_input->max_city[col]) {
        got_input->max_city[col] = row;
    }
    else if (row < got_input->min_city[col]) {
        got_input->min_city[col] = row;
    }
    return;
}

void set_mincity(int index, int city, Inputs *got_input) {
    if (got_input == NULL) return;
    if (index < 0 || index >= got_input->n_cities) {
        free_inputs(got_input);
        error();
    }
    got_input->min_city[index] = city;
    return;
}

void set_maxcity(int index, int city, Inputs *got_input) {
    if (got_input == NULL) return;
    if (index < 0 || index >= got_input->n_cities) {
        free_inputs(got_input);
        error();
    }
    got_input->max_city[index] = city;
    return;
}

void set_n_cities(Inputs *got_input, int n_cities) {
    if (got_input == NULL) return;
    got_input->n_cities = n_cities;
	return;
}

void set_max_value(Inputs *got_input, double max_value) {
    if (got_input == NULL) return;
    got_input->max_value = max_value;
	return;
}

void set_n_edges(Inputs *got_input, int n_edges) {
    if (got_input == NULL) return;
    got_input->n_edges = n_edges;
	return;
}

/********************************************************************************
* Free Inputs
********************************************************************************/

void free_inputs(Inputs *got_input) {
    if (got_input != NULL) {
        free2d(got_input->adj_matrix, got_input->n_cities);
        free_safe(got_input->min1);    
        free_safe(got_input->min2);
        free_safe(got_input);
    }
    
    return;
}

/**********************************************************************************
* create_solution()
*
* Arguments:    n_cities - number of cities of the graph
*
* Returns:      Solution pointer
*
* Side-Effects: Allocates memory
*
* Description:  Creates a solution
*
**********************************************************************************/

Solution *create_solution(int n_cities) {
    Solution *new_solution;
	new_solution = NULL;

    new_solution = (Solution *) malloc(sizeof(Solution));
    if (new_solution == NULL) return NULL;

    new_solution->BestTour = (int *) malloc((n_cities + 1) * sizeof(int));
    if (new_solution->BestTour == NULL) {
        /* All needed frees and exits in error */
        free_solution(new_solution);
        return NULL;
    }

	for (int i = 0; i < n_cities + 1; i++) {
		new_solution->BestTour[i] = -1;
	}

	new_solution->BestTourCost = -1;
    return new_solution;
}


/********************************************************************************
* Getters of Solution
********************************************************************************/

int get_BestTour_item(int index, Solution *got_solution, int n_cities) {
    if (got_solution == NULL) return -1;
    if (index < 0 || index >= n_cities + 1) return -1;
    return got_solution->BestTour[index];
}

int *get_BestTour(Solution *got_solution) {
    if (got_solution == NULL) return NULL;
    return got_solution->BestTour;
}

double get_BestTourCost(Solution *got_solution) {
    if (got_solution == NULL) return -1;
    if (got_solution == NULL) {
        return -1;
    }
    return got_solution->BestTourCost;
}

int valid_BestTour(Solution *got_solution, int n_cities) {
    if (got_solution == NULL) return -1;
    if (got_solution->BestTour[n_cities] == -1) {
        return 0;
    }
    return 1;
}


/********************************************************************************
* Setters of Solution
********************************************************************************/

void set_BestTour_item(int index, Solution *got_solution, int city, int n_cities) {
    if (got_solution == NULL) return;
    if (index < 0 || index >= n_cities + 1) return;
    got_solution->BestTour[index] = city;
    return;
}

void set_BestTour(Solution *got_solution, int *tour, int n_cities) {
    if (got_solution == NULL) return;
    for (int i = 0; i < n_cities + 1; i++) {
        set_BestTour_item(i, got_solution, tour[i], n_cities);
    }
    return;
}

void set_BestTourCost(Solution *got_solution, double BestTourCost) {
    if (got_solution == NULL) return;
    got_solution->BestTourCost = BestTourCost;
	return;
}

/********************************************************************************
* Free Solution
********************************************************************************/

void free_solution(Solution *got_solution) {
    if (got_solution != NULL) {
        free_safe(got_solution->BestTour);
        free_safe(got_solution);
    }
    return;
}


/**********************************************************************************
* create_path()
*
* Arguments:    n_cities - number of cities of the graph
*
* Returns:      Path pointer
*
* Side-Effects: Allocates memory
*
* Description:  Creates a path
*
**********************************************************************************/

Path *create_path(int n_cities) {
    if (n_cities <= 0) return NULL;
    Path *new_path;
	new_path = NULL;

    new_path = (Path *) malloc(sizeof(Path) + (n_cities + 1) * sizeof(int));
    if (new_path == NULL) return NULL;

    // new_path->Tour = (int *) malloc((n_cities + 1) * sizeof(int));
    // if (new_path->Tour == NULL) {
    //     /* All needed frees and exits in error */
    //     free_path(new_path);
    //     return NULL;
    // }

    new_path->Tour[0] = 0;
    
	for (int i = 1; i < n_cities + 1; i++) {
		new_path->Tour[i] = -1;
	}

	new_path->cost = 0;
    new_path->bound = 0;
    new_path->length = 1;
    new_path->node = 0;
    new_path->isInTour = 1;
    return new_path;
}

/********************************************************************************
* Getters of Path
********************************************************************************/

int get_Tour_item(int index, Path *got_path, int n_cities) {
    if (got_path == NULL) return -1;
    if (index < 0 || index >= n_cities + 1) return -1;
    return got_path->Tour[index];
}

int *get_Tour(Path *got_path) {
    if (got_path == NULL) return NULL;
    return got_path->Tour;
}

double get_cost(Path *got_path) {
    if (got_path == NULL) return -1;
    return got_path->cost;
}

double get_bound(Path *got_path) {
    if (got_path == NULL) return -1;
    return got_path->bound;
}

int get_length(Path *got_path) {
    if (got_path == NULL) return -1;
    return got_path->length;
}

int get_node(Path *got_path) {
    if (got_path == NULL) return -1;
    return got_path->node;
}

long get_isInTour(Path *got_path) {
    if (got_path == NULL) return -1;
    return got_path->isInTour;
}

void print_path(Path *got_path, int n_cities) {
    if (got_path == NULL) return;
    printf("Node = %d,\tLength = %d,\tBound = %lf,\tCost = %lf,\tTour = ", got_path->node, got_path->length, got_path->bound, got_path->cost);
    for (int i = 0; i < n_cities + 1; i++) {
         printf("%d ", got_path->Tour[i]);
    }
    printf("\n");
    return;
}


/********************************************************************************
* Setters of Path
********************************************************************************/

void set_Tour_item(int index, Path *got_path, int city, int n_cities) {
    if (got_path == NULL) return;
    if (index < 0 || index >= n_cities + 1) return;
    got_path->Tour[index] = city;
    return;
}

void set_Tour(Path *got_path, int *tour, int n_cities) {
    if (got_path == NULL) return;
    
    for (int i = 0; i < n_cities + 1; i++) {
        set_Tour_item(i, got_path, tour[i], n_cities);
    }
    return;
}

void set_cost(Path *got_path, double cost) {
    if (got_path == NULL) return;
    got_path->cost = cost;
    return;
}

void set_bound(Path *got_path, double bound) {
    if (got_path == NULL) return;
    got_path->bound = bound;
    return;
}

void set_length(Path *got_path, int length) {
    if (got_path == NULL) return;
    got_path->length = length;
    return;
}

void set_node(Path *got_path, int node) {
    if (got_path == NULL) return;
    got_path->node = node;
    return;
}

void set_isInTour(Path *got_path, int index, long isInTour) {
    if (got_path == NULL) return;
    got_path->isInTour = isInTour | (1 << index);

    return;
}

/********************************************************************************
* Free Path
********************************************************************************/

void free_path(Path *got_path) {
    if (got_path != NULL) {
        // free_safe(got_path->Tour);
        free_safe(got_path);
    }
    return;
}
