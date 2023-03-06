/**********************************************************************************
* (h) Tomás Fonseca - 66325
*     Afonso Alemão - 96135
*     Rui Daniel    - 96317
*
* Last Modified: 27/02/2023
*
*
* Name: structs.h
*
* Description: structs.c header
*
* Comments:
*
**********************************************************************************/

#ifndef _STRUCTS_H
#define _STRUCTS_H

/* Structs */
typedef struct _inputs Inputs;
typedef struct _solution Solution;
typedef struct _path Path;

/* Initializers */
Inputs *create_inputs(int);
Solution *create_solution(int);
Path *create_path(int);

/* Frees */
void free_inputs(Inputs *);
void free_solution(Solution *);
void free_path(Path *got_path);

/* Getters for struct Inputs */
double get_ajd_matrix_item(int, int, Inputs *);
double get_max_value(Inputs *);
int get_n_edges(Inputs *got_input);
int get_n_cities(Inputs *);
double get_min1(int index, Inputs *got_input);
double get_min2(int index, Inputs *got_input);

/* Getters for struct Solution */
int get_BestTour_item(int, Solution *, int);
int *get_BestTour(Solution *);
double get_BestTourCost(Solution *);
int valid_BestTour(Solution *, int);

/* Getters for struct Path */
int get_Tour_item(int, Path *, int);
int *get_Tour(Path *got_path);
double get_cost(Path *);
double get_bound(Path *);
int get_length(Path *);
int get_node(Path *);
void print_path(Path *got_path, int n_cities);

/* Setters for struct Inputs */
void set_ajd_matrix_item(int, int, Inputs *, double);
void set_max_value(Inputs *, double);
void set_n_cities(Inputs *, int);
void set_n_edges(Inputs *got_input, int n_edges);

/* Setters for struct Solution */
void set_BestTour_item(int, Solution *, int, int);
void set_Tour(Path *got_path, int *tour, int n_cities);
void set_BestTour(Solution *, int *, int);
void set_BestTourCost(Solution *, double);

/* Setters for struct Path */
void set_Tour_item(int, Path *, int, int);
void set_Tour(Path *, int *, int);
void set_cost(Path *, double);
void set_bound(Path *, double);
void set_length(Path *, int);
void set_node(Path *, int);

#endif