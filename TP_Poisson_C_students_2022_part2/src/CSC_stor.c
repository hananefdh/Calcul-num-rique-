#include "lib_poisson1D.h"
#include <stdlib.h>
#include <stdio.h>

void poisson1D_CSC(int la, double *values, int *row_indices, int *col_ptr)
{
    int nnz = 3 * la - 2; // Nombre total de valeurs non nulles
    int idx = 0;
    col_ptr[0] = 0;

    for (int i = 0; i < la; ++i)
    {
        if (i > 0)
        {
            values[idx] = -1.0;       // Valeur de la diagonale inférieure
            row_indices[idx] = i - 1; // Indice de ligne correspondant
            ++idx;
        }

        values[idx] = 2.0;    // Valeur de la diagonale principale
        row_indices[idx] = i; // Indice de ligne correspondant
        ++idx;

        if (i < la - 1)
        {
            values[idx] = -1.0;       // Valeur de la diagonale supérieure
            row_indices[idx] = i + 1; // Indice de ligne correspondant
            ++idx;
        }

        col_ptr[i + 1] = idx; // Pointeur de colonne pour la prochaine colonne
    }
}
// Produit matrice-vecteur pour CSC
void dcscmv(int la, double *values, int *row_indices, int *col_ptr, double *x, double *y)
{
    for (int j = 0; j < la; ++j)
    {
        y[j] = 0.0;
        for (int i = col_ptr[j]; i < col_ptr[j + 1]; ++i)
        {
            y[row_indices[i]] += values[i] * x[j];
        }
    }
}