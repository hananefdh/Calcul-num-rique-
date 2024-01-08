

#include "lib_poisson1D.h"
#include <stdlib.h>
#include <stdio.h>

// Stockage CSR
void poisson1D_CSR(int la, double *values, int *col_indices, int *row_ptr)
{
    int nnz = 3 * la - 2; // Nombre total de valeurs non nulles
    int idx = 0;
    row_ptr[0] = 0;

    for (int i = 0; i < la; ++i)
    {
        if (i > 0)
        {
            values[idx] = -1.0;       // Valeur de la diagonale inférieure
            col_indices[idx] = i - 1; // Indice de colonne correspondant
            ++idx;
        }

        values[idx] = 2.0;    // Valeur de la diagonale principale
        col_indices[idx] = i; // Indice de colonne correspondant
        ++idx;

        if (i < la - 1)
        {
            values[idx] = -1.0;       // Valeur de la diagonale supérieure
            col_indices[idx] = i + 1; // Indice de colonne correspondant
            ++idx;
        }

        row_ptr[i + 1] = idx; // Pointeur de ligne pour la prochaine ligne
    }
}

void dcsrmv(int la, double *values, int *col_indices, int *row_ptr, double *x, double *y)
{
    for (int i = 0; i < la; ++i)
    {
        double sum = 0.0;
        for (int j = row_ptr[i]; j < row_ptr[i + 1]; ++j)
        {
            sum += values[j] * x[col_indices[j]];
        }
        y[i] = sum;
    }
}
