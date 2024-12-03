////////
// title: PhyiCell/custom_modules/custom.cpp
//
// language: C++
// date: 2015-2024
// license: BSD-3-Clause
// author: Alexandre Bertin, Marcelo Hurtado, Elmar Bucher, Paul Macklin
// original source code: https://github.com/MathCancer/PhysiCell
// modified source code: https://github.com/elmbeech/physicellembedding
// modified source code: https://github.com/VeraPancaldiLab/RL_TME
////////


#include "custom.h"

int get_celltypescount() {
    // cell count per celltype
    parameters.ints("count_nurse_cell") = 0;
    parameters.ints("count_cancer_cell") = 0;

    for (int i=0; i < all_cells->size(); i++) {
        Cell *pCell = (*all_cells)[i];
        if (!(pCell->phenotype.death.dead)) {
            if (pCell->type_name == "cancer_cell") {
                parameters.ints("count_cancer_cell")++;  // hard coded
            }
            else if (pCell->type_name == "nurse_cell") {
                parameters.ints("count_nurse_cell")++;  // hard coded
            }
            else {
                printf("Error: unknown cell type variable! %s", pCell->type_name.c_str());
            }
        }
    }
    return 0;
}
