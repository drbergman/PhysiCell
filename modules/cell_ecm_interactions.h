#ifndef _cell_ecm_interactions_
#define _cell_ecm_interactions_

#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h"
#include "./extracellular_matrix.h"

using namespace BioFVM; 
using namespace PhysiCell;

double sign_function(double number);

void custom_update_cell_velocity(Cell *pCell, Phenotype &phenotype, double dt);

void ecm_to_cell_interactions(Cell *pCell, Phenotype &phenotype, double dt);

// uses cell motility vector for realigning ECM.
void ecm_remodeling_function(Cell *pCell, Phenotype &phenotype, double dt);

void custom_update_motility_vector(Cell *pCell, Phenotype &phenotype, double dt_);

void create_default_ecm_compatible_agent(void);

void SVG_plot_custom(std::string filename, Microenvironment &M, double z_slice, double time, std::vector<std::string> (*cell_coloring_function)(Cell *), std::string line_pattern);
void write_ecm_Data_matlab(std::string filename);

#endif