	/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2025, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/
 
#include <vector>
#include <string>
#include <functional>
#include <memory>
#include <algorithm>

#ifndef __PhysiCell_signal_response__
#define __PhysiCell_signal_response__

#include "./PhysiCell_constants.h" 
#include "./PhysiCell_phenotype.h" 
#include "./PhysiCell_cell.h" 

namespace PhysiCell{

// create the signal and behavior dictionaries 
void setup_signal_behavior_dictionaries( void );

using SignalValue = std::function<double(Cell *)>;
struct Signal
{
	private:
	std::string name; // primary synonym
	SignalValue value;
	
	public:
	std::string get_name() const { return name; }
	Signal() : name("no name"), value(nullptr) {}
	Signal(std::string my_name, SignalValue my_value) : name(std::move(my_name)), value(std::move(my_value)) {}
	
	double get_value(Cell* pCell) const { return value(pCell); }
};

using BehaviorValue = std::function<double&(Cell *)>;
using BehaviorBaseValue = std::function<double(Cell_Definition *)>;
struct Behavior
{
private:
	std::string name; // primary synonym
	BehaviorValue value;
	BehaviorBaseValue base_value;

public:
	std::string get_name() const { return name; }
	Behavior() : name("no name"), value(nullptr), base_value(nullptr) {}
	Behavior(std::string my_name, BehaviorValue my_value, BehaviorBaseValue my_base_value) : name(std::move(my_name)), value(std::move(my_value)), base_value(std::move(my_base_value)) {}

	double get_value(Cell* pCell) const { return value(pCell); }
	void set_value(Cell *pCell, double new_value) const { value(pCell) = new_value; }
	double get_base_value(Cell *) const;
	double get_base_value(Cell_Definition *pCD) const { return base_value(pCD); }
};

double substrate_density(Cell *pCell, int substrate_index);
double internalized_substrate_density(Cell *pCell, int substrate_index);
double substrate_gradient_norm(Cell *pCell, int substrate_index);
double cell_pressure(Cell *pCell);
double cell_volume(Cell *pCell);
double cell_contact(Cell *, int);
double live_cell_contact(Cell *);
double dead_cell_contact(Cell *);
double apoptotic_cell_contact(Cell *);
double necrotic_cell_contact(Cell *);
double other_dead_cell_contact(Cell *);
double contact_with_basement_membrane(Cell *);
double cell_damage(Cell *);
double damage_delivered(Cell *);
double is_attacking(Cell *);
double is_dead(Cell *);
double total_attack_time(Cell *);
double get_current_time(Cell *);
double is_apoptotic(Cell *);
double is_necrotic(Cell *);
double cell_custom_signal(Cell *, int);

double &cell_secretion_rate(Cell *, int);
double &cell_secretion_target(Cell *, int);
double &cell_uptake_rate(Cell *, int);
double &cell_net_export_rate(Cell *, int);
double &cell_cycle_entry_rate(Cell *);
double &phase_exit_rate(Cell *, int);
double &cell_apoptosis_rate(Cell *);
double &cell_necrosis_rate(Cell *);
double &cell_migration_speed(Cell *);
double &cell_migration_bias(Cell *);
double &cell_migration_persistence_time(Cell *);
double &cell_chemotaxis_sensitivity(Cell *, int);
double &cell_cell_adhesion_strength(Cell *);
double &cell_cell_adhesion_elastic_constant(Cell *);
double &cell_adhesion_affinity_to_type(Cell *, int);
double &cell_relative_maximum_adhesion_distance(Cell *);
double &cell_cell_repulsion(Cell *);
double &cell_basement_membrane_adhesion(Cell *);
double &cell_basement_membrane_repulsion(Cell *);
double &cell_phagocytose_apoptotic(Cell *);
double &cell_phagocytose_necrotic(Cell *);
double &cell_phagocytose_other_dead(Cell *);
double &cell_phagocytose_live_cell_type(Cell *, int);
double &cell_attack_type(Cell *, int);
double &cell_fuse_to_type(Cell *, int);
double &cell_transform_to_type(Cell *, int);
double &cell_asymmetric_division_to_type(Cell *, int);
double &cell_is_movable(Cell *);
double &cell_immunogenicity_to_type(Cell *, int);
double &cell_attachment_rate(Cell *);
double &cell_detachment_rate(Cell *);
double &maximum_number_attachments(Cell *);
double &cell_attack_damage_rate(Cell *);
double &cell_attack_duration(Cell *);
double &cell_damage_rate(Cell *);
double &cell_custom_behavior(Cell *, int);

double &cell_secretion_rate_base(Cell_Definition *, int);
double &cell_secretion_target_base(Cell_Definition *, int);
double &cell_uptake_rate_base(Cell_Definition *, int);
double &cell_net_export_rate_base(Cell_Definition *, int);
double &cell_cycle_entry_rate_base(Cell_Definition *);
double &phase_exit_rate_base(Cell_Definition *, int);
double &cell_apoptosis_rate_base(Cell_Definition *);
double &cell_necrosis_rate_base(Cell_Definition *);
double &cell_migration_speed_base(Cell_Definition *);
double &cell_migration_bias_base(Cell_Definition *);
double &cell_migration_persistence_time_base(Cell_Definition *);
double &cell_chemotaxis_sensitivity_base(Cell_Definition *, int);
double &cell_cell_adhesion_strength_base(Cell_Definition *);
double &cell_cell_adhesion_elastic_constant_base(Cell_Definition *);
double &cell_adhesion_affinity_to_type_base(Cell_Definition *, int);
double &cell_relative_maximum_adhesion_distance_base(Cell_Definition *);
double &cell_cell_repulsion_base(Cell_Definition *);
double &cell_basement_membrane_adhesion_base(Cell_Definition *);
double &cell_basement_membrane_repulsion_base(Cell_Definition *);
double &cell_phagocytose_apoptotic_base(Cell_Definition *);
double &cell_phagocytose_necrotic_base(Cell_Definition *);
double &cell_phagocytose_other_dead_base(Cell_Definition *);
double &cell_phagocytose_live_cell_type_base(Cell_Definition *, int);
double &cell_attack_type_base(Cell_Definition *, int);
double &cell_fuse_to_type_base(Cell_Definition *, int);
double &cell_transform_to_type_base(Cell_Definition *, int);
double &cell_asymmetric_division_to_type_base(Cell_Definition *, int);
double &cell_is_movable_base(Cell_Definition *);
double &cell_immunogenicity_to_type_base(Cell_Definition *, int);
double &cell_attachment_rate_base(Cell_Definition *);
double &cell_detachment_rate_base(Cell_Definition *);
double &maximum_number_attachments_base(Cell_Definition *);
double &cell_attack_damage_rate_base(Cell_Definition *);
double &cell_attack_duration_base(Cell_Definition *);
double &cell_damage_rate_base(Cell_Definition *);
double &cell_custom_behavior_base(Cell_Definition *, int);

void add_signal(const std::string &name, SignalValue value);
void add_signal(const std::vector<std::string> &synonyms, SignalValue value);

void add_behavior(const std::string &name, BehaviorValue value, BehaviorBaseValue base_value);
void add_behavior(const std::vector<std::string> &synonyms, BehaviorValue value, BehaviorBaseValue base_value);

// display dictionaries 
void display_signal_dictionary( void );
void display_behavior_dictionary( void );

void display_signal_dictionary( std::ostream& os );
void display_behavior_dictionary( std::ostream& os );

void display_signal_dictionary_with_synonyms( void );
void display_behavior_dictionary_with_synonyms( void );
void display_signal_dictionary_with_synonyms( std::ostream& os );
void display_behavior_dictionary_with_synonyms( std::ostream& os );



/* signal functions */
bool signal_exists(const std::string &signal_name);

// find index for named signal (returns -1 if not found)
// int find_signal_index( std::string signal_name );

// coming soon: 
// std::vector<int> find_signal_indices( std::vector<std::string> signal_names );

// // get the name of a signal index 
// std::string signal_name( int i );

// create a full signal vector 
std::unordered_map<std::string, double> get_signals( Cell* pCell );

// create a subset of the signal vector with the supplied indices
std::vector<double> get_selected_signals(Cell *pCell, const std::vector<std::string> &names);

// grab a single signal by its index or name
double get_single_signal( Cell* pCell, const std::string &name );

/* behavior functions */ 
bool behavior_exists( const std::string &behavior_name );
// find index for named behavior / response / parameter (returns -1 if not found)
// int find_parameter_index( std::string response_name );
// int find_behavior_index( std::string response_name );

// std::vector<int> find_behavior_indices( std::vector<std::string> behavior_names );

// get the name of a behavior index 
std::string behavior_name( int i );

// write a full behavior vector (phenotype parameters) to the cell
void set_behaviors(Cell *pCell, const std::vector<double> &parameters);

// write a selected set of behavior parameters to the cell
void set_selected_behaviors(Cell *pCell, const std::vector<std::string> &names, const std::vector<double> &parameters);

// write a single behavior parameter
void set_single_behavior(Cell *pCell, const std::string &name, const double &parameter);

/* get current behaviors */ 

// get all current behavior
std::unordered_map<std::string, double> get_behaviors( Cell* pCell );

// get selected current behavior
std::vector<double> get_behaviors( Cell* pCell , const std::vector<std::string> &names );

// get single current behavior 
double get_single_behavior( Cell* pCell , const std::string &name );

/* get base behaviors (from cell definition) */ 

// get all base behaviors (from cell's definition) 
std::unordered_map<std::string, double> get_base_behaviors( Cell* pCell );

// get selected base behaviors (from cell's definition)
std::unordered_map<std::string, double> get_base_behaviors( Cell* pCell , const std::vector<std::string> &names );

// get single base behavior (from cell's definition)
double get_single_base_behavior( Cell* pCell , const std::string &name ); 

double get_single_base_behavior( Cell_Definition* pCD , const std::string &name ); 


}; 

#endif 
