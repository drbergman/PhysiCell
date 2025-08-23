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
 
#include "./PhysiCell_signal_behavior.h"

using namespace BioFVM; 

namespace PhysiCell{

// std::vector<double> signal_scales; 
typedef void (Secretion::*Advancer)(Basic_Agent *pCell, Phenotype &phenotype, double dt);

std::unordered_map<std::string, std::shared_ptr<Signal>> all_signals;
std::unordered_map<std::string, std::shared_ptr<Behavior>> all_behaviors;

double substrate_density(Cell *pCell, int substrate_index) { return pCell->nearest_density_vector()[substrate_index]; }
double internalized_substrate_density(Cell *pCell, int substrate_index) { return pCell->phenotype.molecular.internalized_total_substrates[substrate_index] / pCell->phenotype.volume.total; }
double substrate_gradient_norm(Cell *pCell, int substrate_index) { return norm(pCell->nearest_gradient(substrate_index)); }
double cell_pressure(Cell *pCell) { return pCell->state.simple_pressure; }
double cell_volume(Cell *pCell) { return pCell->phenotype.volume.total; }

double cell_contact(Cell *pCell, int cell_type)
{
	double signal_value = 0.0;
	for (auto neighbor : pCell->state.neighbors)
	{
		if (neighbor->type == cell_type)
		{ signal_value += 1.0; }
	}
	return signal_value;
}

double live_cell_contact(Cell *pCell)
{
	double signal_value = 0.0;
	for (auto neighbor : pCell->state.neighbors)
	{
		if (!neighbor->phenotype.death.dead)
		{ signal_value += 1.0; }
	}
	return signal_value;
}

double dead_cell_contact(Cell *pCell)
{
	double signal_value = 0.0;
	for (auto neighbor : pCell->state.neighbors)
	{
		if (neighbor->phenotype.death.dead)
		{ signal_value += 1.0; }
	}
	return signal_value;
}

double apoptotic_cell_contact(Cell *pCell)
{
	double signal_value = 0.0;
	for (auto neighbor : pCell->state.neighbors)
	{
		if (neighbor->phenotype.cycle.current_phase().code == PhysiCell_constants::apoptotic)
		{ signal_value += 1.0; }
	}
	return signal_value;
}

double necrotic_cell_contact(Cell *pCell)
{
	double signal_value = 0.0;
	for (auto neighbor : pCell->state.neighbors)
	{
		if (neighbor->is_necrotic())
		{ signal_value += 1.0; }
	}
	return signal_value;
}

double other_dead_cell_contact(Cell *pCell)
{
	double signal_value = 0.0;
	for (auto neighbor : pCell->state.neighbors)
	{
		if (neighbor->phenotype.death.dead && (neighbor->phenotype.cycle.current_phase().code != PhysiCell_constants::apoptotic) && !neighbor->is_necrotic())
		{ signal_value += 1.0; }
	}
	return signal_value;
}

double contact_with_basement_membrane(Cell *pCell) { return (double)(pCell->state.contact_with_basement_membrane); }
double cell_damage(Cell *pCell) { return pCell->phenotype.cell_integrity.damage; };
double damage_delivered(Cell *pCell) { return pCell->phenotype.cell_interactions.total_damage_delivered; };
double is_attacking(Cell *pCell) { return pCell->phenotype.cell_interactions.pAttackTarget ? 1.0 : 0.0; }
double is_dead(Cell *pCell) { return (double)pCell->phenotype.death.dead ? 0.0 : 1.0; }
double total_attack_time(Cell *pCell) { return pCell->state.total_attack_time; }
double get_current_time(Cell *) {return PhysiCell_globals.current_time; }
double cell_custom_signal(Cell *pCell, int i) { return pCell->custom_data.variables[i].value; };
double is_apoptotic(Cell *pCell) { return (double)(pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::apoptotic); }
double is_necrotic(Cell *pCell) { return (double)(pCell->is_necrotic()); }

double &cell_secretion_rate(Cell *pCell, int i) { return pCell->phenotype.secretion.secretion_rates[i]; }
double &cell_secretion_target(Cell *pCell, int i) { return pCell->phenotype.secretion.saturation_densities[i]; }
double &cell_uptake_rate(Cell *pCell, int i) { return pCell->phenotype.secretion.uptake_rates[i]; }
double &cell_net_export_rate(Cell *pCell, int i) { return pCell->phenotype.secretion.net_export_rates[i]; }
double &cell_cycle_entry_rate(Cell *pCell) { return pCell->phenotype.cycle.data.exit_rate(0); }

double &phase_exit_rate(Cell *pCell, int i)
{
	auto &phases = pCell->phenotype.cycle.model().phases;
	if (i >= phases.size())
	{
		std::cerr << "ERROR: Attempting to access cycle phase " << i << " exit rate..." << std::endl
				  << "...but cells of type " << pCell->type_name << " only have 0-"
				  << phases.size() - 1 << " phases." << std::endl;
		exit(-1);
	}
	return pCell->phenotype.cycle.data.exit_rate(i);
}

double &cell_apoptosis_rate(Cell *pCell)
{
	static int apoptosis_index = pCell->phenotype.death.find_death_model_index(PhysiCell_constants::apoptosis_death_model);
	return pCell->phenotype.death.rates[apoptosis_index];
}

double &cell_necrosis_rate(Cell *pCell)
{
	static int necrosis_index = pCell->phenotype.death.find_death_model_index(PhysiCell_constants::necrosis_death_model);
	return pCell->phenotype.death.rates[necrosis_index];
}

double &cell_migration_speed(Cell *pCell) { return pCell->phenotype.motility.migration_speed; }
double &cell_migration_bias(Cell *pCell) { return pCell->phenotype.motility.migration_bias; }
double &cell_migration_persistence_time(Cell *pCell) { return pCell->phenotype.motility.persistence_time; }
double &cell_chemotaxis_sensitivity(Cell *pCell, int i) { return pCell->phenotype.motility.chemotactic_sensitivities[i]; }
double &cell_cell_adhesion_strength(Cell *pCell) { return pCell->phenotype.mechanics.cell_cell_adhesion_strength; }
double &cell_cell_adhesion_elastic_constant(Cell *pCell) { return pCell->phenotype.mechanics.attachment_elastic_constant; }
double &cell_adhesion_affinity_to_type(Cell *pCell, int i) { return pCell->phenotype.mechanics.cell_adhesion_affinities[i]; }
double &cell_relative_maximum_adhesion_distance(Cell *pCell) { return pCell->phenotype.mechanics.relative_maximum_adhesion_distance; }
double &cell_cell_repulsion(Cell *pCell) { return pCell->phenotype.mechanics.cell_cell_repulsion_strength; }
double &cell_basement_membrane_adhesion(Cell *pCell) { return pCell->phenotype.mechanics.cell_BM_adhesion_strength; }
double &cell_basement_membrane_repulsion(Cell *pCell) { return pCell->phenotype.mechanics.cell_BM_repulsion_strength; }
double &cell_phagocytose_apoptotic(Cell *pCell) { return pCell->phenotype.cell_interactions.apoptotic_phagocytosis_rate; }
double &cell_phagocytose_necrotic(Cell *pCell) { return pCell->phenotype.cell_interactions.necrotic_phagocytosis_rate; }
double &cell_phagocytose_other_dead(Cell *pCell) { return pCell->phenotype.cell_interactions.other_dead_phagocytosis_rate; }
double &cell_phagocytose_live_cell_type(Cell *pCell, int i) { return pCell->phenotype.cell_interactions.live_phagocytosis_rates[i]; }
double &cell_attack_type(Cell *pCell, int i) { return pCell->phenotype.cell_interactions.attack_rates[i]; }
double &cell_fuse_to_type(Cell *pCell, int i) { return pCell->phenotype.cell_interactions.fusion_rates[i]; }
double &cell_transform_to_type(Cell *pCell, int i) { return pCell->phenotype.cell_transformations.transformation_rates[i]; }
double &cell_asymmetric_division_to_type(Cell *pCell, int i) { return pCell->phenotype.cycle.asymmetric_division.asymmetric_division_probabilities[i]; }
double &cell_is_movable(Cell *pCell) { return (double &)pCell->is_movable; }
double &cell_immunogenicity_to_type(Cell *pCell, int i) { return pCell->phenotype.cell_interactions.immunogenicities[i]; };
double &cell_attachment_rate(Cell *pCell) { return pCell->phenotype.mechanics.attachment_rate; }
double &cell_detachment_rate(Cell *pCell) { return pCell->phenotype.mechanics.detachment_rate; }
double &maximum_number_attachments(Cell *pCell) { return (double &)pCell->phenotype.mechanics.maximum_number_of_attachments; }
double &cell_attack_damage_rate(Cell *pCell) { return pCell->phenotype.cell_interactions.attack_damage_rate; }
double &cell_attack_duration(Cell *pCell) { return pCell->phenotype.cell_interactions.attack_duration; }
double &cell_damage_rate(Cell *pCell) { return pCell->phenotype.cell_integrity.damage_rate; }
double &cell_damage_repair_rate(Cell *pCell) { return pCell->phenotype.cell_integrity.damage_repair_rate; }
double &cell_custom_behavior(Cell *pCell, int i) { return pCell->custom_data.variables[i].value; };


void setup_signal_behavior_dictionaries( void )
{
	extern std::unordered_map<std::string,int> cell_definition_indices_by_name; 
	extern std::unordered_map<int,int> cell_definition_indices_by_type; 
	extern std::unordered_map<int,Cell_Definition*> cell_definitions_by_type; 
	
	create_base_cells();
	
	// set key parameters on number of signals, etc. 
	// make registry of signals 
	// make registry of responses 
	
	static bool setup_done = false; 
	if( setup_done == true )
	{ return; }
	setup_done = true; 
	
	int m = microenvironment.number_of_densities(); 
	int n = cell_definition_indices_by_name.size(); 
	
	all_signals.clear(); 
	all_behaviors.clear();
	
	// construct signals 
	
	std::vector<std::string> signal_synonyms;
	SignalValue signal_value;
	BehaviorValue behavior_value;
	
	// substrate densities 
	for( int i=0; i < m ; i++ )
	{
		signal_value = [i](Cell *pCell) -> double
		{ return substrate_density(pCell, i); };
		add_signal(microenvironment.density_names[i], signal_value);
	}

    // internalized substrates 
    for( int i=0; i < m ; i++ )
	{
		signal_synonyms = {"intracellular " + microenvironment.density_names[i], "internalized " + microenvironment.density_names[i]};
		signal_value = [i](Cell *pCell) -> double
		{ return internalized_substrate_density(pCell, i); };
		add_signal(signal_synonyms, signal_value);
	}

    // substrate gradients 
	for( int i=0; i < m ; i++ )
	{
		signal_synonyms = {microenvironment.density_names[i] + " gradient",
						   "grad(" + microenvironment.density_names[i] + ")",
						   "gradient of " + microenvironment.density_names[i]};
		signal_value = [i](Cell *pCell) -> double
		{ return substrate_gradient_norm(pCell, i); };
		add_signal(signal_synonyms, signal_value);
	}

	// mechanical pressure 
	add_signal("pressure", cell_pressure);

	// total volume
	add_signal("volume", cell_volume);

	// contact with each cell type
	for (int i = 0; i < n; i++)
	{
		Cell_Definition* pCD = cell_definitions_by_type[i]; 
		signal_synonyms =  {"contact with " + pCD->name};
		signal_value = [i](Cell *pCell) -> double
		{ return cell_contact(pCell, i); };;
		add_signal(signal_synonyms, signal_value);
	}
	
	// contact with (any) live cell 
	signal_synonyms = {"contact with live cell", "contact with live cells"};
	add_signal(signal_synonyms, live_cell_contact);

	// contact with (any) dead cell 
	signal_synonyms = {"contact with dead cell", "contact with dead cells"};
	add_signal(signal_synonyms, dead_cell_contact);

	// contact with apoptotic cell 
	signal_synonyms = {"contact with apoptotic cell", "contact with apoptotic cells"};
	add_signal(signal_synonyms, apoptotic_cell_contact);

	// contact with necrotic cell 
	signal_synonyms = {"contact with necrotic cell", "contact with necrotic cells"};
	add_signal(signal_synonyms, necrotic_cell_contact);

	// contact with other dead cell 
	signal_synonyms = {"contact with other dead cell", "contact with other dead cells"};
	add_signal(signal_synonyms, other_dead_cell_contact);

	// contact with basement membrane 
	signal_synonyms = {"contact with basement membrane", "contact with BM"};
	add_signal(signal_synonyms, contact_with_basement_membrane);

	// damage state 
	add_signal("damage", cell_damage);

	// damage delivered
	signal_synonyms = {"damage delivered", "total damage delivered"};
	add_signal(signal_synonyms, damage_delivered);

	// attacking yes/no?  
	signal_synonyms = {"attacking", "is attacking"};
	add_signal(signal_synonyms, is_attacking);

	// live / dead state 
	signal_synonyms = {"dead", "is dead"};
	add_signal(signal_synonyms, is_dead);

	// total attack time 
	add_signal("total_attack_time", total_attack_time);

	// current time
	signal_synonyms = {"time", "current time", "global time"};
	add_signal(signal_synonyms, get_current_time);

	// custom signals
	std::vector<std::string> base_custom_names = {"custom:", "custom: ", "custom "};
	for( int nc=0 ; nc < cell_defaults.custom_data.variables.size() ; nc++ )
	{
		signal_synonyms.resize(base_custom_names.size());
		for (size_t i = 0; i < base_custom_names.size(); i++)
		{
			signal_synonyms[i] = base_custom_names[i] + cell_defaults.custom_data.variables[nc].name;
		}
		signal_value = [nc](Cell *pCell) -> double
		{ return cell_custom_signal(pCell, nc); };
		add_signal(signal_synonyms, signal_value);

		behavior_value = [nc](Cell *pCell) -> double&
		{ return cell_custom_behavior(pCell, nc); };
		add_behavior(signal_synonyms, behavior_value);
	}

	// is apoptotic
	signal_synonyms = {"apoptotic", "is_apoptotic"};
	add_signal(signal_synonyms, is_apoptotic);

	signal_synonyms = {"necrotic", "is_necrotic"};
	add_signal(signal_synonyms, is_necrotic);

/*
	// immunogenicity to each cell type 
	for( int i=0; i < n ; i++ )
	{
		map_index++; 
		Cell_Definition* pCD = cell_definitions_by_type[i]; 
		std::string temp =  "immunogenicity to " + pCD->name; 
		signal_to_int[temp] = map_index; 
		int_to_signal[map_index] = temp; 		
				// synonyms 
		std::string temp1 = "immunogenicity to cell type " + std::to_string( pCD->type ); 
		signal_to_int[temp1] = map_index; 
	}
*/



	/* add new signals above this line */

	// construct behaviors 

	std::vector<std::string> behavior_synonyms;

	for( int i=0; i < m ; i++ )
	{
		std::string name = microenvironment.density_names[i];

		// secretion rate 
		behavior_value = [i](Cell *pCell) -> double&
		{ return cell_secretion_rate(pCell, i); };
		add_behavior(name + " secretion", behavior_value);

		// secretion target 
		behavior_synonyms = {name + " secretion target", name + " secretion saturation density"};
		behavior_value = [i](Cell *pCell) -> double&
		{ return cell_secretion_target(pCell, i); };
		add_behavior(behavior_synonyms, behavior_value);

		// uptake rate 
		behavior_value = [i](Cell *pCell) -> double&
		{ return cell_uptake_rate(pCell, i); };
		add_behavior(name + " uptake", behavior_value);

		// net export rate 
		behavior_value = [i](Cell *pCell) -> double&
		{ return cell_net_export_rate(pCell, i); };
		add_behavior(name + " export", behavior_value);
	}

	behavior_synonyms = {"cycle entry", "exit from cycle phase 0"};
	add_behavior(behavior_synonyms, cell_cycle_entry_rate);

	// other cyle phases 
	for( int i=1; i < 6; i++ )
	{
		behavior_value = [i](Cell *pCell) -> double&
		{ return phase_exit_rate(pCell, i); };
		add_behavior("exit from cycle phase " + std::to_string(i), behavior_value);
	}

	// apoptosis
	add_behavior("apoptosis", cell_apoptosis_rate);

	// necrosis
	add_behavior("necrosis", cell_necrosis_rate);

	// migration speed
	add_behavior("migration speed", cell_migration_speed);

	// migration bias
	add_behavior("migration bias", cell_migration_bias);

	// migration persistence time
	add_behavior("migration persistence time", cell_migration_persistence_time);

	// chemotactic sensitivities 
	for( int i=0; i < m ; i++ )
	{
		behavior_synonyms = {"chemotactic response to " + microenvironment.density_names[i],
		                     "chemotactic sensitivity to " + microenvironment.density_names[i]};
		behavior_value = [i](Cell *pCell) -> double&
		{ return cell_chemotaxis_sensitivity(pCell, i); };
		add_behavior(behavior_synonyms, behavior_value);
	}
	
	// cell-cell adhesion 
	add_behavior("cell-cell adhesion", cell_cell_adhesion_strength);

	// cell-cell adhesion elastic constant
	add_behavior("cell-cell adhesion elastic constant", cell_cell_adhesion_elastic_constant);

    // cell adhesion affinities 
	// cell-type specific adhesion 
	for( int i=0; i < n ; i++ )
	{
		behavior_value = [i](Cell *pCell) -> double&
		{ return cell_adhesion_affinity_to_type(pCell, i); };
		add_behavior("adhesive affinity to " + cell_definitions_by_type[i]->name,
					 behavior_value);
	}

	// max adhesion distance
	add_behavior("relative maximum adhesion distance", cell_relative_maximum_adhesion_distance);

	// cell-cell repulsion
	add_behavior("cell-cell repulsion", cell_cell_repulsion);

	// cell-BM adhesion
	behavior_synonyms = {"cell-BM adhesion", "cell-membrane adhesion"};
	add_behavior(behavior_synonyms, cell_basement_membrane_adhesion);

	// cell-BM repulsion 
	behavior_synonyms = {"cell-BM repulsion", "cell-membrane repulsion"};
	add_behavior(behavior_synonyms, cell_basement_membrane_repulsion);

	// phagocytosis of apoptotic cell
	behavior_synonyms = {"phagocytose apoptotic cell",
						 "phagocytosis of apoptotic cell",
						 "phagocytosis of apoptotic cells"};
	add_behavior(behavior_synonyms, cell_phagocytose_apoptotic);

	// phagocytosis of necrotic cell
	behavior_synonyms = {"phagocytose necrotic cell",
						 "phagocytosis of necrotic cell",
						 "phagocytosis of necrotic cells"};
	add_behavior(behavior_synonyms, cell_phagocytose_necrotic);

	// phagocytosis of other dead cell
	behavior_synonyms = {"phagocytose other dead cell",
						 "phagocytosis of other dead cell",
						 "phagocytosis of other dead cells"};
	add_behavior(behavior_synonyms, cell_phagocytose_other_dead);

	// phagocytosis of each live cell type 
	for( int i=0; i < n ; i++ )
	{
		behavior_synonyms = {"phagocytose " + cell_definitions_by_type[i]->name,
							 "phagocytosis of " + cell_definitions_by_type[i]->name};
		behavior_value = [i](Cell *pCell) -> double&
		{ return cell_phagocytose_live_cell_type(pCell, i); };
		add_behavior(behavior_synonyms, behavior_value);
	}

	// attack of each live cell type 
	for( int i=0; i < n ; i++ )
	{
		behavior_value = [i](Cell *pCell) -> double&
		{ return cell_attack_type(pCell, i); };
		add_behavior("attack " + cell_definitions_by_type[i]->name, behavior_value);
	}

	// fusion 
	for( int i=0; i < n ; i++ )
	{
		behavior_value = [i](Cell *pCell) -> double&
		{ return cell_fuse_to_type(pCell, i); };
		behavior_synonyms = {"fuse to " + cell_definitions_by_type[i]->name};
		add_behavior(behavior_synonyms, behavior_value);
	}

	// transformation / transition 
	for( int i=0; i < n ; i++ )
	{
		behavior_synonyms = {"transition to " + cell_definitions_by_type[i]->name,
							 "transform to " + cell_definitions_by_type[i]->name};
		behavior_value = [i](Cell *pCell) -> double&
		{ return cell_transform_to_type(pCell, i); };
		add_behavior(behavior_synonyms, behavior_value);
	}

	// asymmetic division
	for( int i=0; i < n ; i++ )
	{
		behavior_value = [i](Cell *pCell) -> double&
		{ return cell_asymmetric_division_to_type(pCell, i); };
		add_behavior("asymmetric division to " + cell_definitions_by_type[i]->name, behavior_value);
	}

	// is movable
	behavior_synonyms = {"is_movable", "is movable", "movable"};
	add_behavior(behavior_synonyms, cell_is_movable);

	// immunogenicity to each cell type
	for (int i = 0; i < n; i++)
	{
		behavior_value = [i](Cell *pCell) -> double&
		{ return cell_immunogenicity_to_type(pCell, i); };
		add_behavior("immunogenicity to " + cell_definitions_by_type[i]->name, behavior_value);
	}

	// cell attachment rate
	add_behavior("cell attachment rate", cell_attachment_rate);

	// cell detachment rate
	add_behavior("cell detachment rate", cell_detachment_rate);

	// maximum number of cell attachments
	add_behavior("maximum number of cell attachments", maximum_number_attachments);

	// attack damage rate
	add_behavior("attack damage rate", cell_attack_damage_rate);

	// attack duration
	add_behavior("attack duration", cell_attack_duration);

	// damage rate
	add_behavior("damage rate", cell_damage_rate);

	// damage repair rate
	add_behavior("damage repair rate", cell_damage_repair_rate);

	/* add new behaviors above this line */

    // resize scales; 
    // signal_scales.resize( int_to_signal.size() , 1.0 ); 

    display_signal_dictionary(); 
    display_behavior_dictionary(); 
/*
	// now create empty SR models for each cell definition 

	for( int i=0 ; i < cell_definitions_by_index.size() ; i++ )
	{ create_SR_model( *cell_definitions_by_index[i] ); }
*/    

	return;
}

void add_signal( std::string name, SignalValue value )
{ return add_signal(std::vector<std::string>{std::move(name)}, std::move(value)); }

void add_signal(const std::vector<std::string> &synonyms, SignalValue value)
{
	// allocate the callable once
	auto signal_ptr = std::make_shared<Signal>(synonyms[0], std::move(value));
	
	// assign the same pointer to all synonyms
	for (const auto &name : synonyms)
	{ all_signals[name] = signal_ptr; }
	return;
}

std::unordered_map<std::string,Cell*> base_cells_by_name;

void create_base_cells()
{
	base_cells_by_name.clear();
	for( auto& cell_definition : cell_definitions_by_name )
	{
		Cell_Definition* pCD = cell_definition.second;
		Cell* pCell = new Cell;
		pCell->set_total_volume( pCell->phenotype.volume.total );

		pCell->type = pCD->type;
		pCell->type_name = pCD->name;

		pCell->custom_data = pCD->custom_data;
		pCell->parameters = pCD->parameters;
		pCell->functions = pCD->functions;

		pCell->phenotype = pCD->phenotype;
		if (pCell->phenotype.intracellular)
			pCell->phenotype.intracellular->start();

		pCell->is_movable = false; //  true;
		pCell->is_out_of_domain = true;
		pCell->displacement.resize(3, 0.0); // state?

		pCell->set_total_volume(pCell->phenotype.volume.total);

		// store the base cell by name
		base_cells_by_name[pCD->name] = std::move(pCell);
	}
}

void add_behavior( std::string name, BehaviorValue value )
{ return add_behavior({std::vector<std::string>{std::move(name)}}, std::move(value)); }

void add_behavior(const std::vector<std::string> &synonyms, BehaviorValue value)
{
	auto behavior_ptr = std::make_shared<Behavior>(synonyms[0], std::move(value));
	for (auto name : synonyms)
	{ all_behaviors[name] = behavior_ptr; }
	return;
}

// double& signal_scale( std::string signal_name )
// { return signal_scales[ find_signal_index(signal_name) ]; }

// double& signal_scale( int signal_index  )
// { return signal_scales[signal_index]; }

void display_signal_dictionary( std::ostream& os )
{
	os << "Signals: " << std::endl 
	   << "=======" << std::endl; 
	std::vector<std::string> displayed_signals = {};
	std::string signal_name;
	for (auto signal : all_signals)
	{
		signal_name = signal.second.get()->get_name();
		for (auto displayed_signal : displayed_signals)
		{
			if (displayed_signal == signal_name)
			{ continue; }
		}
		os << "- " << signal_name << std::endl;
		displayed_signals.push_back(signal_name);
	}
	os << std::endl;
    return; 
}

void display_signal_dictionary( void )
{ display_signal_dictionary( std::cout); std::cout << std::endl; }


void display_signal_dictionary_with_synonyms( void )
{ display_signal_dictionary_with_synonyms( std::cout ); }
/*
	std::cout << "Signals (with synonyms): " << std::endl 
			  << "=======================" << std::endl; 
	for( auto it = signal_to_int.begin() ; it != signal_to_int.end() ; it++ )
	{ std::cout << it->second << " : " << it->first << std::endl; }
	std::cout << std::endl << std::endl;  	
    return; 
*/

void display_signal_dictionary_with_synonyms( std::ostream& os )
{
	os << "Signals (with synonyms): " << std::endl 
	   << "=======================" << std::endl;
	for (auto signal : all_signals)
	{
		os << "- " << signal.first << std::endl;
	}
	os << std::endl << std::endl;  	
    return; 
}

void display_behavior_dictionary( std::ostream& os )
{
	os << "Behaviors: " << std::endl 
	   << "=========" << std::endl; 
	std::vector<std::string> displayed_behaviors = {};
	std::string behavior_name;
	for (auto behavior : all_behaviors)
	{
		behavior_name = behavior.second.get()->get_name();
		for (auto displayed_behavior : displayed_behaviors)
		{
			if (displayed_behavior == behavior_name)
			{
				continue;
			}
		}
		os << "- " << behavior_name << std::endl;
		displayed_behaviors.push_back(behavior_name);
	}
	os << std::endl;
    return; 
}

void display_behavior_dictionary( void )
{
	display_behavior_dictionary( std::cout );
	std::cout << std::endl; 
	return; 
}

void display_behavior_dictionary_with_synonyms( std::ostream& os )
{
	os << "Behaviors (with synonyms): " << std::endl
	   << "=========================" << std::endl;
	for (auto behavior : all_behaviors)
	{
		os << "- " << behavior.first << std::endl;
	}
	os << std::endl << std::endl;
	return;
}

void display_behavior_dictionary_with_synonyms( void )
{ display_behavior_dictionary_with_synonyms( std::cout ); return; }
/*
	std::cout << "Behaviors (with synonyms): " << std::endl 
			  << "=========================" << std::endl; 
	for( auto it = behavior_to_int.begin() ; it != behavior_to_int.end() ; it++ )
	{ std::cout << it->second << " : " << it->first << std::endl; }
	std::cout << std::endl << std::endl;  	
    return; 
*/	

bool signal_exists( std::string signal_name )
{ return all_signals.find(signal_name) != all_signals.end(); }

bool behavior_exists( std::string behavior_name )
{ return all_behaviors.find(behavior_name) != all_behaviors.end(); }

// int find_signal_index( std::string signal_name )
// {
// 	auto search = signal_to_int.find( signal_name );
// 	// safety first! 
// 	if( search != signal_to_int.end() )
//     { return search->second; }   

// 	std::cout << "having trouble finding " << signal_name << std::endl; 

//     return -1; 
// }

// std::vector<int> find_signal_indices( std::vector<std::string> signal_names )
// {
// 	std::vector<int> output( signal_names.size() , 0 ); 
// 	for( int n=0; n < signal_names.size(); n++ )
// 	{ output[n] = find_signal_index(signal_names[n]); }
// 	return output; 
// }

// std::string signal_name( int i )
// {
// 	if( i >= 0 && i < int_to_signal.size() )
// 	{ return int_to_signal[i]; }	
// 	return "not found"; 
// }

// int find_parameter_index( std::string response_name )
// {
// 	auto search = behavior_to_int.find( response_name );
// 	if( search != behavior_to_int.end() )
//     { return search->second; }   
//     return -1; 
// }

// int find_behavior_index( std::string response_name )
// { return find_parameter_index(response_name); }

// std::vector<int> find_behavior_indices( std::vector<std::string> behavior_names )
// {
// 	std::vector<int> output( behavior_names.size() , 0 ); 
// 	for( int n=0; n < behavior_names.size(); n++ )
// 	{ output[n] = find_behavior_index(behavior_names[n]); }
// 	return output; 
// }

std::unordered_map<std::string, double> get_signals( Cell* pCell )
{
	std::unordered_map<std::string, double> output;
	std::vector<std::string> signals_got = {};
	std::string signal_name;
	for (auto& signal : all_signals)
	{
		signal_name = signal.second.get()->get_name();
		if (std::find(signals_got.begin(), signals_got.end(), signal_name) != signals_got.end())
		{ continue; }
		output[signal_name] = signal.second.get()->get_value(pCell);
		signals_got.push_back(signal_name);
	}
	return output;
}

// create a signal vector of only the cell contacts 
std::vector<double> get_cell_contact_signals( Cell* pCell )
{
	static int m = microenvironment.number_of_densities(); 
	static int n = cell_definition_indices_by_name.size(); 

	std::vector<double> output( n+2+3 , 0.0 ); 
	// process all neighbors 
	int dead_cells = 0; 
	int live_cells = 0; 
    int apop_cells = 0;
    int necro_cells = 0; 
    int other_dead_cells = 0; 

	for( int i=0; i < pCell->state.neighbors.size(); i++ )
	{
		Cell* pC = pCell->state.neighbors[i]; 
		if( pC->phenotype.death.dead == true )
		{
			dead_cells++; 
            if(pC->phenotype.cycle.current_phase().code == PhysiCell_constants::apoptotic )
            { apop_cells++; }

            if( pC->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_swelling || 
                pC->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_lysed || 
                pC->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic )
            { necro_cells++; }  	
		} 
		else
		{ live_cells++; } 
		int nCT = cell_definition_indices_by_type[pC->type]; 
		output[nCT] += 1; 
	}
    other_dead_cells = dead_cells - apop_cells - necro_cells; 

	output[n] = live_cells; 
	output[n+1] = dead_cells; 

	output[n+2] = apop_cells; 
	output[n+3] = necro_cells; 
	output[n+4] = other_dead_cells; 

	// rescale 
	// std::string search_for = "contact with " + cell_definitions_by_type[0]->name; 
	// static int scaling_start_index = find_signal_index( search_for ); 
	// for( int i=0; i < n+2 ; i++ )
	// { output[i] /= signal_scales[scaling_start_index+i]; }

	return output; 
}

std::vector<double> get_selected_signals( Cell* pCell , const std::vector<std::string> &signal_names )
{
	std::vector<double> signals( signal_names.size() , 0.0 );
	for( int i=0; i < signal_names.size(); i++ )
	{
		signals[i] = get_single_signal(pCell, signal_names[i]);
	}
	return signals; 
}

double get_single_signal( Cell* pCell, std::string name )
{
	return all_signals[name]->get_value(pCell);
}

// behaviors 
void set_behaviors( Cell* pCell , std::unordered_map<std::string, double> behavior_mappings )
{
	for (auto& behavior_mapping : behavior_mappings)
	{ set_single_behavior(pCell, behavior_mapping.first, behavior_mapping.second); }
	return; 
}

void set_single_behavior( Cell* pCell, std::string name , double parameter )
{ return all_behaviors[name]->set_value(pCell, parameter); }

std::unordered_map<std::string, double> get_behaviors( Cell* pCell )
{
	std::unordered_map<std::string, double> behavior_values;
	std::vector<std::string> behaviors_got = {};
	std::string behavior_name;
	for (const auto& behavior : all_behaviors)
	{
		behavior_name = behavior.second.get()->get_name();
		if (std::find(behaviors_got.begin(), behaviors_got.end(), behavior_name) != behaviors_got.end())
		{ continue; }
		behavior_values[behavior_name] = behavior.second->get_value(pCell);
		behaviors_got.push_back(behavior_name);
	}
	return behavior_values; 
}

double get_single_behavior( Cell* pCell , std::string name )
{ return all_behaviors[name]->get_value(pCell); }

std::vector<double> get_behaviors( Cell* pCell , std::vector<std::string> names )
{
	std::vector<double> parameters( names.size() , 0.0 ); 
	for( int n=0; n < names.size(); n++ )
	{ parameters[n] = get_single_behavior(pCell,names[n]); }
	return parameters; 
}

void set_selected_behaviors( Cell* pCell , std::vector<std::string> names , std::vector<double> parameters )
{
	for( int i=0 ; i < names.size() ; i++ )
	{ set_single_behavior(pCell, names[i], parameters[i]); }
	return; 
}

std::unordered_map<std::string, double> get_base_behaviors( Cell* pCell )
{
	std::unordered_map<std::string, double> base_behavior_map;
	std::vector<std::string> behaviors_got = {};
	for (const auto& behavior : all_behaviors)
	{
		std::string behavior_name = behavior.second->get_name();
		if (std::find(behaviors_got.begin(), behaviors_got.end(), behavior_name) != behaviors_got.end())
		{ continue; }
		base_behavior_map[behavior_name] = behavior.second->get_value(base_cells_by_name[pCell->type_name]);
		behaviors_got.push_back(behavior_name);
	}
	return base_behavior_map;
}

double get_single_base_behavior( Cell* pCell , std::string name )
{
	return all_behaviors[name]->get_value(base_cells_by_name[pCell->type_name]);
}

double get_single_base_behavior( Cell_Definition* pCD , std::string name )
{
	return all_behaviors[name]->get_value(base_cells_by_name[pCD->name]);
}

std::unordered_map<std::string, double> get_base_behaviors( Cell* pCell , std::vector<std::string> names )
{
	std::unordered_map<std::string, double> base_behavior_map;
	for( const auto& name : names )
	{ base_behavior_map[name] = get_single_base_behavior(pCell,name); }
	return base_behavior_map;
}

};
