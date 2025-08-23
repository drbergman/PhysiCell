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

typedef void (Secretion::*Advancer)(Basic_Agent *pCell, Phenotype &phenotype, double dt);

std::unordered_map<std::string, std::shared_ptr<Signal>> all_signals;
std::unordered_map<std::string, std::shared_ptr<Behavior>> all_behaviors;

double Behavior::get_base_value(Cell* pCell) const
{ return base_value(&get_cell_definition(pCell->type)); }

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
double is_dead(Cell *pCell) { return (double)pCell->phenotype.death.dead; }
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
double &cell_immunogenicity_to_type(Cell *pCell, int i) { return pCell->phenotype.cell_interactions.immunogenicities[i]; }
double &cell_attachment_rate(Cell *pCell) { return pCell->phenotype.mechanics.attachment_rate; }
double &cell_detachment_rate(Cell *pCell) { return pCell->phenotype.mechanics.detachment_rate; }
double &maximum_number_attachments(Cell *pCell) { return (double &)pCell->phenotype.mechanics.maximum_number_of_attachments; }
double &cell_attack_damage_rate(Cell *pCell) { return pCell->phenotype.cell_interactions.attack_damage_rate; }
double &cell_attack_duration(Cell *pCell) { return pCell->phenotype.cell_interactions.attack_duration; }
double &cell_damage_rate(Cell *pCell) { return pCell->phenotype.cell_integrity.damage_rate; }
double &cell_damage_repair_rate(Cell *pCell) { return pCell->phenotype.cell_integrity.damage_repair_rate; }
double &cell_custom_behavior(Cell *pCell, int i) { return pCell->custom_data.variables[i].value; };


double &cell_secretion_rate_base(Cell_Definition *pCD, int i) { return pCD->phenotype.secretion.secretion_rates[i]; }
double &cell_secretion_target_base(Cell_Definition *pCD, int i) { return pCD->phenotype.secretion.saturation_densities[i]; }
double &cell_uptake_rate_base(Cell_Definition *pCD, int i) { return pCD->phenotype.secretion.uptake_rates[i]; }
double &cell_net_export_rate_base(Cell_Definition *pCD, int i) { return pCD->phenotype.secretion.net_export_rates[i]; }
double &cell_cycle_entry_rate_base(Cell_Definition *pCD) { return pCD->phenotype.cycle.data.exit_rate(0); }

double &phase_exit_rate_base(Cell_Definition *pCD, int i)
{
	auto &phases = pCD->phenotype.cycle.model().phases;
	if (i >= phases.size())
	{
		std::cerr << "ERROR: Attempting to access cycle phase " << i << " exit rate..." << std::endl
				  << "...but cells of type " << pCD->name << " only have 0-"
				  << phases.size() - 1 << " phases." << std::endl;
		exit(-1);
	}
	return pCD->phenotype.cycle.data.exit_rate(i);
}

double &cell_apoptosis_rate_base(Cell_Definition *pCD)
{
	static int apoptosis_index = pCD->phenotype.death.find_death_model_index(PhysiCell_constants::apoptosis_death_model);
	return pCD->phenotype.death.rates[apoptosis_index];
}

double &cell_necrosis_rate_base(Cell_Definition *pCD)
{
	static int necrosis_index = pCD->phenotype.death.find_death_model_index(PhysiCell_constants::necrosis_death_model);
	return pCD->phenotype.death.rates[necrosis_index];
}

double &cell_migration_speed_base(Cell_Definition *pCD) { return pCD->phenotype.motility.migration_speed; }
double &cell_migration_bias_base(Cell_Definition *pCD) { return pCD->phenotype.motility.migration_bias; }
double &cell_migration_persistence_time_base(Cell_Definition *pCD) { return pCD->phenotype.motility.persistence_time; }
double &cell_chemotaxis_sensitivity_base(Cell_Definition *pCD, int i) { return pCD->phenotype.motility.chemotactic_sensitivities[i]; }
double &cell_cell_adhesion_strength_base(Cell_Definition *pCD) { return pCD->phenotype.mechanics.cell_cell_adhesion_strength; }
double &cell_cell_adhesion_elastic_constant_base(Cell_Definition *pCD) { return pCD->phenotype.mechanics.attachment_elastic_constant; }
double &cell_adhesion_affinity_to_type_base(Cell_Definition *pCD, int i) { return pCD->phenotype.mechanics.cell_adhesion_affinities[i]; }
double &cell_relative_maximum_adhesion_distance_base(Cell_Definition *pCD) { return pCD->phenotype.mechanics.relative_maximum_adhesion_distance; }
double &cell_cell_repulsion_base(Cell_Definition *pCD) { return pCD->phenotype.mechanics.cell_cell_repulsion_strength; }
double &cell_basement_membrane_adhesion_base(Cell_Definition *pCD) { return pCD->phenotype.mechanics.cell_BM_adhesion_strength; }
double &cell_basement_membrane_repulsion_base(Cell_Definition *pCD) { return pCD->phenotype.mechanics.cell_BM_repulsion_strength; }
double &cell_phagocytose_apoptotic_base(Cell_Definition *pCD) { return pCD->phenotype.cell_interactions.apoptotic_phagocytosis_rate; }
double &cell_phagocytose_necrotic_base(Cell_Definition *pCD) { return pCD->phenotype.cell_interactions.necrotic_phagocytosis_rate; }
double &cell_phagocytose_other_dead_base(Cell_Definition *pCD) { return pCD->phenotype.cell_interactions.other_dead_phagocytosis_rate; }
double &cell_phagocytose_live_cell_type_base(Cell_Definition *pCD, int i) { return pCD->phenotype.cell_interactions.live_phagocytosis_rates[i]; }
double &cell_attack_type_base(Cell_Definition *pCD, int i) { return pCD->phenotype.cell_interactions.attack_rates[i]; }
double &cell_fuse_to_type_base(Cell_Definition *pCD, int i) { return pCD->phenotype.cell_interactions.fusion_rates[i]; }
double &cell_transform_to_type_base(Cell_Definition *pCD, int i) { return pCD->phenotype.cell_transformations.transformation_rates[i]; }
double &cell_asymmetric_division_to_type_base(Cell_Definition *pCD, int i) { return pCD->phenotype.cycle.asymmetric_division.asymmetric_division_probabilities[i]; }
double &cell_is_movable_base(Cell_Definition *pCD) { return (double &)pCD->is_movable; }
double &cell_immunogenicity_to_type_base(Cell_Definition *pCD, int i) { return pCD->phenotype.cell_interactions.immunogenicities[i]; }
double &cell_attachment_rate_base(Cell_Definition *pCD) { return pCD->phenotype.mechanics.attachment_rate; }
double &cell_detachment_rate_base(Cell_Definition *pCD) { return pCD->phenotype.mechanics.detachment_rate; }
double &maximum_number_attachments_base(Cell_Definition *pCD) { return (double &)pCD->phenotype.mechanics.maximum_number_of_attachments; }
double &cell_attack_damage_rate_base(Cell_Definition *pCD) { return pCD->phenotype.cell_interactions.attack_damage_rate; }
double &cell_attack_duration_base(Cell_Definition *pCD) { return pCD->phenotype.cell_interactions.attack_duration; }
double &cell_damage_rate_base(Cell_Definition *pCD) { return pCD->phenotype.cell_integrity.damage_rate; }
double &cell_damage_repair_rate_base(Cell_Definition *pCD) { return pCD->phenotype.cell_integrity.damage_repair_rate; }
double &cell_custom_behavior_base(Cell_Definition *pCD, int i) { return pCD->custom_data.variables[i].value; };

void setup_signal_behavior_dictionaries( void )
{
	extern std::unordered_map<std::string,int> cell_definition_indices_by_name; 
	extern std::unordered_map<int,int> cell_definition_indices_by_type; 
	extern std::unordered_map<int,Cell_Definition*> cell_definitions_by_type; 
	
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
	
	std::vector<std::string> synonyms;
	SignalValue signal_value;
	BehaviorValue behavior_value;
	BehaviorBaseValue behavior_base_value;

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
		synonyms = {"intracellular " + microenvironment.density_names[i],
						   "internalized " + microenvironment.density_names[i]};
		signal_value = [i](Cell *pCell) -> double
		{ return internalized_substrate_density(pCell, i); };
		add_signal(synonyms, signal_value);
	}

    // substrate gradients 
	for( int i=0; i < m ; i++ )
	{
		synonyms = {microenvironment.density_names[i] + " gradient",
						   "grad(" + microenvironment.density_names[i] + ")",
						   "gradient of " + microenvironment.density_names[i]};
		signal_value = [i](Cell *pCell) -> double
		{ return substrate_gradient_norm(pCell, i); };
		add_signal(synonyms, signal_value);
	}

	// mechanical pressure 
	add_signal("pressure", cell_pressure);

	// total volume
	add_signal("volume", cell_volume);

	// contact with each cell type
	for (int i = 0; i < n; i++)
	{
		Cell_Definition* pCD = cell_definitions_by_type[i]; 
		signal_value = [i](Cell *pCell) -> double
		{ return cell_contact(pCell, i); };
		add_signal("contact with " + pCD->name, signal_value);
	}
	
	// contact with (any) live cell 
	synonyms = {"contact with live cell", "contact with live cells"};
	add_signal(synonyms, live_cell_contact);

	// contact with (any) dead cell 
	synonyms = {"contact with dead cell", "contact with dead cells"};
	add_signal(synonyms, dead_cell_contact);

	// contact with apoptotic cell 
	synonyms = {"contact with apoptotic cell", "contact with apoptotic cells"};
	add_signal(synonyms, apoptotic_cell_contact);

	// contact with necrotic cell 
	synonyms = {"contact with necrotic cell", "contact with necrotic cells"};
	add_signal(synonyms, necrotic_cell_contact);

	// contact with other dead cell 
	synonyms = {"contact with other dead cell", "contact with other dead cells"};
	add_signal(synonyms, other_dead_cell_contact);

	// contact with basement membrane 
	synonyms = {"contact with basement membrane", "contact with BM"};
	add_signal(synonyms, contact_with_basement_membrane);

	// damage state 
	add_signal("damage", cell_damage);

	// damage delivered
	synonyms = {"damage delivered", "total damage delivered"};
	add_signal(synonyms, damage_delivered);

	// attacking yes/no?  
	synonyms = {"attacking", "is attacking"};
	add_signal(synonyms, is_attacking);

	// live / dead state 
	synonyms = {"dead", "is dead"};
	add_signal(synonyms, is_dead);

	// total attack time 
	add_signal("total_attack_time", total_attack_time);

	// current time
	synonyms = {"time", "current time", "global time"};
	add_signal(synonyms, get_current_time);

	// custom signals
	std::vector<std::string> base_custom_names = {"custom:", "custom: ", "custom "};
	for( int nc=0 ; nc < cell_defaults.custom_data.variables.size() ; nc++ )
	{
		synonyms.resize(base_custom_names.size());
		for (size_t i = 0; i < base_custom_names.size(); i++)
		{
			synonyms[i] = base_custom_names[i] + cell_defaults.custom_data.variables[nc].name;
		}
		signal_value = [nc](Cell *pCell) -> double
		{ return cell_custom_signal(pCell, nc); };
		add_signal(synonyms, signal_value);

		behavior_value = [nc](Cell *pCell) -> double&
		{ return cell_custom_behavior(pCell, nc); };
		behavior_base_value = [nc](Cell_Definition *pCD) -> double
		{ return cell_custom_behavior_base(pCD, nc); };
		add_behavior(synonyms, behavior_value, behavior_base_value);
	}

	// is apoptotic
	synonyms = {"apoptotic", "is_apoptotic"};
	add_signal(synonyms, is_apoptotic);

	synonyms = {"necrotic", "is_necrotic"};
	add_signal(synonyms, is_necrotic);

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

	for( int i=0; i < m ; i++ )
	{
		std::string name = microenvironment.density_names[i];

		// secretion rate 
		behavior_value = [i](Cell *pCell) -> double&
		{ return cell_secretion_rate(pCell, i); };
		behavior_base_value = [i](Cell_Definition *pCD) -> double
		{ return cell_secretion_rate_base(pCD, i); };
		add_behavior(name + " secretion", behavior_value, behavior_base_value);

		// secretion target 
		synonyms = {name + " secretion target", name + " secretion saturation density"};
		behavior_value = [i](Cell *pCell) -> double&
		{ return cell_secretion_target(pCell, i); };
		behavior_base_value = [i](Cell_Definition *pCD) -> double
		{ return cell_secretion_target_base(pCD, i); };
		add_behavior(synonyms, behavior_value, behavior_base_value);

		// uptake rate 
		behavior_value = [i](Cell *pCell) -> double&
		{ return cell_uptake_rate(pCell, i); };
		behavior_base_value = [i](Cell_Definition *pCD) -> double
		{ return cell_uptake_rate_base(pCD, i); };
		add_behavior(name + " uptake", behavior_value, behavior_base_value);

		// net export rate 
		behavior_value = [i](Cell *pCell) -> double&
		{ return cell_net_export_rate(pCell, i); };
		behavior_base_value = [i](Cell_Definition *pCD) -> double
		{ return cell_net_export_rate_base(pCD, i); };
		add_behavior(name + " export", behavior_value, behavior_base_value);
	}

	synonyms = {"cycle entry", "exit from cycle phase 0"};
	add_behavior(synonyms, cell_cycle_entry_rate, cell_cycle_entry_rate_base);

	// other cyle phases 
	for( int i=1; i < 6; i++ )
	{
		behavior_value = [i](Cell *pCell) -> double&
		{ return phase_exit_rate(pCell, i); };
		behavior_base_value = [i](Cell_Definition *pCD) -> double
		{ return phase_exit_rate_base(pCD, i); };
		add_behavior("exit from cycle phase " + std::to_string(i), behavior_value, behavior_base_value);
	}

	// apoptosis
	add_behavior("apoptosis", cell_apoptosis_rate, cell_apoptosis_rate_base);

	// necrosis
	add_behavior("necrosis", cell_necrosis_rate, cell_necrosis_rate_base);

	// migration speed
	add_behavior("migration speed", cell_migration_speed, cell_migration_speed_base);

	// migration bias
	add_behavior("migration bias", cell_migration_bias, cell_migration_bias_base);

	// migration persistence time
	add_behavior("migration persistence time", cell_migration_persistence_time, cell_migration_persistence_time_base);

	// chemotactic sensitivities 
	for( int i=0; i < m ; i++ )
	{
		synonyms = {"chemotactic response to " + microenvironment.density_names[i],
					"chemotactic sensitivity to " + microenvironment.density_names[i]};
		behavior_value = [i](Cell *pCell) -> double&
		{ return cell_chemotaxis_sensitivity(pCell, i); };
		behavior_base_value = [i](Cell_Definition *pCD) -> double
		{ return cell_chemotaxis_sensitivity_base(pCD, i); };
		add_behavior(synonyms, behavior_value, behavior_base_value);
	}
	
	// cell-cell adhesion 
	add_behavior("cell-cell adhesion", cell_cell_adhesion_strength, cell_cell_adhesion_strength_base);

	// cell-cell adhesion elastic constant
	add_behavior("cell-cell adhesion elastic constant", cell_cell_adhesion_elastic_constant, cell_cell_adhesion_elastic_constant_base);

    // cell adhesion affinities 
	// cell-type specific adhesion 
	for( int i=0; i < n ; i++ )
	{
		behavior_value = [i](Cell *pCell) -> double&
		{ return cell_adhesion_affinity_to_type(pCell, i); };
		behavior_base_value = [i](Cell_Definition *pCD) -> double
		{ return cell_adhesion_affinity_to_type_base(pCD, i); };
		add_behavior("adhesive affinity to " + cell_definitions_by_type[i]->name,
					 behavior_value, behavior_base_value);
	}

	// max adhesion distance
	add_behavior("relative maximum adhesion distance", cell_relative_maximum_adhesion_distance, cell_relative_maximum_adhesion_distance_base);

	// cell-cell repulsion
	add_behavior("cell-cell repulsion", cell_cell_repulsion, cell_cell_repulsion_base);

	// cell-BM adhesion
	synonyms = {"cell-BM adhesion", "cell-membrane adhesion"};
	add_behavior(synonyms, cell_basement_membrane_adhesion, cell_basement_membrane_adhesion_base);

	// cell-BM repulsion 
	synonyms = {"cell-BM repulsion", "cell-membrane repulsion"};
	add_behavior(synonyms, cell_basement_membrane_repulsion, cell_basement_membrane_repulsion_base);

	// phagocytosis of apoptotic cell
	synonyms = {"phagocytose apoptotic cell",
						 "phagocytosis of apoptotic cell",
						 "phagocytosis of apoptotic cells"};
	add_behavior(synonyms, cell_phagocytose_apoptotic, cell_phagocytose_apoptotic_base);

	// phagocytosis of necrotic cell
	synonyms = {"phagocytose necrotic cell",
						 "phagocytosis of necrotic cell",
						 "phagocytosis of necrotic cells"};
	add_behavior(synonyms, cell_phagocytose_necrotic, cell_phagocytose_necrotic_base);

	// phagocytosis of other dead cell
	synonyms = {"phagocytose other dead cell",
						 "phagocytosis of other dead cell",
						 "phagocytosis of other dead cells"};
	add_behavior(synonyms, cell_phagocytose_other_dead, cell_phagocytose_other_dead_base);

	// phagocytosis of each live cell type 
	for( int i=0; i < n ; i++ )
	{
		synonyms = {"phagocytose " + cell_definitions_by_type[i]->name,
					"phagocytosis of " + cell_definitions_by_type[i]->name};
		behavior_value = [i](Cell *pCell) -> double&
		{ return cell_phagocytose_live_cell_type(pCell, i); };
		behavior_base_value = [i](Cell_Definition *pCD) -> double
		{ return cell_phagocytose_live_cell_type_base(pCD, i); };
		add_behavior(synonyms, behavior_value, behavior_base_value);
	}

	// attack of each live cell type 
	for( int i=0; i < n ; i++ )
	{
		behavior_value = [i](Cell *pCell) -> double&
		{ return cell_attack_type(pCell, i); };
		behavior_base_value = [i](Cell_Definition *pCD) -> double
		{ return cell_attack_type_base(pCD, i); };
		add_behavior("attack " + cell_definitions_by_type[i]->name, behavior_value, behavior_base_value);
	}

	// fusion 
	for( int i=0; i < n ; i++ )
	{
		behavior_value = [i](Cell *pCell) -> double&
		{ return cell_fuse_to_type(pCell, i); };
		behavior_base_value = [i](Cell_Definition *pCD) -> double
		{ return cell_fuse_to_type_base(pCD, i); };
		add_behavior("fuse to " + cell_definitions_by_type[i]->name, behavior_value, behavior_base_value);
	}

	// transformation / transition 
	for( int i=0; i < n ; i++ )
	{
		synonyms = {"transition to " + cell_definitions_by_type[i]->name,
					"transform to " + cell_definitions_by_type[i]->name};
		behavior_value = [i](Cell *pCell) -> double&
		{ return cell_transform_to_type(pCell, i); };
		behavior_base_value = [i](Cell_Definition *pCD) -> double
		{ return cell_transform_to_type_base(pCD, i); };
		add_behavior(synonyms, behavior_value, behavior_base_value);
	}

	// asymmetic division
	for( int i=0; i < n ; i++ )
	{
		behavior_value = [i](Cell *pCell) -> double&
		{ return cell_asymmetric_division_to_type(pCell, i); };
		behavior_base_value = [i](Cell_Definition *pCD) -> double
		{ return cell_asymmetric_division_to_type_base(pCD, i); };
		add_behavior("asymmetric division to " + cell_definitions_by_type[i]->name, behavior_value, behavior_base_value);
	}

	// is movable
	synonyms = {"is_movable", "is movable", "movable"};
	add_behavior(synonyms, cell_is_movable, cell_is_movable_base);

	// immunogenicity to each cell type
	for (int i = 0; i < n; i++)
	{
		behavior_value = [i](Cell *pCell) -> double&
		{ return cell_immunogenicity_to_type(pCell, i); };
		behavior_base_value = [i](Cell_Definition *pCD) -> double
		{ return cell_immunogenicity_to_type_base(pCD, i); };
		add_behavior("immunogenicity to " + cell_definitions_by_type[i]->name, behavior_value, behavior_base_value);
	}

	// cell attachment rate
	add_behavior("cell attachment rate", cell_attachment_rate, cell_attachment_rate_base);

	// cell detachment rate
	add_behavior("cell detachment rate", cell_detachment_rate, cell_detachment_rate_base);

	// maximum number of cell attachments
	add_behavior("maximum number of cell attachments", maximum_number_attachments, maximum_number_attachments_base);

	// attack damage rate
	add_behavior("attack damage rate", cell_attack_damage_rate, cell_attack_damage_rate_base);

	// attack duration
	add_behavior("attack duration", cell_attack_duration, cell_attack_duration_base);

	// damage rate
	add_behavior("damage rate", cell_damage_rate, cell_damage_rate_base);

	// damage repair rate
	add_behavior("damage repair rate", cell_damage_repair_rate, cell_damage_repair_rate_base);

	/* add new behaviors above this line */

    display_signal_dictionary(); 
    display_behavior_dictionary(); 
/*
	// now create empty SR models for each cell definition 

	for( int i=0 ; i < cell_definitions_by_index.size() ; i++ )
	{ create_SR_model( *cell_definitions_by_index[i] ); }
*/    

	return;
}

void add_signal( const std::string& name, SignalValue value )
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

void add_behavior(const std::string &name, BehaviorValue value, BehaviorBaseValue base_value)
{ return add_behavior(std::vector<std::string>{std::move(name)}, std::move(value), std::move(base_value)); }

void add_behavior(const std::vector<std::string> &synonyms, BehaviorValue value, BehaviorBaseValue base_value)
{
	auto behavior_ptr = std::make_shared<Behavior>(synonyms[0], std::move(value), std::move(base_value));
	for (auto name : synonyms)
	{ all_behaviors[name] = behavior_ptr; }
	return;
}

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

bool signal_exists(const std::string &signal_name)
{ return all_signals.find(signal_name) != all_signals.end(); }

bool behavior_exists(const std::string &behavior_name)
{ return all_behaviors.find(behavior_name) != all_behaviors.end(); }

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

std::vector<double> get_selected_signals( Cell* pCell , const std::vector<std::string> &signal_names )
{
	std::vector<double> signals( signal_names.size() , 0.0 );
	for( int i=0; i < signal_names.size(); i++ )
	{
		signals[i] = get_single_signal(pCell, signal_names[i]);
	}
	return signals; 
}

double get_single_signal(Cell *pCell, const std::string &name)
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

void set_single_behavior( Cell* pCell, const std::string& name , const double& parameter )
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

double get_single_behavior( Cell* pCell , const std::string& name )
{
	auto it = all_behaviors.find(name);
	if (it == all_behaviors.end())
	{
		std::cerr << "ERROR: Behavior '" << name << "' not found!" << std::endl;
		exit(-1);
	}
	return it->second->get_value(pCell);
}

std::vector<double> get_behaviors( Cell* pCell , std::vector<std::string> names )
{
	std::vector<double> parameters( names.size() , 0.0 ); 
	for( int n=0; n < names.size(); n++ )
	{ parameters[n] = get_single_behavior(pCell, names[n]); }
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
		base_behavior_map[behavior_name] = get_single_base_behavior(pCell, behavior_name);
		behaviors_got.push_back(behavior_name);
	}
	return base_behavior_map;
}

double get_single_base_behavior( Cell* pCell , const std::string& name )
{
	auto it = all_behaviors.find(name);
	if (it == all_behaviors.end())
	{
		std::cerr << "ERROR: Behavior '" << name << "' not found!" << std::endl;
		exit(-1);
	}
	return it->second->get_base_value(pCell);
}

double get_single_base_behavior(Cell_Definition *pCD, const std::string &name)
{
	auto it = all_behaviors.find(name);
	if (it == all_behaviors.end())
	{
		std::cerr << "ERROR: Behavior '" << name << "' not found!" << std::endl;
		exit(-1);
	}
	return it->second->get_base_value(pCD);
}

std::unordered_map<std::string, double> get_base_behaviors( Cell* pCell , std::vector<std::string> names )
{
	std::unordered_map<std::string, double> base_behavior_map;
	for( const auto& name : names )
	{ base_behavior_map[name] = get_single_base_behavior(pCell,name); }
	return base_behavior_map;
}

};
