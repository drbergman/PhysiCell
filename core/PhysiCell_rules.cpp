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

#include "./PhysiCell_rules.h"

namespace PhysiCell{

#ifndef __PhysiCell_rules_cpp__
#define __PhysiCell_rules_cpp__
#endif 

Hypothesis_Rule::Hypothesis_Rule()
{
	behavior = nullptr; 
	base_value = 1.0; 
	max_value = 10.0; 
	min_value = 0.1; 

	up_signals.resize(0);  
	up_half_maxes.resize(0);  
	up_hill_powers.resize(0);  
	up_applies_to_dead_cells.resize(0); 

	down_signals.resize(0);  
	down_half_maxes.resize(0);  
	down_hill_powers.resize(0);  
	down_applies_to_dead_cells.resize(0); 

	cell_type = "none"; 
	pCell_Definition = NULL; 

	return; 
}

std::string convert_bool_to_response( bool input )
{
	if( input )
	{ return "increases"; }
	return "decreases"; 
}

void Hypothesis_Rule::display( std::ostream& os )
{
	os << "For cell type " << cell_type << ": " << std::endl; 
	reduced_display(os);
	return;
}

void Hypothesis_Rule::reduced_display( std::ostream& os )
{
	std::string behavior_name = behavior->get_name();
	for( int j=0; j < down_signals.size(); j++ )
	for( const auto& down_signal : down_signals )
	{ os << down_signal->get_name() << " decreases " << behavior_name << std::endl; }
	for( const auto& up_signal : up_signals )
	{ os << up_signal->get_name() << " increases " << behavior_name << std::endl; }
	return; 
}

void Hypothesis_Rule::detailed_display( std::ostream& os )
{
	// os << "For cell type " << cell_type << ": " << std::endl; 
	std::string behavior_name = behavior->get_name();
	os << behavior_name << " is modulated from " << min_value << " to " << max_value << " with a base value of " << base_value << std::endl; 
	os << "--------------------------------------------------------" << std::endl; 
	for( int j=0; j < down_signals.size(); j++ )
	{
		os << "\t" << down_signals[j] << " decreases " << behavior_name
			<< " with half-max " << down_half_maxes[j] << " and Hill power " << down_hill_powers[j] << "."; 
		if( down_applies_to_dead_cells[j] == true )
		{ os << " Rule applies to dead cells."; }
		os << std::endl; 
	}
	for( int j=0; j < up_signals.size(); j++ )
	{
		os << "\t" << up_signals[j] << " increases " << behavior_name
			<< " with half-max " << up_half_maxes[j] << " and Hill power " << up_hill_powers[j] << "."; 
		if( up_applies_to_dead_cells[j] == true )
		{ os << " Rule applies to dead cells."; }
		os << std::endl; 
	}
	return; 
}

void Hypothesis_Rule::English_detailed_display( std::ostream& os )
{
	for( int j=0 ; j < down_signals.size(); j++ )
	{
		os << down_signals[j] << " decreases ";
		os << behavior->get_name() << " from " << base_value << " towards " ;
		os << min_value;
		os << " with a Hill response, with half-max " << down_half_maxes[j] ; 
		os << " and Hill power " << down_hill_powers[j] << ".";
		if( down_applies_to_dead_cells[j] == true )
		{ os << " Rule applies to dead cells."; }
		os << std::endl; 
	}
	for( int j=0 ; j < up_signals.size(); j++ )
	{
		os << up_signals[j] << " increases ";
		os << behavior->get_name() << " from " << base_value << " towards " ;
		os << max_value;
		os << " with a Hill response, with half-max " << up_half_maxes[j] ; 
		os << " and Hill power " << up_hill_powers[j] << ".";
		if( up_applies_to_dead_cells[j] == true )
		{ os << " Rule applies to dead cells."; }
		os << std::endl; 
	}
}

void Hypothesis_Rule::English_detailed_display_HTML( std::ostream& os )
{
	for( int j=0 ; j < down_signals.size(); j++ )
	{
		os << "<li>" << down_signals[j] << " decreases ";
		os << behavior->get_name() << " from " << base_value << " towards " ;
		os << min_value;
		os << " with a Hill response, with half-max " << down_half_maxes[j] ;
		os << " and Hill power " << down_hill_powers[j] << ".";
		if( down_applies_to_dead_cells[j] == true )
		{ os << " Rule applies to dead cells."; }
		os << "</li>" << std::endl; 
	}
	for( int j=0 ; j < up_signals.size(); j++ )
	{
		os << "<li>" << up_signals[j] << " increases ";
		os << behavior->get_name() << " from " << base_value << " towards " ;
		os << max_value;
		os << " with a Hill response, with half-max " << up_half_maxes[j] ;
		os << " and Hill power " << up_hill_powers[j] << ".";
		if( up_applies_to_dead_cells[j] == true )
		{ os << " Rule applies to dead cells."; }
		os << "</li>" << std::endl; 
	}
}

void Hypothesis_Rule::English_display( std::ostream& os )
{
	for( int j=0 ; j < down_signals.size(); j++ )
	{
		os << down_signals[j] << " decreases "; 
		os << behavior->get_name() << std::endl; 
	}
	for( int j=0 ; j < up_signals.size(); j++ )
	{
		os << up_signals[j] << " increases "; 
		os << behavior->get_name() << std::endl; 
	}
}

void Hypothesis_Rule::English_display_HTML( std::ostream& os )
{
	for( int j=0 ; j < down_signals.size(); j++ )
	{
		os << "<li>" << down_signals[j] << " decreases "; 
		os << behavior->get_name() << "</li>" << std::endl; 
	}
	for( int j=0 ; j < up_signals.size(); j++ )
	{
		os << "<li>" << up_signals[j] << " increases "; 
		os << behavior->get_name() << "</li>" << std::endl; 
	}
}

void Hypothesis_Rule::add_signal(std::string signal, std::string response, double half_max, double hill_power, bool use_for_dead)
{
    // check: is this a valid signal? (is it in the dictionary?)
	if (!signal_exists(signal))
    {
        std::cout << "Error! Attempted to add signal " << signal << " which is not in the dictionary." << std::endl; 
        std::cout << "Either fix your model or add the missing signal to the simulation." << std::endl; 

		std::cout << "\t\tSee possible fixes at https://github.com/physicell-training/PhysiCell_common_errors\n\n"; 

        exit(-1); 
    }

	// check to see if the signal and response already there 
	bool is_up = false; // true if up-regulate, false if down
	if( response == "increase" || response == "increases" || response == "promotes" )
	{ is_up = true; }
	int n = find_signal(signal, is_up); 

	// if so, then just warn and exit.  
	if( n > -1 )
	{
		std::cout << "Error! Signal " << signal << " and Response " << response << " were already part of the rule." << 
		std::endl; 

		std::cout << "\t\tSee possible fixes at https://github.com/physicell-training/PhysiCell_common_errors\n\n"; 

		exit(-1); 
	}

	// add the signal; 
	if( is_up == true )
	{
		up_signals.push_back( all_signals[signal].get() ); 
		up_half_maxes.push_back( half_max ); 
		up_hill_powers.push_back( hill_power ); 
		up_applies_to_dead_cells.push_back( use_for_dead ); 
	}
	else
	{
		down_signals.push_back( all_signals[signal].get() ); 
		down_half_maxes.push_back( half_max ); 
		down_hill_powers.push_back( hill_power ); 
		down_applies_to_dead_cells.push_back( use_for_dead ); 
	}
	return; 
}

double Hypothesis_Rule::evaluate(std::vector<double> down_signal_values, std::vector<double> up_signal_values)
{
	// up-regulation part
	double HU = multivariate_Hill_response_function(up_signal_values, up_half_maxes, up_hill_powers);
	double U = base_value + (max_value - base_value) * HU;

	// then the down-regulation part
	double DU = multivariate_Hill_response_function(down_signal_values, down_half_maxes, down_hill_powers);
	double output = U + (min_value - U) * DU;

	return output;
}

double Hypothesis_Rule::evaluate( Cell* pCell )
{
	// construct signal vector 
	bool dead = pCell->phenotype.death.dead;
	bool apply = false;
	std::vector<double> down_signal_values( down_signals.size() , 0.0 ); 
	for (size_t i = 0; i < down_signals.size(); i++)
	{
		if(dead && !down_applies_to_dead_cells[i]) { continue; }
		down_signal_values.push_back(down_signals[i]->get_value(pCell));
		apply = true;
	}
	std::vector<double> up_signal_values( up_signals.size() , 0.0 ); 
	for (size_t i = 0; i < up_signals.size(); i++)
	{
		if(dead && !up_applies_to_dead_cells[i]) { continue; }
		up_signal_values.push_back(up_signals[i]->get_value(pCell));
		apply = true;
	}

	if (!apply)
	{ return behavior->get_value(pCell); } // leave the behavior value unchanged

	return evaluate( down_signal_values, up_signal_values ); 
}

void Hypothesis_Rule::apply( Cell* pCell )
{
	// evaluate the rule 
	double param = evaluate( pCell ); 

	// apply it ot the appropriate behavior 
	behavior->set_value(pCell, param);

	return; 
}

void Hypothesis_Rule::sync_to_cell_definition( Cell_Definition* pCD )
{
	if( pCD == NULL )
	{ return; }

	cell_type = pCD->name; 
	pCell_Definition = pCD; 

	// sync base behavior 
	base_value = behavior->get_base_value(pCD);

	return; 
}

void Hypothesis_Rule::sync_to_cell_definition( std::string cell_name )
{ return sync_to_cell_definition( find_cell_definition(cell_name) ); }

int Hypothesis_Rule::find_signal( std::string name, bool is_up )
{
	std::vector<Signal*> &signals = is_up ? up_signals : down_signals;
	for (size_t i = 0; i < signals.size(); i++) {
		if (signals[i]->get_name() == name) {
			return i; // Return the index of the found signal
		}
	}
	return -1; // Signal not found
}

Hypothesis_Ruleset::Hypothesis_Ruleset()
{
	cell_type = "none"; 
	pCell_Definition = NULL; 

	rules.resize(0); 
	rules_map.clear(); 

	return; 
}

void Hypothesis_Ruleset::display( std::ostream& os )
{
	os << "Behavioral rules for cell type " << cell_type << ":" << std::endl; 
	os << "===================================================" << std::endl; 
	for( int i=0; i < rules.size() ; i++ )
	{ rules[i]->reduced_display(os); }
	os << std::endl; 
	return; 
}

void Hypothesis_Ruleset::detailed_display( std::ostream& os )
{
	os << "Behavioral rules for cell type " << cell_type << ":" << std::endl; 
	os << "===================================================" << std::endl; 
	for( int i=0; i < rules.size() ; i++ )
	{ rules[i]->detailed_display(os); }
	os << std::endl; 
	return; 
}


void Hypothesis_Ruleset::sync_to_cell_definition( Cell_Definition* pCD )
{
	pCell_Definition = pCD; 
	cell_type = pCD->name; 

	for( int i=0; i < rules.size(); i++ )
	{ rules[i]->sync_to_cell_definition(pCD); }

	return; 
}

Hypothesis_Rule* Hypothesis_Ruleset::add_behavior( std::string behavior , double min_behavior, double max_behavior )
{
    // check: is this a valid signal? (is it in the dictionary?)
    if( !behavior_exists(behavior) )
    {
        std::cout << "ERROR! Attempted to add behavior " << behavior << " which is not in the dictionary." << std::endl; 
        std::cout << "Either fix your model or add the missing behavior to the simulation." << std::endl; 

		std::cout << "\t\tSee possible fixes at https://github.com/physicell-training/PhysiCell_common_errors\n\n";

        exit(-1); 
    }

	// first, check. Is there already a ruleset? 
	auto search = rules_map.find( behavior ); 

		// if not, add it 
	if( search == rules_map.end() )
	{
		Hypothesis_Rule *pHR = new Hypothesis_Rule;

		pHR->set_behavior(all_behaviors[behavior]);

		std::cout << "Adding behavior " << behavior << " to ruleset for cell type " << cell_type << std::endl;
		std::cout << "  It's name is set to " << pHR->get_behavior()->get_name() << std::endl;

		pHR->sync_to_cell_definition( pCell_Definition ); 

		pHR->min_value = min_behavior; 
		pHR->max_value = max_behavior;

		rules.push_back(pHR);
		rules_map[ behavior ] = pHR; 

		return pHR; 
	}

		// otherwise, edit it 
	Hypothesis_Rule* pHR = search->second; 

	/*
		// March 28 2023 fix  : let's not overwrite eixsting values
	pHR->min_value = min_behavior; 
	pHR->max_value = max_behavior; 
	*/

	return pHR; 
}

Hypothesis_Rule* Hypothesis_Ruleset::add_behavior( std::string behavior )
{ 
	double min_behavior = 9e99; // Min behaviour high value
	double max_behavior = -9e99; // Max behaviour low value
	return Hypothesis_Ruleset::add_behavior( behavior, min_behavior, max_behavior );
}

void Hypothesis_Ruleset::sync_to_cell_definition( std::string cell_name )
{ return sync_to_cell_definition( find_cell_definition(cell_name) ); }

Hypothesis_Rule* Hypothesis_Ruleset::find_behavior( std::string name )
{
    auto search = rules_map.find( name); 
	if( search == rules_map.end() )
	{
		// std::cout << "Warning! Ruleset does not contain " << name << std::endl; 
		// std::cout << "         Returning NULL." << std::endl; 
		return NULL; 
	}

	return search->second; 
}

Hypothesis_Rule& Hypothesis_Ruleset::operator[]( std::string name )
{
	Hypothesis_Rule* pHR = find_behavior(name);
	return *pHR; 
} 

void Hypothesis_Ruleset::apply( Cell* pCell )
{
	for( int n=0; n < rules.size() ; n++ )
	{ rules[n]->apply( pCell );  }
	return; 
}

std::unordered_map< Cell_Definition* , Hypothesis_Ruleset > hypothesis_rulesets; 

void add_hypothesis_ruleset( Cell_Definition* pCD )
{
	auto search = hypothesis_rulesets.find( pCD );
	if( search == hypothesis_rulesets.end() )
	{
		Hypothesis_Ruleset HRS; 
		HRS.sync_to_cell_definition( pCD ); 
		hypothesis_rulesets[pCD] = HRS; 
	}
	return; 
}

void intialize_hypothesis_rulesets( void )
{
	hypothesis_rulesets.clear(); // empty(); 

	for( int n; n < cell_definitions_by_index.size() ; n++ )
	{
		Cell_Definition* pCD = cell_definitions_by_index[n]; 
		add_hypothesis_ruleset(pCD); 
	}

	return; 
}

Hypothesis_Ruleset& access_ruleset( Cell_Definition* pCD )
{ return hypothesis_rulesets[pCD]; }

Hypothesis_Ruleset* find_ruleset( Cell_Definition* pCD )
{ return &(hypothesis_rulesets[pCD]); }

void display_hypothesis_rulesets( std::ostream& os )
{
	for( int n=0 ; n < cell_definitions_by_index.size() ; n++ )
	{ hypothesis_rulesets[ cell_definitions_by_index[n] ].display( os ); }

	return; 
}

void detailed_display_hypothesis_rulesets( std::ostream& os )
{
	for( int n=0 ; n < cell_definitions_by_index.size() ; n++ )
	{ hypothesis_rulesets[ cell_definitions_by_index[n] ].detailed_display( os ); }

	return; 
}

void add_rule( std::string cell_type, std::string signal, std::string behavior , std::string response , double half_max, double hill_power, bool use_for_dead )
{
	Cell_Definition* pCD = find_cell_definition(cell_type); 
    if( !pCD )
    {
        std::cout << "Warning: Attempted to add rule for " << cell_type  
            << ", but no cell definition found for this type." << std::endl; 
        exit(-1); 
    }

	Hypothesis_Ruleset* pHRS = find_ruleset( pCD ); 
    if( !pHRS )
    {
        std::cout << "Warning: Attempted to add rule for " << cell_type  
            << ", but no hypothesis ruleset found for this type." << std::endl; 
        exit(-1); 
    }

	if( pHRS->find_behavior(behavior) )
	{
		if( (*pHRS)[behavior].get_behavior()->get_name() != behavior )
		{ (*pHRS)[behavior].set_behavior(all_behaviors[behavior]); std::cout << "wha?" << std::endl; }
	}

	pHRS->add_behavior(behavior);

	(*pHRS)[behavior].add_signal(signal, response, half_max, hill_power, use_for_dead);

	return; 
}

void set_behavior_parameters( std::string cell_type, std::string behavior, 
   double min_value, double max_value )
{
	Cell_Definition* pCD = find_cell_definition( cell_type ); 
    if( !pCD )
    {
        std::cout << "Warning: Attempted to set parameters for " 
          << behavior << " in " << cell_type  
            << ", but the cell definition is not found." << std::endl; 
        exit(-1);         
    }

	if( find_ruleset(pCD) == NULL )
	{
        std::cout << "Warning: Attempted to set parameters for " 
          << behavior << " in " << cell_type  
            << ", but there is no hypothesis ruleset for this cell type." << std::endl; 
        exit(-1);              
    }

	if( hypothesis_rulesets[pCD].find_behavior(behavior) == NULL )
	{
        std::cout << "Warning: Attempted to set parameters for " 
          << behavior << " in " << cell_type  
            << ", but there is no rules for this behavior for this cell type." << std::endl; 
        exit(-1);              
    }

	hypothesis_rulesets[pCD][behavior].min_value = min_value; 
	hypothesis_rulesets[pCD][behavior].max_value = max_value; 

	return;
}

void set_behavior_parameters( std::string cell_type, std::string behavior, 
   double min_value, double base_value , double max_value )
{
	Cell_Definition* pCD = find_cell_definition( cell_type ); 
    if( !pCD )
    {
        std::cout << "Warning: Attempted to set parameters for " 
          << behavior << " in " << cell_type  
            << ", but the cell definition is not found." << std::endl; 
        exit(-1);         
    }

	if( find_ruleset(pCD) == NULL )
	{
        std::cout << "Warning: Attempted to set parameters for " 
          << behavior << " in " << cell_type  
            << ", but there is no hypothesis ruleset for this cell type." << std::endl; 
        exit(-1);              
    }


	if( hypothesis_rulesets[pCD].find_behavior(behavior) == NULL )
	{
        std::cout << "Warning: Attempted to set parameters for " 
          << behavior << " in " << cell_type  
            << ", but there is no rules for this behavior for this cell type." << std::endl; 
        exit(-1);              
    }

	if ( min_value < hypothesis_rulesets[pCD][behavior].min_value )
	{ hypothesis_rulesets[pCD][behavior].min_value = min_value; }
	if ( max_value > hypothesis_rulesets[pCD][behavior].max_value )
	{ hypothesis_rulesets[pCD][behavior].max_value = max_value; } 
	hypothesis_rulesets[pCD][behavior].base_value = base_value; 

	return;
}


void set_behavior_base_value( std::string cell_type, std::string behavior, double base_value )
{
	Cell_Definition* pCD = find_cell_definition( cell_type ); 
    if( !pCD )
    {
        std::cout << "Warning: Attempted to set base parameter for " 
          << behavior << " in " << cell_type  
            << ", but the cell definition is not found." << std::endl; 
        exit(-1);         
    }

	if( find_ruleset(pCD) == NULL )
    {
        std::cout << "Warning: Attempted to set base parameter for " 
          << behavior << " in " << cell_type  
            << ", but no hypothesis ruleset found for this cell type." << std::endl; 
        exit(-1);         
    }

	if( hypothesis_rulesets[pCD].find_behavior(behavior) == NULL )
    {
        std::cout << "Warning: Attempted to set base parameter for " 
          << behavior << " in " << cell_type  
            << ", but no rule for this behavior found for this cell type." << std::endl; 
        exit(-1);         
    }

	hypothesis_rulesets[pCD][behavior].base_value = base_value; 

	return;
}

void set_behavior_min_value( std::string cell_type, std::string behavior, double min_value )
{
	Cell_Definition* pCD = find_cell_definition( cell_type ); 
    if( !pCD )
    {
        std::cout << "Warning: Attempted to set base parameter for " 
          << behavior << " in " << cell_type  
            << ", but the cell definition is not found." << std::endl; 
        exit(-1);         
    }

	if( find_ruleset(pCD) == NULL )
    {
        std::cout << "Warning: Attempted to set base parameter for " 
          << behavior << " in " << cell_type  
            << ", but no hypothesis ruleset found for this cell type." << std::endl; 
        exit(-1);         
    }

	if( hypothesis_rulesets[pCD].find_behavior(behavior) == NULL )
    {
        std::cout << "Warning: Attempted to set base parameter for " 
          << behavior << " in " << cell_type  
            << ", but no rule for this behavior found for this cell type." << std::endl; 
        exit(-1);         
    }
	
	if ( min_value < hypothesis_rulesets[pCD][behavior].min_value )
	{ hypothesis_rulesets[pCD][behavior].min_value = min_value; }

	return;
}

void set_behavior_max_value( std::string cell_type, std::string behavior, double max_value )
{
	Cell_Definition* pCD = find_cell_definition( cell_type ); 
    if( !pCD )
    {
        std::cout << "Warning: Attempted to set base parameter for " 
          << behavior << " in " << cell_type  
            << ", but the cell definition is not found." << std::endl; 
        exit(-1);         
    }

	if( find_ruleset(pCD) == NULL )
    {
        std::cout << "Warning: Attempted to set base parameter for " 
          << behavior << " in " << cell_type  
            << ", but no hypothesis ruleset found for this cell type." << std::endl; 
        exit(-1);         
    }

	if( hypothesis_rulesets[pCD].find_behavior(behavior) == NULL )
    {
        std::cout << "Warning: Attempted to set base parameter for " 
          << behavior << " in " << cell_type  
            << ", but no rule for this behavior found for this cell type." << std::endl; 
        exit(-1);         
    }

	if ( max_value > hypothesis_rulesets[pCD][behavior].max_value )
	{ hypothesis_rulesets[pCD][behavior].max_value = max_value; }
	
	return;
}



void apply_ruleset( Cell* pCell )
{
	Cell_Definition* pCD = find_cell_definition( pCell->type_name ); 
	hypothesis_rulesets[pCD].apply( pCell );
	return; 
}

void rule_phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	apply_ruleset( pCell );

	// by default, rules only apply to live cells
	// 

/*
	// safety checks for dead cells 
	if( get_single_signal(pCell,"dead") > 0.11 )
	{
		// can't die twice
		set_single_behavior(pCell,"apoptosis",0.0);
		set_single_behavior(pCell,"necrosis",0.0);

		// can't cycle 
		set_single_behavior(pCell,"cycle entry",0.0);

		// can't crawl 
		set_single_behavior(pCell,"migration speed",0.0);
	}
*/

	return; 	
}

/* add these to core */ 
/*
std::vector<double> linear_response_to_Hill_parameters( double s0, double s1 )
{
	static double tol = 0.1; 
	static double param1 = (1-tol)/tol; 
	static double param2 = log(param1); 

	// half max, then hill power 
	double hm = 0.5* (s0+s1); 

	// hp so that H(s1) ~ (1-tol)
	double hp = round( param2 / log(s1/hm) ); 

	std::vector<double> output = { hm , hp }; 

	return output; 
}

std::vector<double> Hill_response_to_linear_parameters( double half_max , double Hill_power )
{
	static double tol = 0.1; 
	static double param1 = (1-tol)/tol; 
	double param2 = pow( param1 , 1.0/ Hill_power ); 

	// s1 such that H(s1) ~ (1-tol)
	double s1 = half_max * param2; 

	// s0 for symmetry
	double s0 = 2*half_max -s1; 
	if( s0 < 0 )
	{ s0 = 0.0; }

	std::vector<double> output = {s0,s1}; 

	return output; 
}
*/

void split_csv( std::string input , std::vector<std::string>& output , char delim )
{
	output.resize(0); 

	std::istringstream is(input);
	std::string part;
	while( getline(is, part, delim ) )
	{ output.push_back(part); }

	return; 
}

std::string csv_strings_to_English_v1( std::vector<std::string> strings , bool include_cell_header )
{
	std::string output = ""; 

	if( include_cell_header )
	{
		output += "In "; 
		output += strings[0]; 
		output += " cells:\n\t"; // In {cell type X} cells: 
	}

// malignant epithelial,oxygen,decreases,necrosis,2.80E-03,0,decreases,3.75,8,0

	output += strings[1] ; // {signal}
	output += " ";

	output += strings[2] ; // {increases/decreases}
	output += " ";

	output += strings[3] ; // {behavior}


	output += " from "; // {base}
	output += strings[4]; 

	output += " towards ";
	output += strings[5]; 

	output += " with a Hill response, with half-max "; 
	output += strings[6]; 

	output += " and Hill power "; 
	output += strings[7]; 

	output += "."; 
	bool use_when_dead = false; 
	char start_char = toupper( strings[8][0] ); 
	if( start_char == 'T' || start_char == '1' )
	{ output += " Rule applies to dead cells."; }
	
	return output; 
}

std::string csv_strings_to_English_v3( std::vector<std::string> strings , bool include_cell_header )
{
	std::string output = ""; 

	if( include_cell_header )
	{
		output += "In "; 
		output += strings[0]; 
		output += " cells:\n\t"; // In {cell type X} cells: 
	}

	output += strings[1] ; // {signal}
	output += " ";

	output += strings[2] ; // {increases/decreases}
	output += " ";

	output += strings[3] ; // {behavior}


//	output += " from "; // {base}
//	output += strings[4]; 

	output += " towards ";
	output += strings[4]; 

	output += " with a Hill response, with half-max "; 
	output += strings[5]; 

	output += " and Hill power "; 
	output += strings[6]; 

	output += "."; 
	bool use_when_dead = false; 
	char start_char = toupper( strings[7][0] ); 
	if( start_char == 'T' || start_char == '1' )
	{ output += " Rule applies to dead cells."; }
	
	return output; 
}


std::string csv_strings_to_English( std::vector<std::string> strings , bool include_cell_header )
{
	std::string output = ""; 

	if( include_cell_header )
	{
		output += "In "; 
		output += strings[0]; 
		output += " cells:\n\t"; // In {cell type X} cells: 
	}

	output += strings[5] ; // {signal}
	output += " ";

	output += strings[6] ; // {increases/decreases}
	output += " ";

	output += strings[1] ; // {behavior}


	output += " from "; // {base}
	output += strings[3]; 

	output += " towards ";
	if( strings[6][0] == 'i' || strings[6][0] == 'I' )
	{ output+= strings[4];  }
	else
	{ output+= strings[2];  }

	output += " with a Hill response, with half-max "; 
	output += strings[7]; 

	output += " and Hill power "; 
	output += strings[8]; 

	output += "."; 
	bool use_when_dead = false; 
	char start_char = toupper( strings[9][0] ); 
	if( start_char == 'T' || start_char == '1' )
	{ output += " Rule applies to dead cells."; }
	
	return output; 
}


std::string csv_strings_to_English_HTML( std::vector<std::string> strings , bool include_cell_header )
{
	std::string output = "<p>"; 

	if( include_cell_header )
	{
		output += "In "; 
		output += strings[0]; 
		output += " cells:<br>\n"; // In {cell type X} cells: 
	}

	output += "&nbsp;";
	output += strings[5] ; // {signal}
	output += " ";

	output += strings[6] ; // {increases/decreases}
	output += " ";

	output += strings[1] ; // {behavior}

	output += " from "; // {base}
	output += strings[3]; 

	output += " towards ";
	if( strings[6][0] == 'i' || strings[6][0] == 'I' )
	{ output+= strings[4];  }
	else
	{ output+= strings[2];  }

	output += " with a Hill response, with half-max "; 
	output += strings[7]; 

	output += " and Hill power "; 
	output += strings[8]; 

	output += "."; 
	bool use_when_dead = false; 
	char start_char = toupper( strings[9][0] ); 
	if( start_char == 'T' || start_char == '1' )
	{ output += " Rule applies to dead cells."; }
	output += "\n</p>\n"; 
	
	return output; 
}


/*
0              1          2        3           4         5       6          7         8           9 
 Cell type, behavior, min value, base value, max value, signal, direction, half-max, Hill power, dead
*/ 

void parse_csv_rule_v0( std::vector<std::string> input )
{
	// if it's wrong length, skip 
	bool skip = false; 
	if( input.size() != 10 )
	{ skip = true; }
	// if any empty strings, skip
	for( int n=0 ; n < input.size() ; n++ )
	{
		if( input[n].empty() == true )
		{ skip = true; }
	}
	if( skip == true )
	{
		std::cout << "Warning: Misformed rule (likely from an empty rules file). Skipping." << std::endl; 
		return; 
	}

	std::string temp = csv_strings_to_English( input , false ); 

	// string portions of the rule
	std::string cell_type = input[0]; 
	std::string behavior = input[1]; 
	std::string signal = input[5]; 
	std::string response = input[6]; 

	// numeric portions of the rule 
	double min_value  = std::atof( input[2].c_str() );
	double base_value = std::atof( input[3].c_str() );
	double max_value  = std::atof( input[4].c_str() ); 

	double half_max  = std::atof( input[7].c_str() );
	double hill_power = std::atof( input[8].c_str() );
	bool use_for_dead = (bool) std::atof( input[9].c_str() ); 

	std::cout << "Adding rule for " << cell_type << " cells:\n\t"; 
	std::cout << temp << std::endl;

	add_rule(cell_type, signal, behavior, response, half_max, hill_power, use_for_dead);

	// compare to base behavior value in cell def for discrepancies 
	Cell_Definition* pCD = find_cell_definition(cell_type); 
	double ref_base_value = get_single_base_behavior(pCD,behavior); 
	if( fabs(ref_base_value-base_value) > 1e-15 )
	{
		std::cout << "Error: Base value for " << behavior << " in cell type " << cell_type  << std::endl 
				  << "       has base value " << base_value << " in the rule, " // << std::endl 
		 		  << "but base value " << ref_base_value << " in the cell definition." << std::endl
				  << "       Fix this discrepancy to continue the model." << std::endl << std::endl; 
				  exit(-1); 
	}

	set_behavior_parameters(cell_type,behavior,min_value,base_value,max_value); 

	return;  
}

void parse_csv_rule_v0( std::string input )
{
	std::vector<std::string> tokenized_string; 
	split_csv( input , tokenized_string , ','); 

	// Make sure it was truly comma-separated. 
	// If not, try tab.
	if( tokenized_string.size() != 10 )
	{ split_csv( input , tokenized_string , '\t');  }

	return parse_csv_rule_v0( tokenized_string ); 
}

void parse_csv_rules_v0( std::string filename )
{
	std::fstream fs( filename, std::ios::in );
	if( !fs )
	{
		std::cout << "Warning: Rules file " << filename << " failed to open." << std::endl; 
		return; 
	}

	std::cout << "Processing rules in file " << filename << " ... " << std::endl; 

	while( fs.eof() == false )
	{
		std::string line; 	
		std::getline( fs , line, '\n'); 
		if( line.size() > 0 )
		{ parse_csv_rule_v0(line); }
	}

	fs.close(); 

	std::cout << "Done!" << std::endl << std::endl; 

	return; 
}

/* v1 work */

/*
 v0:::
0              1          2        3           4         5       6          7         8           9 
 Cell type, behavior, min value, base value, max value, signal, direction, half-max, Hill power, dead


 v1:::
 0          1       2          3         4           5                   6         7           8
 Cell type, signal, direction, behavior, base value, max response value, half-max, Hill power, applies to dead? 
*/ 

void parse_csv_rule_v1( std::vector<std::string> input )
{
	// if it's wrong length, skip 
	bool skip = false; 
	if( input.size() != 9 )
	{ skip = true; }
	// if any empty strings, skip
	for( int n=0 ; n < input.size() ; n++ )
	{
		if( input[n].empty() == true )
		{ skip = true; }
	}
	if( skip == true )
	{
		std::cout << "Warning: Misformed rule (likely from an empty rules file). Skipping." << std::endl; 

		for( int n=0 ; n < input.size(); n++ )
		{
			std::cout << n << " : " << input[n] << std::endl; 
		}

		return; 
	}

	std::string temp = csv_strings_to_English_v1( input , false ); // need a v1 version of this

	// string portions of the rule
	std::string cell_type = input[0]; 
	std::string signal = input[1]; 
	std::string response = input[2]; 
	std::string behavior = input[3]; 

	// numeric portions of the rule 
	// double min_value  = std::atof( input[2].c_str() );

	double base_value = std::atof( input[4].c_str() );
	double max_response = std::atof( input[5].c_str() ); 

	// hmm from here 
	// double max_value  = std::atof( input[4].c_str() ); 

	double half_max  = std::atof( input[6].c_str() );
	double hill_power = std::atof( input[7].c_str() );
	bool use_for_dead = (bool) std::atof( input[8].c_str() ); 

	std::cout << "Adding rule for " << cell_type << " cells:\n\t"; 
	std::cout << temp << std::endl;

	add_rule(cell_type, signal, behavior, response, half_max, hill_power, use_for_dead);

	// compare to base behavior value in cell def for discrepancies 
	Cell_Definition* pCD = find_cell_definition(cell_type); 
	double ref_base_value = get_single_base_behavior(pCD,behavior); 
	if( fabs(ref_base_value-base_value) > 1e-15 )
	{
		std::cout << "Error: Base value for " << behavior << " in cell type " << cell_type  << std::endl 
				  << "       has base value " << base_value << " in the rule, " // << std::endl 
		 		  << "but base value " << ref_base_value << " in the cell definition." << std::endl
				  << "       Fix this discrepancy to continue the model." << std::endl << std::endl; 
				  exit(-1); 
	}

	set_behavior_base_value(cell_type,behavior,base_value);

	if( response == "increases")
	{ 
		set_behavior_min_value(cell_type,behavior,ref_base_value); 
		set_behavior_max_value(cell_type,behavior,max_response);
	}
	else
	{ 
		set_behavior_min_value(cell_type,behavior,max_response); 
		set_behavior_max_value(cell_type,behavior,ref_base_value);
	}

	return;  
}

void parse_csv_rule_v1( std::string input )
{
	std::vector<std::string> tokenized_string; 
	split_csv( input , tokenized_string , ','); 

	// Make sure it was truly comma-separated. 
	// If not, try tab.
	if( tokenized_string.size() != 9 )
	{ split_csv( input , tokenized_string , '\t');  }

	// check for comment 
	if(tokenized_string[0][0] == '/' && tokenized_string[0][1] == '/' )
	{ std::cout << "Skipping commented rule (" << input << ")" << std::endl; return; }	

	return parse_csv_rule_v1( tokenized_string ); 
}

void parse_csv_rules_v1( std::string filename )
{
	std::fstream fs( filename, std::ios::in );
	if( !fs )
	{
		std::cout << "Warning: Rules file " << filename << " failed to open." << std::endl; 
		return; 
	}

	std::cout << "Processing rules in file " << filename << " ... " << std::endl; 

	while( fs.eof() == false )
	{
		std::string line; 	
		std::getline( fs , line, '\n'); 
		if( line.size() > 0 )
		{ parse_csv_rule_v1(line); }
	}

	fs.close(); 

	std::cout << "Done!" << std::endl << std::endl; 

	return; 
}


/* end of v1 work */

/* v2 work */

/*
 v0:::
0              1          2        3           4         5       6          7         8           9 
 Cell type, behavior, min value, base value, max value, signal, direction, half-max, Hill power, dead


 v1:::
 0          1       2          3         4           5                   6         7           8
 Cell type, signal, direction, behavior, base value, max response value, half-max, Hill power, applies to dead? 

 v2:::
 0          1       2          3         4                   5         6           7           
 Cell type, signal, direction, behavior, max response value, half-max, Hill power, applies to dead?  
*/ 

void parse_csv_rule_v3( std::vector<std::string> input )
{
	// if it's wrong length, skip 
	bool skip = false; 
	if( input.size() != 8 )
	{ skip = true; }
	// if any empty strings, skip
	for( int n=0 ; n < input.size() ; n++ )
	{
		if( input[n].empty() == true )
		{ skip = true; }
	}
	if( skip == true )
	{
		std::cout << "Warning: Misformed rule (likely from an empty rules file). Skipping." << std::endl; 

		for( int n=0 ; n < input.size(); n++ )
		{
			std::cout << n << " : " << input[n] << std::endl; 
		}

		return; 
	}

	std::string temp = csv_strings_to_English_v3( input , false ); // need a v1 version of this

	// string portions of the rule
	std::string cell_type = input[0]; 
	std::string signal = input[1]; 
	std::string response = input[2]; 
	std::string behavior = input[3]; 

	// numeric portions of the rule 
	// double min_value  = std::atof( input[2].c_str() );

	// double base_value = std::atof( input[4].c_str() );
	double max_response = std::atof( input[4].c_str() ); 

	// hmm from here 
	// double max_value  = std::atof( input[4].c_str() ); 

	double half_max  = std::atof( input[5].c_str() );
	double hill_power = std::atof( input[6].c_str() );
	bool use_for_dead = (bool) std::atof( input[7].c_str() ); 

	std::cout << "Adding rule for " << cell_type << " cells:\n\t"; 
	std::cout << temp << std::endl;

	add_rule(cell_type, signal, behavior, response, half_max, hill_power, use_for_dead);

	// compare to base behavior value in cell def for discrepancies 
	Cell_Definition* pCD = find_cell_definition(cell_type); 
	double ref_base_value = get_single_base_behavior(pCD,behavior); 

	set_behavior_base_value(cell_type,behavior,ref_base_value);
	if( response == "increases")
	{ 
		set_behavior_min_value(cell_type,behavior,ref_base_value); 
		set_behavior_max_value(cell_type,behavior,max_response);
	}
	else
	{ 
		set_behavior_min_value(cell_type,behavior,max_response); 
		set_behavior_max_value(cell_type,behavior,ref_base_value);
	}
	return;  
}

void parse_csv_rule_v3( std::string input )
{
	std::vector<std::string> tokenized_string; 
	split_csv( input , tokenized_string , ','); 

	// Make sure it was truly comma-separated. 
	// If not, try tab.
	if( tokenized_string.size() != 8 )
	{ split_csv( input , tokenized_string , '\t');  }

	// check for comment 
	if(tokenized_string[0][0] == '/' && tokenized_string[0][1] == '/' )
	{ std::cout << "Skipping commented rule (" << input << ")" << std::endl; return; }	

	return parse_csv_rule_v3( tokenized_string ); 
}

void parse_csv_rules_v3( std::string filename )
{
	std::fstream fs( filename, std::ios::in );
	if( !fs )
	{
		std::cout << "Warning: Rules file " << filename << " failed to open." << std::endl; 
		return; 
	}

	std::cout << "Processing rules in file " << filename << " ... " << std::endl; 

	while( fs.eof() == false )
	{
		std::string line; 	
		std::getline( fs , line, '\n'); 
		if( line.size() > 0 )
		{ parse_csv_rule_v3(line); }
	}

	fs.close(); 

	std::cout << "Done!" << std::endl << std::endl; 

	return; 
}


/* end of v2 work */

// needs fixing
void parse_rules_from_pugixml( void )
{
	pugi::xml_node node = physicell_config_root.child( "cell_rules" ); 
	if( !node )
	{ 
		std::cout << "Error: Could not find <cell_rules> section of XML config file." << std::endl 
				 <<  "       Cannot parse cell rules, so disabling." << std::endl; 

		PhysiCell_settings.rules_enabled = false; 
		return; 
	}

	// find the set of rulesets 
	node = node.child( "rulesets" ); 
	if( !node )
	{ 
		std::cout << "Error: Could not find <rulesets> in the <cell_rules> section of XML config file." << std::endl 
				 <<  "       Cannot parse cell rules, so disabling." << std::endl; 

		PhysiCell_settings.rules_enabled = false; 
		return; 
	}
	// find the first ruleset 
	node = node.child( "ruleset");
	if( !node )
	{ 
		std::cout << "Error: Could not find any <ruleset> in the <rulesets> section of XML config file." << std::endl 
				 <<  "       Cannot parse cell rules, so disabling." << std::endl; 

		PhysiCell_settings.rules_enabled = false; 
		return; 
	}

	while( node )
	{
		std::cout << node.name() << std::endl;
		if( node.attribute("enabled").as_bool() == true )
		{ 
			std::string folder = xml_get_string_value( node, "folder" ); 
			std::string filename = xml_get_string_value( node, "filename" ); 
			std::string input_filename = folder + "/" + filename; 

			std::cout << "\tProcessing ruleset in " << input_filename << " ... " << std::endl; 
			std::string format = node.attribute("format").as_string(); 
			std::string protocol = node.attribute("protocol").as_string(); 
			double version = node.attribute("version").as_double(); 

			bool done = false; 

			if( format == "CSV" || format == "csv" )  
			{ 			
				if( version < 1.0 )
				{
					std::cout << "\tFormat: CSV (prototype version)" << std::endl; 

					// parse_csv_rules_v0( input_filename ); // parse all rules in a CSV file 

					PhysiCell_settings.rules_enabled = true; 

					std::cout << "\t\t**Error: Version < 3 no longer supported.\n\n"; 
					std::cout << "\t\tSee possible fixes at https://github.com/physicell-training/PhysiCell_common_errors\n\n"; 
					exit(-1); 

					done = true; 
				}
			
				if(version >= 1.0 - 1e-10 && version < 2.0 - 1e-10 && protocol == "CBHG" && done == false )
				{
					std::cout << "\tFormat: CSV (version " << version << ")" << std::endl; 

					// parse_csv_rules_v1( input_filename ); // parse all rules in a CSV file 

					PhysiCell_settings.rules_enabled = true; 

					std::cout << "\t\t**Error: Version < 3 no longer supported.\n\n"; 
					std::cout << "\t\tSee possible fixes at https://github.com/physicell-training/PhysiCell_common_errors\n\n"; 
					exit(-1); 

					done = true; 
				}

				if(version >= 2.0 - 1e-10 && version < 3.0 - 1e-10 && protocol == "CBHG" && done == false )
				{
					std::cout << "\tFormat: CSV (preprint version " << version << ")" << std::endl; 

					// parse_csv_rules_v2( input_filename ); // parse all rules in a CSV file 

					PhysiCell_settings.rules_enabled = true; 

					std::cout << "\t\t**Error: Version < 3 no longer supported.\n\n"; 
					std::cout << "\t\tSee possible fixes at https://github.com/physicell-training/PhysiCell_common_errors\n\n"; 
					exit(-1); 

					done = true; 
				}

				if(version >= 3.0 - 1e-10 && protocol == "CBHG" && done == false )
				{
					std::cout << "\tFormat: CSV (current version " << version << ")" << std::endl; 

					parse_csv_rules_v3( input_filename ); // parse all rules in a CSV file 

					PhysiCell_settings.rules_enabled = true; 

					done = true; 
				}

			}  


			if( done == false )
			{ std::cout << "\tWarning: Ruleset had unknown format (" << format << "). Skipping!" << std::endl; }
			else
			{ copy_file_to_output( input_filename ); }

		}
		else
		{ std::cout << "\tRuleset disabled ... " << std::endl; }
		node = node.next_sibling( "ruleset"); 		
	}
	return; 
}

void parse_rules_from_parameters_v0( void )
{
	bool enabled = parameters.bools( "rules_enabled" ); 

	// enabled? 
	if( enabled == false )
	{ return; }

	// get filename 

	std::string folder = parameters.strings( "rules_folder" ); 
	std::string filename = parameters.strings( "rules_filename" ); 
	std::string input_filename = folder + "/" + filename; 

	std::string filetype = parameters.strings( "rules_type" );  ; 

	// what kind? 
	if( filetype == "csv" || filetype == "CSV" )
	{
		std::cout << "Loading rules from CSV file " << input_filename << " ... " << std::endl; 
		// load_cells_csv( input_filename );
		parse_csv_rules_v0( input_filename ); 
		return; 
	}

	return; 
}

void spaces_to_underscore( std::string& str )
{
	for( int n=0 ; n < str.size(); n++ )
	{ if( str[n] == ' ' ){ str[n] = '_'; } }
}

/*
submit this as bugfix in PhysiCell (PhysiCell_settings.cpp)

template <class T>
int Parameters<T>::find_index( std::string search_name )
{
	auto search = name_to_index_map.find( search_name );
	if( search != name_to_index_map.end() )
	{ return search->second; }
	return -1; 
	// return name_to_index_map[ search_name ]; 
}

*/

void stream_annotated_English_rules( std::ostream& os )
{
	os << "Cell Hypothesis Rules" << std::endl << std::endl; 
	for( int n=0; n < cell_definitions_by_index.size(); n++ )
	{
		Cell_Definition* pCD = cell_definitions_by_index[n]; 
		Hypothesis_Ruleset* pHRS = find_ruleset( pCD ); 
		os << "In " << pHRS->cell_type << " cells:" << std::endl; 

		for( int k=0 ; k < pHRS->rules.size(); k++ ) 
		{ pHRS->rules[k]->English_display(os); }
		os << std::endl; 
	}
	return; 
} 

void stream_annotated_English_rules_HTML( std::ostream& os )
{
	os << "<html>\n<body><h1>Cell Hypothesis Rules</h1>" << std::endl; 
	for( int n=0; n < cell_definitions_by_index.size(); n++ )
	{
		os << "<p>"; 
		Cell_Definition* pCD = cell_definitions_by_index[n]; 
		Hypothesis_Ruleset* pHRS = find_ruleset( pCD ); 
		os << "In " << pHRS->cell_type << " cells:" << std::endl; 
		os << "<ul>" << std::endl; 
		for( int k=0 ; k < pHRS->rules.size(); k++ ) 
		{ pHRS->rules[k]->English_display_HTML(os); }
		os << "</ul>\n</p>" << std::endl; 
	}
	os << "</body>\n</html>" << std::endl; 
	return; 
	
} 

void save_annotated_English_rules( void )
{
	std::string filename = PhysiCell_settings.folder + "/rules.txt";
	std::ofstream of( filename , std::ios::out );
	stream_annotated_English_rules( of ); 
	of.close(); 
}

void stream_annotated_detailed_English_rules( std::ostream& os )
{
	os << "Cell Hypothesis Rules (detailed)" << std::endl << std::endl; 
	for( int n=0; n < cell_definitions_by_index.size(); n++ )
	{
		Cell_Definition* pCD = cell_definitions_by_index[n]; 
		Hypothesis_Ruleset* pHRS = find_ruleset( pCD ); 
		os << "In " << pHRS->cell_type << " cells:" << std::endl; 
		for( int k=0 ; k < pHRS->rules.size(); k++ ) 
		{ pHRS->rules[k]->English_detailed_display(os); }
		os << std::endl; 
	}
	return; 
} 

void stream_annotated_detailed_English_rules_HTML( std::ostream& os )
{
	os << "<html>\n<body><h1>Cell Hypothesis Rules (detailed)</h1>" << std::endl; 
	for( int n=0; n < cell_definitions_by_index.size(); n++ )
	{
		os << "<p>"; 
		Cell_Definition* pCD = cell_definitions_by_index[n]; 
		Hypothesis_Ruleset* pHRS = find_ruleset( pCD ); 
		os << "In " << pHRS->cell_type << " cells:" << std::endl; 
		os << "<ul>" << std::endl; 
		for( int k=0 ; k < pHRS->rules.size(); k++ ) 
		{ pHRS->rules[k]->English_detailed_display_HTML(os); }
		os << "</ul>\n</p>" << std::endl; 
	}
	os << "</body>\n</html>" << std::endl; 
	return; 
} 

void save_annotated_detailed_English_rules( void )
{
	std::string filename = PhysiCell_settings.folder + "/detailed_rules.txt";
	std::ofstream of( filename , std::ios::out );
	stream_annotated_detailed_English_rules( of ); 
	of.close(); 
}

void save_annotated_detailed_English_rules_HTML( void )
{
	std::string filename = PhysiCell_settings.folder + "/detailed_rules.html";
	std::ofstream of( filename , std::ios::out );
	stream_annotated_detailed_English_rules_HTML( of ); 
	of.close(); 
}

void save_annotated_English_rules_HTML( void )
{
	std::string filename = PhysiCell_settings.folder + "/rules.html";
	std::ofstream of( filename , std::ios::out );
	stream_annotated_English_rules_HTML( of ); 
	of.close(); 
}

// v0 version 

void export_rules_csv_v0( std::string filename )
{
	std::fstream fs( filename, std::ios::out );
	if( !fs )
	{
		std::cout << "Warning: Rules export file " << filename << " failed to open." << std::endl; 
		return; 
	}

	std::cout << "Exporting rules to file " << filename << " (v0 format) ... " << std::endl; 

	for( int n=0; n < cell_definitions_by_index.size(); n++ )
	{
		Cell_Definition* pCD = cell_definitions_by_index[n]; 
		Hypothesis_Ruleset* pHRS = find_ruleset( pCD ); 

		std::string cell_type = pHRS->cell_type; 

		std::cout << cell_type << " :: "; 
		for( int k=0 ; k < pHRS->rules.size(); k++ ) 
		{
			std::string behavior = pHRS->rules[k]->get_behavior()->get_name(); 
			std::cout << behavior << " : "; 

			double min_value = pHRS->rules[k]->min_value; 
			double max_value = pHRS->rules[k]->max_value; 
			double base_value = pHRS->rules[k]->base_value; 
			for( int i=0; i < pHRS->rules[k]->down_signals.size() ; i++ )
			{
				std::string signal = pHRS->rules[k]->down_signals[i]->get_name(); 
				std::cout << signal << " decreases "; 
				double half_max = pHRS->rules[k]->down_half_maxes[i];
				double hill_power = pHRS->rules[k]->down_hill_powers[i];
				bool use_for_dead = false; 

				// output the rule 
				fs << cell_type << "," << behavior << "," 
					<< min_value << "," << base_value << "," << max_value << ","
					<< signal << ",decreases," 
					<< half_max << "," << hill_power << "," 
					<< use_for_dead << std::endl; 
			}
			for( int i=0; i < pHRS->rules[k]->up_signals.size() ; i++ )
			{
				std::string signal = pHRS->rules[k]->up_signals[i]->get_name(); 
				std::cout << signal << " increases "; 
				double half_max = pHRS->rules[k]->up_half_maxes[i];
				double hill_power = pHRS->rules[k]->up_hill_powers[i];
				bool use_for_dead = false; 

				// output the rule 
				fs << cell_type << "," << behavior << "," 
					<< min_value << "," << base_value << "," << max_value << ","
					<< signal << ",increases," 
					<< half_max << "," << hill_power << "," 
					<< use_for_dead << std::endl; 
			}
			std::cout << std::endl; 
		}
	}

/*
0              1          2        3           4          5       6          7         8           9 
 Cell type, behavior, min value, base value, max value,   signal, direction, half-max, Hill power, dead
*/ 
	fs.close(); 

	std::cout << "Done!" << std::endl << std::endl; 

	return; 
}

// v1 

void export_rules_csv_v1( std::string filename )
{
	std::fstream fs( filename, std::ios::out );
	if( !fs )
	{
		std::cout << "Warning: Rules export file " << filename << " failed to open." << std::endl; 
		return; 
	}

	std::cout << "Exporting rules to file " << filename << " (v1 format) ... " << std::endl; 

	for( int n=0; n < cell_definitions_by_index.size(); n++ )
	{
		Cell_Definition* pCD = cell_definitions_by_index[n]; 
		Hypothesis_Ruleset* pHRS = find_ruleset( pCD ); 

		std::string cell_type = pHRS->cell_type; 

		for( int k=0 ; k < pHRS->rules.size(); k++ ) 
		{
			std::string behavior = pHRS->rules[k]->get_behavior()->get_name(); 

			double min_value = pHRS->rules[k]->min_value; 
			double max_value = pHRS->rules[k]->max_value; 
			double base_value = pHRS->rules[k]->base_value; 
			for( int i=0; i < pHRS->rules[k]->down_signals.size() ; i++ )
			{
				std::string signal = pHRS->rules[k]->down_signals[i]->get_name(); 
				std::string response = "decreases";
				double max_response = min_value;
				double half_max = pHRS->rules[k]->down_half_maxes[i];
				double hill_power = pHRS->rules[k]->down_hill_powers[i];
				bool use_for_dead = pHRS->rules[k]->down_applies_to_dead_cells[i];

				// output the rule 
				fs << cell_type << "," << signal << "," << response << "," << behavior << "," // 0,1,2,3
					<< base_value << "," << max_response << "," << half_max << "," << hill_power << "," // 4,5, 6,7, 
					<< use_for_dead << std::endl; // 8 
			}
			for( int i=0; i < pHRS->rules[k]->up_signals.size() ; i++ )
			{
				std::string signal = pHRS->rules[k]->up_signals[i]->get_name(); 
				std::string response = "increases";
				double max_response = max_value;
				double half_max = pHRS->rules[k]->up_half_maxes[i];
				double hill_power = pHRS->rules[k]->up_hill_powers[i];
				bool use_for_dead = pHRS->rules[k]->up_applies_to_dead_cells[i];

				// output the rule 
				fs << cell_type << "," << signal << "," << response << "," << behavior << "," // 0,1,2,3
					<< base_value << "," << max_response << "," << half_max << "," << hill_power << "," // 4,5, 6,7, 
					<< use_for_dead << std::endl; // 8 
			}
		}
	}

/*
 Cell type, signal, direcxtion, behavior, base, max_response, half-max, hill , dead 
*/ 
	fs.close(); 

	std::cout << "Done!" << std::endl << std::endl; 

	return; 
}

void export_rules_csv_v3( std::string filename )
{
	std::fstream fs( filename, std::ios::out );
	if( !fs )
	{
		std::cout << "Warning: Rules export file " << filename << " failed to open." << std::endl; 
		return; 
	}

	std::cout << "Exporting rules to file " << filename << " (v3 format) ... " << std::endl; 

	for( int n=0; n < cell_definitions_by_index.size(); n++ )
	{
		Cell_Definition* pCD = cell_definitions_by_index[n]; 
		Hypothesis_Ruleset* pHRS = find_ruleset( pCD ); 

		std::string cell_type = pHRS->cell_type; 

		for( int k=0 ; k < pHRS->rules.size(); k++ ) 
		{
			std::string behavior = pHRS->rules[k]->get_behavior()->get_name(); 

			double min_value = pHRS->rules[k]->min_value; 
			double max_value = pHRS->rules[k]->max_value; 
			double base_value = pHRS->rules[k]->base_value; 
			for( int i=0; i < pHRS->rules[k]->down_signals.size() ; i++ )
			{
				std::string signal = pHRS->rules[k]->down_signals[i]->get_name(); 
				std::string response = "decreases";
				double max_response = min_value; 
				double half_max = pHRS->rules[k]->down_half_maxes[i];
				double hill_power = pHRS->rules[k]->down_hill_powers[i];
				bool use_for_dead = pHRS->rules[k]->down_applies_to_dead_cells[i];

				// output the rule 
				fs << cell_type << "," << signal << "," << response << "," << behavior << "," // 0,1,2,3
					// << base_value << "," 
					<< max_response << "," << half_max << "," << hill_power << "," // 4,5, 6,7, 
					<< use_for_dead << std::endl; // 8 
			}
			for( int i=0; i < pHRS->rules[k]->up_signals.size() ; i++ )
			{
				std::string signal = pHRS->rules[k]->up_signals[i]->get_name(); 
				std::string response = "increases";
				double max_response = max_value;
				double half_max = pHRS->rules[k]->up_half_maxes[i];
				double hill_power = pHRS->rules[k]->up_hill_powers[i];
				bool use_for_dead = pHRS->rules[k]->up_applies_to_dead_cells[i];

				// output the rule 
				fs << cell_type << "," << signal << "," << response << "," << behavior << "," // 0,1,2,3
					// << base_value << "," 
					<< max_response << "," << half_max << "," << hill_power << "," // 4,5,6,
					<< use_for_dead << std::endl; // 7
			}
		}
	}

/*
 Cell type, signal, direcxtion, behavior, base, max_response, half-max, hill , dead 
*/ 
	fs.close(); 

	std::cout << "Done!" << std::endl << std::endl; 

	return; 
}


std::vector<double> UniformInUnitDisc( void )
{
	static double two_pi = 6.283185307179586; 
	double theta = UniformRandom(); // U(0,1)
	theta *= two_pi; // U(0,2*pi)
	double r = sqrt( UniformRandom() );  // sqrt( U(0,1) )
	return { r*cos(theta), r*sin(theta), 0.0 }; 
}

std::vector<double> UniformInUnitSphere( void )
{
	// reference: https://doi.org/10.1063/1.168311, adapting equation 13

	static double two_pi = 6.283185307179586; 

    double T = UniformRandom(); 
	double sqrt_T = sqrt(T); 
	double sqrt_one_minus_T = 1.0;
	sqrt_one_minus_T -= T; 
	sqrt_one_minus_T = sqrt( sqrt_one_minus_T ); 

	double param1 = pow( UniformRandom() , 0.33333333333333333333333333333333333333 );  //  xi^(1/3), 
    double param2 = param1; // xi^(1/3)
	param2 *= 2.0; // 2 * xi^(1/3)
	param2 *= sqrt_T; // 2 * xi(1) * T^(1/2)
	param2 *= sqrt_one_minus_T; //  2 * xi(1) * T^(1/2) * (1-T)^(1/2)
	
    double theta = UniformRandom(); // U(0,1)
	theta *= two_pi; // U(0,2*pi)
	
	return { param2*sin(theta) , param2*cos(theta), param1*(1-2*T) }; 
}

std::vector<double> UniformInAnnulus( double r1, double r2 )
{
	static double two_pi = 6.283185307179586; 

    double theta = UniformRandom(); 
	theta *= two_pi; 
	double r1_2 = r1*r1; 
	double r2_2 = r2*r2; 

    double r = sqrt( r1_2 + (r2_2-r1_2) * UniformRandom() ); 
    double x = r*cos(theta); 
    double y = r*sin(theta); 
    return {x,y,0.0}; 
}

std::vector<double> UniformInShell( double r1, double r2 )
{
	static double two_pi = 6.283185307179586; 

    double T = UniformRandom(); 
	double sqrt_T = sqrt(T); 
	double sqrt_one_minus_T = 1.0;
	sqrt_one_minus_T -= T; 
	sqrt_one_minus_T = sqrt( sqrt_one_minus_T ); 

	double param1 = pow( UniformRandom() , 0.33333333333333333333333333333333333333 );  //  xi^(1/3), 
	// param1 *= (r2-r1); 
	// param1 += r1; 
    double param2 = param1; // xi^(1/3)
	param2 *= 2.0; // 2 * xi^(1/3)
	param2 *= sqrt_T; // 2 * xi(1) * T^(1/2)
	param2 *= sqrt_one_minus_T; //  2 * xi(1) * T^(1/2) * (1-T)^(1/2)
	
    double theta = UniformRandom(); // U(0,1)
	theta *= two_pi; // U(0,2*pi)
	
	return { param2*sin(theta) , param2*cos(theta), param1*(1-2*T) }; 
}

void setup_cell_rules( void )
{
	// setup 
	intialize_hypothesis_rulesets(); 

	// load rules 
	parse_rules_from_pugixml(); 

	// display rules to screen
	display_hypothesis_rulesets( std::cout );

	// save annotations 
	save_annotated_detailed_English_rules(); 
	save_annotated_detailed_English_rules_HTML(); 
	save_annotated_English_rules(); 
	save_annotated_English_rules_HTML(); 

	// save dictionaries 
	std::string dictionary_file = "./" + PhysiCell_settings.folder + "/dictionaries.txt";
	std::ofstream dict_of( dictionary_file , std::ios::out ); 

	// display_signal_dictionary( dict_of ); // done 
	display_signal_dictionary_with_synonyms( dict_of ); // 
	// display_behavior_dictionary( dict_of ); // done 
	display_behavior_dictionary_with_synonyms( dict_of ); // done 
	dict_of.close(); 

	// save rules (v3)
	std::string rules_file = PhysiCell_settings.folder + "/cell_rules_parsed.csv"; 
	export_rules_csv_v3( rules_file ); 

	return; 
}


};
