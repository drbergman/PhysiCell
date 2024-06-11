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
# Copyright (c) 2015-2023, Paul Macklin and the PhysiCell Project             #
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


#ifndef __PhysiCell_signal_response__
#define __PhysiCell_signal_response__

#include "../core/PhysiCell.h"
// #include "./PhysiCell_constants.h" 
// #include "./PhysiCell_phenotype.h" 
// #include "./PhysiCell_cell.h" 

namespace PhysiCell{
double sum_fn( std::vector<double> signals_in );
double diff_fn( std::vector<double> signals_in );

class AbstractSignal
{
public:
    virtual double evaluate( Cell* pCell ) = 0;

    virtual ~AbstractSignal() {}
};

class AbstractSignalAggregator : public AbstractSignal
{
public:
    virtual std::vector<double> get_signals( Cell* pCell ) = 0;
    std::function<double(std::vector<double>)> aggregator;

    double evaluate( Cell* pCell ) {
        return aggregator(get_signals(pCell));
    }
};

class SignalAggregator : public AbstractSignalAggregator
{
private:
    std::vector<AbstractSignal *> signals;

public:
    std::vector<double> get_signals( Cell* pCell ) {
        std::vector<double> signal_values;
        for (auto& signal : signals) {
            signal_values.push_back(signal->evaluate(pCell));
        }
        return signal_values;
    }

    double aggregate_signals( Cell* pCell ) {
        std::vector<double> signal_values = get_signals(pCell);
        return aggregator(signal_values);
    }

    void add_signal(AbstractSignal *pSignal) {
        signals.push_back(pSignal);
    }

    SignalAggregator() {
        aggregator = sum_fn;
    }
};

class SignalMediator : public AbstractSignalAggregator
{
private:
    AbstractSignal *decreasing_signal;
    AbstractSignal *increasing_signal;

public:
    // values that the aggregator can use to calculate the output
    double min_value = 0.1;
    double base_value = 1;
    double max_value = 10;

    std::vector<double> get_signals( Cell* pCell ) {
        return {increasing_signal->evaluate(pCell), decreasing_signal->evaluate(pCell)};
    }

    double default_mediator_aggregator(std::vector<double> signals_in)
    {
    std::cout << "Default mediator aggregator" << std::endl;
    std::cout << "\tSignal 1: " << signals_in[0] << ", Signal 2: " << signals_in[1] << ", Min: " << min_value << ", Base: " << base_value << ", Max: " << max_value << std::endl;
    double out = base_value;
    std::cout << "\tBase value: " << out << std::endl;
    out += (max_value - base_value) * signals_in[1];
    std::cout << "\tAfter adding signal 2: " << out << std::endl;
    out *= 1-signals_in[0];
    std::cout << "\tAfter multiplying by signal 1: " << out << std::endl;
    out += min_value * signals_in[0];
    std::cout << "\tAfter adding signal 1: " << out << std::endl;
    return out;
    // std::cout << "Mediating signals" << std::endl;
    // std::cout << "\tSignal 1: " << signals_in[0] << std::endl;
    // std::cout << "\tSignal 2: " << signals_in[1] << std::endl << std::endl;
    // return min_value * signals_in[0] + (base_value + (max_value - base_value) * signals_in[1]) * (1 - signals_in[0]);
    };

    double aggregate_signals( Cell* pCell ) {
        std::vector<double> signal_values = get_signals(pCell);
        return aggregator(signal_values);
    }

    SignalMediator() {
        aggregator = [this](std::vector<double> signals_in) {
            return this->default_mediator_aggregator(signals_in);
        };
    }
    SignalMediator(AbstractSignal *pIncreasingSignal, AbstractSignal *pDecreasingSignal)
    {
        if (pIncreasingSignal == nullptr || pDecreasingSignal == nullptr)
        {
            throw std::invalid_argument("Null pointer passed to SignalMediator constructor");
        }
        std::cout << "Creating mediator" << std::endl;
        increasing_signal = pIncreasingSignal;
        decreasing_signal = pDecreasingSignal;
        aggregator = [this](std::vector<double> signals_in)
        { return this->default_mediator_aggregator(signals_in); };
    }
    SignalMediator(AbstractSignal *pIncreasingSignal, AbstractSignal *pDecreasingSignal, double min, double base, double max)
    {
        if (pIncreasingSignal == nullptr || pDecreasingSignal == nullptr)
        {
            throw std::invalid_argument("Null pointer passed to SignalMediator constructor");
        }
        std::cout << "Creating mediator" << std::endl;
        increasing_signal = pIncreasingSignal;
        decreasing_signal = pDecreasingSignal;
        min_value = min;
        base_value = base;
        max_value = max;
        aggregator = [this](std::vector<double> signals_in)
        { return this->default_mediator_aggregator(signals_in); };
    }
};

class BehaviorRule
{
public:
    double min_value;
    double base_value;
    double max_value;

    AbstractSignal* signal;

    void apply(Cell* pCell) {
        double return_value = signal->evaluate(pCell);
        std::cout << "Cell behavior set to " << return_value << std::endl;
        return;
        // double signal_value = signal->evaluate(pCell);
        // double return_value;
        // if (signal_value < 0) {
        //     return_value = min_value * signal_value + base_value * (1 - signal_value);
        // }
        // else if (signal_value > 0) {
        //     return_value = max_value * signal_value + base_value * (1 - signal_value);
        // }
        // else {
        //     return_value = base_value;
        // }
        // std::cout << "Cell behavior set to " << return_value << std::endl;
    }
};

// scales for the signals 
extern std::vector<double> signal_scales; 
// easy access to get or set scales 
double& signal_scale( std::string signal_name ); // done 
double& signal_scale( int signal_index ); // done 

// create the signal and behavior dictionaries 
void setup_signal_behavior_dictionaries( void ); // done 

// display dictionaries 
void display_signal_dictionary( void ); // done 
void display_behavior_dictionary( void ); // done 

void display_signal_dictionary( std::ostream& os ); // done 
void display_behavior_dictionary( std::ostream& os ); // done 

void display_signal_dictionary_with_synonyms( void ); // done 
void display_behavior_dictionary_with_synonyms( void ); // done 

/* signal functions */ 

// find index for named signal (returns -1 if not found)
int find_signal_index( std::string signal_name ); // done 

// coming soon: 
std::vector<int> find_signal_indices( std::vector<std::string> signal_names ); // done 

// get the name of a signal index 
std::string signal_name( int i ); // done 

// create a full signal vector 
std::vector<double> get_signals( Cell* pCell ); // done 

// create a signal vector of only the cell contacts 
std::vector<double> get_cell_contact_signals( Cell* pCell ); // done 

// create a subset of the signal vector with the supplied indicies 
std::vector<double> get_selected_signals( Cell* pCell , std::vector<int> indices ); // done 
std::vector<double> get_selected_signals( Cell* pCell , std::vector<std::string> names );  // done 

// grab a single signal by its index or name 
double get_single_signal( Cell* pCell, int index ); // done 
double get_single_signal( Cell* pCell, std::string name ); // done 

/* behavior functions */ 

// find index for named behavior / response / parameter (returns -1 if not found)
int find_parameter_index( std::string response_name ); // done
int find_behavior_index( std::string response_name ); // done 

std::vector<int> find_behavior_indices( std::vector<std::string> behavior_names ); // done 

// get the name of a behavior index 
std::string behavior_name( int i ); // done 

// make a properly sized behavior vector 
std::vector<double> create_empty_behavior_vector(); // done 

// write a full behavior vector (phenotype parameters) to the cell 
void set_behaviors( Cell* pCell , std::vector<double> parameters ); // done 

// write a selected set of behavior parameters to the cell 
void set_selected_behaviors( Cell* pCell , std::vector<int> indices , std::vector<double> parameters ); // done 
void set_selected_behaviors( Cell* pCell , std::vector<std::string> names , std::vector<double> parameters ); // done 

// write a single behavior parameter 
void set_single_behavior( Cell* pCell, int index , double parameter ); // done  
void set_single_behavior( Cell* pCell, std::string name , double parameter ); // done 

/* get current behaviors */ 

// get all current behavior
std::vector<double> get_behaviors( Cell* pCell ); // done 

// get selected current behavior
std::vector<double> get_behaviors( Cell* pCell , std::vector<int> indices ); // doen 
std::vector<double> get_behaviors( Cell* pCell , std::vector<std::string> names ); // done 

// get single current behavior 
double get_single_behavior( Cell* pCell , int index ); // done 
double get_single_behavior( Cell* pCell , std::string name ); // done 

/* get base behaviors (from cell definition) */ 

// get all base behaviors (from cell's definition) 
std::vector<double> get_base_behaviors( Cell* pCell );  // done 

// get selected base behaviors (from cell's definition)
std::vector<double> get_base_behaviors( Cell* pCell , std::vector<int> indices ); // done 
std::vector<double> get_base_behaviors( Cell* pCell , std::vector<std::string> names ); // done 

// get single base behavior (from cell's definition)
double get_single_base_behavior( Cell* pCell , int index ); // done 
double get_single_base_behavior( Cell* pCell , std::string name ); // done 

double get_single_base_behavior( Cell_Definition* pCD , std::string name ); 


class ElementarySignal : public AbstractSignal
{
public:
    std::string signal;
    bool applies_to_dead_cells;

    virtual double transformer( double signal ) {
        return signal;
    }

    double evaluate( Cell* pCell );

    virtual ~ElementarySignal() {}
};

class HillSignal : public ElementarySignal
{
public:
    double half_max;
    double hill_power;

    double transformer(double signal)
    {
        std::cout << "Evaluating Hill signal" << std::endl;
        signal /= half_max;
        signal = pow(signal, hill_power);
        std::cout << "\tReturning " << signal / (1 + signal) << "\n" << std::endl;
        return signal / (1 + signal);
    }
};

class LinearSignal : public ElementarySignal
{
public:
    double signal_min;
    double signal_max;

    double transformer(double signal)
    {
        std::cout << "Evaluating linear signal" << std::endl;
        if (signal < signal_min) {
            std::cout << "\tReturning 0\n" << std::endl;
            return 0;
        }
        else if (signal > signal_max) {
            std::cout << "\tReturning 1\n" << std::endl;
            return 1;
        }
        else {
            std::cout << "\tReturning " << (signal - signal_min) / (signal_max - signal_min) << "\n" << std::endl;
            return (signal - signal_min) / (signal_max - signal_min);
        }
    }
};

class HeavisideSignal : public ElementarySignal
{
public:
    double threshold;

    double transformer(double signal)
    {
        if (signal < threshold) {
            return 0;
        }
        else {
            return 1;
        }
    }
};


}; 

#endif 
