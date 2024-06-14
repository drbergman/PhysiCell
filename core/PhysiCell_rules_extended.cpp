#include "PhysiCell_rules_extended.h"
#include "PhysiCell_signal_behavior.h"

using namespace BioFVM; 

namespace PhysiCell{
double ElementarySignal::evaluate(Cell *pCell)
{
	if (!applies_to_dead_cells && pCell->phenotype.death.dead)
	{
		return 0;
	}
	return transformer(get_single_signal(pCell, signal_name));
}

BehaviorRule::BehaviorRule(std::string cell_type, std::string behavior_name, double min_behavior, double max_behavior)
{
    pCell_Definition = find_cell_definition(cell_type);
    behavior = behavior_name;
    min_value = min_behavior;
    max_value = max_behavior;
    base_value = get_single_base_behavior(pCell_Definition, behavior_name);
}

void BehaviorRule::add_signal(std::string, std::string response)
{
	PartialHillSignal *pDecreasingSignal = new PartialHillSignal();
	PartialHillSignal *pIncreasingSignal = new PartialHillSignal();
	signal = new MediatorSignal(pDecreasingSignal, pIncreasingSignal);
	return;
}

void BehaviorRule::apply(Cell *pCell)
{
	double param = evaluate(pCell);
	set_single_behavior(pCell, behavior, param);
	return;
}

void BehaviorRule::sync_to_cell_definition( Cell_Definition* pCD )
{
	if( pCD == NULL )
	{ return; }

	cell_type = pCD->name; 
	pCell_Definition = pCD; 

	// sync base behavior 
	base_value = get_single_base_behavior(pCD,behavior); 

	return; 
}

void BehaviorRuleset::add_behavior(std::string behavior, double min_behavior, double max_behavior)
{
    // check: is this a valid signal? (is it in the dictionary?)
    if (find_behavior_index(behavior) < 0)
    {
        std::cout << "Warning! Attempted to add behavior " << behavior << " which is not in the dictionary." << std::endl;
        std::cout << "Either fix your model or add the missing behavior to the simulation." << std::endl;

        exit(-1);
    }

    // first, check. Is there already a ruleset?
    for (auto &rule : rules)
    {
        if (rule->behavior == behavior)
        {
            std::cout << "ERROR: Attempted to add behavior " << behavior << " which is already in the ruleset." << std::endl;
            std::cout << "\tIf you want to change the min and max values, use the set_min_max_behavior function." << std::endl;
            exit(-1);
        }
    }
    // if not, add it
    BehaviorRule *pBR = new BehaviorRule(cell_type, behavior, min_behavior, max_behavior);
}

void BehaviorRuleset::sync_to_cell_definition(std::string cell_type)
{
    pCell_Definition = find_cell_definition(cell_type);
    sync_to_cell_definition_finish(cell_type);
}

// BehaviorRule* Hypothesis_Ruleset::add_behavior( std::string behavior , double min_behavior, double max_behavior )
// {
//     // check: is this a valid signal? (is it in the dictionary?)
//     if( find_behavior_index(behavior) < 0 )
//     {
//         std::cout << "Warning! Attempted to add behavior " << behavior << " which is not in the dictionary." << std::endl; 
//         std::cout << "Either fix your model or add the missing behavior to the simulation." << std::endl; 

//         exit(-1); 
//     }

// 	// first, check. Is there already a ruleset? 
// 	auto search = rules_map.find( behavior ); 

// 		// if not, add it 
// 	if( search == rules_map.end() )
// 	{
// 		BehaviorRule *pBR = new BehaviorRule;

// 		pBR->behavior = behavior; 

// 		pBR->sync_to_cell_definition( pCell_Definition ); 

// 		pBR->min_value = min_behavior; 
// 		pBR->max_value = max_behavior;

// 		rules.push_back(pBR);
// 		rules_map[ behavior ] = pBR; 

// 		return pBR; 
// 	}

// 		// otherwise, edit it 
// 	BehaviorRule* pBR = search->second; 

// 	/*
// 		// March 28 2023 fix  : let's not overwrite eixsting values
// 	pBR->min_value = min_behavior; 
// 	pBR->max_value = max_behavior; 
// 	*/

// 	return pBR; 
// }

// BehaviorRule* Hypothesis_Ruleset::add_behavior( std::string behavior )
// { 
// 	double min_behavior = 9e99; // Min behaviour high value
// 	double max_behavior = -9e99; // Max behaviour low value
// 	return Hypothesis_Ruleset::add_behavior( behavior, min_behavior, max_behavior );
// }

// BehaviorRule* Hypothesis_Ruleset::find_behavior( std::string name )
// {
//     auto search = rules_map.find( name); 
// 	if( search == rules_map.end() )
// 	{
// 		// std::cout << "Warning! Ruleset does not contain " << name << std::endl; 
// 		// std::cout << "         Returning NULL." << std::endl; 
// 		return NULL; 
// 	}

// 	return search->second; 
// }

// BehaviorRule& Hypothesis_Ruleset::operator[]( std::string name )
// {
// 	BehaviorRule* pBR = find_behavior(name);
// 	return *pBR; 
// } 

void parse_xml_rules_extended(std::string filename)
{
	bool physicell_rules_dom_initialized = false;
	pugi::xml_document physicell_rules_doc;
	pugi::xml_node physicell_rules_root;
	std::cout << "Using rules file " << filename << " ... " << std::endl;
	pugi::xml_parse_result result = physicell_rules_doc.load_file(filename.c_str());

	if (result.status != pugi::xml_parse_status::status_ok)
	{
		std::cout << "Error loading " << filename << "!" << std::endl;
		exit(-1);
	}

	physicell_rules_root = physicell_rules_doc.child("hypothesis_rulesets");
	physicell_rules_dom_initialized = true;

	pugi::xml_node ruleset_node;
	ruleset_node = xml_find_node(physicell_rules_root, "hypothesis_ruleset");
	std::vector<std::string> cell_definitions_ruled;

	while (ruleset_node)
	{
		std::string cell_type = ruleset_node.attribute("name").value();
		// check if cell_type has already had a ruleset defined for it
		for (int i = 0; i < cell_definitions_ruled.size(); i++)
		{
			if (cell_type == cell_definitions_ruled[i])
			{
				std::cout << "XML Rules ERROR: Two rulesets for " << cell_type << " found." << std::endl
						  << "\tCombine them into a single ruleset please :)" << std::endl
						  << "\tSupport for rules across multiple files for the same cell type not yet supported." << std::endl;
				exit(-1);
			}
		}
		cell_definitions_ruled.push_back(cell_type);

		std::vector<std::string> behaviors_ruled;

		pugi::xml_node behavior_node = ruleset_node.child("behavior");
		while (behavior_node)
		{
			std::string behavior = behavior_node.attribute("name").value();
			for (int i = 0; i < behaviors_ruled.size(); i++)
			{
				if (behavior == behaviors_ruled[i])
				{
					std::cout << "XML Rules ERROR: The behavior " << behavior << " is being set again for " << cell_type << "." << std::endl
							  << "\tSelect only one to keep. :)" << std::endl;
					exit(-1);
				}
			}
			behaviors_ruled.push_back(behavior);
			AbstractSignal *pAS = parse_abstract_signal(behavior_node);
			behavior_node = behavior_node.next_sibling("behavior");
		}
		ruleset_node = ruleset_node.next_sibling("hypothesis_ruleset");
	}
	return;
}

AbstractSignal *parse_abstract_signal(const pugi::xml_node node)
{
	// idea: also start with mediator. if only one of increasing or decreasing provided, then can scrap the mediator and use either an elementary signal or an aggregator
	AbstractSignal* pAS;
	if (signal_is_mediator(node))
	{
		pAS = parse_mediator_signal(node);
	}
	else if (signal_is_aggregator(node))
	{
		pAS = parse_aggregator_signal(node);
	}
	else if (signal_is_elementary(node))
	{
		pAS = parse_elementary_signal(node);
	}
	else
	{
		std::cout << "ERROR: Signal type not recognized." << std::endl;
		exit(-1);
	}
	return pAS;
}

AbstractSignal *parse_mediator_signal(const pugi::xml_node mediator_node)
{
	pugi::xml_node decreasing_signals_node = mediator_node.child("decreasing_signals");
	pugi::xml_node increasing_signals_node = mediator_node.child("increasing_signals");

	AbstractSignal *pDecreasingSignal;
	AbstractSignal *pIncreasingSignal;
	if (decreasing_signals_node)
	{
		pDecreasingSignal = parse_aggregator_signal(decreasing_signals_node);
		// process_decreasing_signals(behavior_node, cell_type, behavior);
	}
	if (increasing_signals_node)
	{
		pIncreasingSignal = parse_aggregator_signal(increasing_signals_node);
		// process_increasing_signals(behavior_node, cell_type, behavior);
	}
	MediatorSignal *pSM = new MediatorSignal();
	return;
}

AbstractSignal *parse_aggregator_signal(const pugi::xml_node aggregator_node)
{
	pugi::xml_node signal_node = xml_find_node(aggregator_node, "signal");
	while(signal_node)
	{
		parse_abstract_signal(signal_node);
		signal_node = signal_node.next_sibling("signal");
	}
}

AbstractSignal *parse_elementary_signal(const pugi::xml_node elementary_node)
{
	// parse the signal
	// add the signal to the signal dictionary
	// return;
}

bool signal_is_mediator(pugi::xml_node node)
{
	bool is_mediator = node.name() == "behavior"; // by default, assume that the top-level signal is a mediator
	is_mediator |= (node.attribute("type")) && node.attribute("type").value() == "mediator";
	is_mediator |= xml_find_node(node, "decreasing_signals"); // if the type was not declared but there are decreasing signals, then it is a mediator
	is_mediator |= xml_find_node(node, "increasing_signals"); // if the type was not declared but there are increasing signals, then it is a mediator
	return is_mediator;
}

bool signal_is_aggregator(pugi::xml_node node)
{
	bool is_aggregator = (node.attribute("type")) && node.attribute("type").value() == "aggregator";
	pugi::xml_node signal_node = xml_find_node(node, "signal");
	if (signal_node && signal_node.next_sibling("signal"))
	{
		return true; // if there is more than one signal, then it is an aggregator
	}
	return is_aggregator;
}

bool signal_is_elementary(pugi::xml_node node)
{
	if (node.name() == "signal")
	{
		return true; // if the node is a signal node, then it is an elementary signal by process of elimination
	}
	return false;
}

};