#include "PhysiCell_rules_extended.h"

namespace PhysiCell{

#ifndef __PhysiCell_rules_extended_cpp__
#define __PhysiCell_rules_extended_cpp__
#endif

double sum_aggregator(std::vector<double> signals_in)
{
    double signal_sum = 0;
    for (auto &signal : signals_in)
    {
        signal_sum += signal;
    }
    return signal_sum;
}

double multivariate_hill_aggregator(std::vector<double> partial_hill_signals)
{
    double out = sum_aggregator(partial_hill_signals);
    return out / (1 + out);
}

double MediatorSignal::decreasing_dominant_mediator(std::vector<double> signals_in)
{
    return min_value * signals_in[0] + (base_value + (max_value - base_value) * signals_in[1]) * (1 - signals_in[0]);
}

double MediatorSignal::increasing_dominant_mediator(std::vector<double> signals_in)
{
    return max_value * signals_in[1] + (base_value + (min_value - base_value) * signals_in[0]) * (1 - signals_in[1]);
}

double MediatorSignal::neutral_mediator(std::vector<double> signals_in)
{
    return base_value + (min_value - base_value) * signals_in[0] + (max_value - base_value) * signals_in[1];
}

void MediatorSignal::use_increasing_dominant_mediator()
{
	aggregator = [this](std::vector<double> signals_in)
	{
		return this->increasing_dominant_mediator(signals_in);
	};
}

void MediatorSignal::use_neutral_mediator()
{
	aggregator = [this](std::vector<double> signals_in)
	{
		return this->neutral_mediator(signals_in);
	};
}

void ElementarySignal::add_reference(std::unique_ptr<SignalReference> pSR)
{
	signal_reference = std::move(pSR);	   // set the reference
	if (add_signal_reference)
	{
		add_signal_reference(signal_reference.get()); // depending on the class of ElementarySignal, do what else is needed when adding a reference
	}
}

double ElementarySignal::evaluate(Cell *pCell)
{
	if (!applies_to_dead_cells && pCell->phenotype.death.dead)
	{
		return 0;
	}
	double signal = get_single_signal(pCell, signal_name);
	if (signal_reference != nullptr)
	{
		signal = signal_reference->coordinate_transform(signal);
	}
	return transformer(signal);
}

void BehaviorRule::apply(Cell *pCell)
{
	double param = signal->evaluate(pCell);
	set_single_behavior(pCell, behavior, param);
	return;
}

void BehaviorRuleset::apply(Cell *pCell)
{
	for (auto &rule : rules)
	{
		rule->apply(pCell);
	}
	return;
}

std::unordered_map< Cell_Definition* , std::unique_ptr<BehaviorRuleset> > behavior_rulesets; 

void add_behavior_ruleset( Cell_Definition* pCD )
{
	auto search = behavior_rulesets.find( pCD );
	if (search == behavior_rulesets.end())
	{
		// I don't think I need to sync to the cell definition
        behavior_rulesets[pCD] = std::unique_ptr<BehaviorRuleset>(new BehaviorRuleset()); // cannot use make_unique since that is only introduced in C++14
	}
	return; 
}

void intialize_behavior_rulesets( void )
{
	behavior_rulesets.clear(); // empty(); 
	for( int n; n < cell_definitions_by_index.size() ; n++ )
	for (auto &pCD : cell_definitions_by_index)
	{
		add_behavior_ruleset(pCD); 
	}
	return; 
}

void setup_behavior_rules( void )
{
	intialize_behavior_rulesets();

	parse_behavior_rules_from_pugixml();

	// display_behavior_rules( std::cout );

	// save_annotated_detailed_English_behavior_rules(); 
	// save_annotated_detailed_English_behavior_rules_HTML(); 
	// save_annotated_English_behavior_rules(); 
	// save_annotated_English_behavior_rules_HTML(); 

	std::string dictionary_file = "./" + PhysiCell_settings.folder + "/dictionaries.txt";
	std::ofstream dict_of( dictionary_file , std::ios::out ); 

	display_signal_dictionary( dict_of ); // done 
	display_behavior_dictionary( dict_of ); // done 
	dict_of.close(); 

	// // save rules (v1)
	// std::string rules_file = PhysiCell_settings.folder + "/cell_rules.csv"; 
	// export_rules_csv_v1( rules_file ); 
	return; 
}

void apply_behavior_ruleset( Cell* pCell )
{
	Cell_Definition* pCD = find_cell_definition( pCell->type_name ); 
	behavior_rulesets[pCD]->apply( pCell );
	return; 
}

// BehaviorRule::BehaviorRule(std::string cell_type, std::string behavior_name, double min_behavior, double max_behavior)
// {
//     pCell_Definition = find_cell_definition(cell_type);
//     behavior = behavior_name;
//     min_value = min_behavior;
//     max_value = max_behavior;
//     base_value = get_single_base_behavior(pCell_Definition, behavior_name);
// }

// void BehaviorRule::append_signal(std::string, std::string response)
// {
// 	PartialHillSignal *pDecreasingSignal = new PartialHillSignal();
// 	PartialHillSignal *pIncreasingSignal = new PartialHillSignal();
// 	signal = new MediatorSignal(pDecreasingSignal, pIncreasingSignal);
// 	return;
// }

// void BehaviorRule::apply(Cell *pCell)
// {
// 	double param = evaluate(pCell);
// 	set_single_behavior(pCell, behavior, param);
// 	return;
// }

// void BehaviorRule::sync_to_cell_definition( Cell_Definition* pCD )
// {
// 	if( pCD == NULL )
// 	{ return; }

// 	cell_type = pCD->name; 
// 	pCell_Definition = pCD; 

// 	// sync base behavior 
// 	base_value = get_single_base_behavior(pCD,behavior); 

// 	return; 
// }

// void BehaviorRuleset::add_behavior(std::string behavior, double min_behavior, double max_behavior)
// {
//     // check: is this a valid signal? (is it in the dictionary?)
//     if (find_behavior_index(behavior) < 0)
//     {
//         std::cout << "Warning! Attempted to add behavior " << behavior << " which is not in the dictionary." << std::endl;
//         std::cout << "Either fix your model or add the missing behavior to the simulation." << std::endl;

//         exit(-1);
//     }

//     // first, check. Is there already a ruleset?
//     for (auto &rule : rules)
//     {
//         if (rule->behavior == behavior)
//         {
//             std::cout << "ERROR: Attempted to add behavior " << behavior << " which is already in the ruleset." << std::endl;
//             std::cout << "\tIf you want to change the min and max values, use the set_min_max_behavior function." << std::endl;
//             exit(-1);
//         }
//     }
//     // if not, add it
//     BehaviorRule *pBR = new BehaviorRule(cell_type, behavior, min_behavior, max_behavior);
// }

// void BehaviorRuleset::sync_to_cell_definition(std::string cell_type)
// {
//     pCell_Definition = find_cell_definition(cell_type);
//     sync_to_cell_definition_finish(cell_type);
// }

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

void parse_xml_behavior_rules(const std::string filename)
{
	bool physicell_rules_dom_initialized = false;
	pugi::xml_document physicell_rules_doc;
	pugi::xml_node physicell_rules_root;
	std::cout << "Using rules file " << filename << " ... " << std::endl;
	pugi::xml_parse_result result = physicell_rules_doc.load_file(filename.c_str());

	if (result.status != pugi::xml_parse_status::status_ok)
	{
		std::cerr << "Error loading " << filename << "!" << std::endl;
		exit(-1);
	}

	physicell_rules_root = physicell_rules_doc.child("behavior_rulesets");
	physicell_rules_dom_initialized = true;

	pugi::xml_node ruleset_node;
	ruleset_node = xml_find_node(physicell_rules_root, "behavior_ruleset");
	std::vector<std::string> cell_definitions_ruled;

	while (ruleset_node)
	{
		std::string cell_type = ruleset_node.attribute("name").value();
		// check if cell_type has already had a ruleset defined for it
		for (int i = 0; i < cell_definitions_ruled.size(); i++)
		{
			if (cell_type == cell_definitions_ruled[i])
			{
				std::cerr << "XML Rules ERROR: Two rulesets for " << cell_type << " found." << std::endl
						  << "\tCombine them into a single ruleset please :)" << std::endl
						  << "\tSupport for rules across multiple files for the same cell type not yet supported." << std::endl;
				exit(-1);
			}
		}
		cell_definitions_ruled.push_back(cell_type);

		std::vector<std::string> behaviors_ruled;
		Cell_Definition *pCD = find_cell_definition(cell_type);

		pugi::xml_node behavior_node = ruleset_node.child("behavior");
		while (behavior_node)
		{
			std::string behavior = behavior_node.attribute("name").value();
			for (int i = 0; i < behaviors_ruled.size(); i++)
			{
				if (behavior == behaviors_ruled[i])
				{
					std::cerr << "XML Rules ERROR: The behavior " << behavior << " is being set again for " << cell_type << "." << std::endl
							  << "\tSelect only one to keep. :)" << std::endl;
					exit(-1);
				}
			}
			behaviors_ruled.push_back(behavior);

			std::unique_ptr<AbstractSignal> pAS = parse_abstract_signal(behavior_node); // make this a unique ptr
			std::unique_ptr<BehaviorRule> pBR = std::unique_ptr<BehaviorRule>(new BehaviorRule(behavior, std::move(pAS))); // make sure this is initalized correctly (and probably using unique_ptrs and releasing the pAS pointer, too)
			behavior_rulesets[pCD]->add_behavior_rule(std::move(pBR)); // make sure add_behavior_rule is defined, unique ptrs are used correctly, and pointers are released if need be

			behavior_node = behavior_node.next_sibling("behavior");
		}
		ruleset_node = ruleset_node.next_sibling("behavior_ruleset");
	}
	return;
}

std::unique_ptr<AbstractSignal> parse_abstract_signal(pugi::xml_node node)
{
	// idea: also start with mediator. if only one of increasing or decreasing provided, then can scrap the mediator and use either an elementary signal or an aggregator
	if (signal_is_mediator(node))
	{
		return parse_mediator_signal(node);
	}
	else if (signal_is_aggregator(node))
	{
		return parse_aggregator_signal(node);
	}
	else if (signal_is_elementary(node))
	{
		return parse_elementary_signal(node);
	}
	std::cerr << "ERROR: Failed to parse node in behavior rulesets:" << std::endl;
	node.print(std::cerr);
	exit(-1);
}

bool signal_is_mediator(pugi::xml_node node)
{
	bool is_mediator = std::string(node.name()) == "behavior"; // by default, assume that the top-level signal is a mediator
	is_mediator |= (node.attribute("type")) && node.attribute("type").value() == "mediator";
	is_mediator |= xml_find_node(node, "decreasing_signals") != nullptr; // if the type was not declared but there are decreasing signals, then it is a mediator
	is_mediator |= xml_find_node(node, "increasing_signals") != nullptr; // if the type was not declared but there are increasing signals, then it is a mediator
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
	if (std::string(node.name()) == "signal")
	{
		return true; // if the node is a signal node, then it is an elementary signal by process of elimination
	}
	return false;
}

std::unique_ptr<AbstractSignal> parse_mediator_signal(pugi::xml_node mediator_node)
{
	double base_value = 1.0;
	if (std::string(mediator_node.name()) == "behavior")
	{
		std::string behavior = mediator_node.attribute("name").value();
		std::string cell_type = mediator_node.parent().attribute("name").value();
		base_value = get_single_base_behavior(find_cell_definition(cell_type), behavior);
	}
	else if (mediator_node.child("base_value"))
	{
		base_value = xml_get_double_value(mediator_node, "base_value");
	}
	pugi::xml_node decreasing_signals_node = mediator_node.child("decreasing_signals");
	pugi::xml_node increasing_signals_node = mediator_node.child("increasing_signals");

	double min_value = base_value; // default to base value unless the decreasing signals node has a max_response
	double max_value = base_value; // default to base value unless the increasing signals node has a max_response

	std::unique_ptr<AbstractSignal> pDecreasingSignal;
	std::unique_ptr<AbstractSignal> pIncreasingSignal;
	if (decreasing_signals_node)
	{
		pDecreasingSignal = parse_aggregator_signal(decreasing_signals_node);
		if (decreasing_signals_node.child("max_response"))
		{
			min_value = xml_get_double_value(decreasing_signals_node, "max_response");
		}
	}
	else
	{
		// later can check that at least one signal present, and if not, switch to an aggregator on the increasing signals
		pDecreasingSignal = std::unique_ptr<AggregatorSignal>(new AggregatorSignal());
	}
	if (increasing_signals_node)
	{
		pIncreasingSignal = parse_aggregator_signal(increasing_signals_node);
		if (increasing_signals_node.child("max_response"))
		{
			max_value = xml_get_double_value(increasing_signals_node, "max_response");
		}
	}
	else
	{
		// later can check that at least one signal present, and if not, switch to an aggregator on the decreasing signals
		pIncreasingSignal = std::unique_ptr<AggregatorSignal>(new AggregatorSignal());
	}

	std::unique_ptr<MediatorSignal> pMS = std::unique_ptr<MediatorSignal>(new MediatorSignal(std::move(pDecreasingSignal), std::move(pIncreasingSignal), min_value, base_value, max_value));

	pugi::xml_node mediator_fn_node = xml_find_node(mediator_node, "mediator");
	if (mediator_fn_node)
	{
		std::string mediator_fn = xml_get_my_string_value(mediator_fn_node);
		if (mediator_fn == "decreasing dominant")
		{
			// decreasing dominant mediator is the default
		}
		else if (mediator_fn == "increasing dominant")
		{
			pMS->use_increasing_dominant_mediator();
		}
		else if (mediator_fn == "neutral")
		{
			pMS->use_neutral_mediator();
		}
		else
		{
			std::cerr << "ERROR: Mediator not recognized: " << mediator_fn << std::endl;
			std::cerr << "Must be one of: decreasing dominant, increasing dominant, neutral" << std::endl;
			exit(-1);
		}
	}
	return pMS;
}

std::unique_ptr<AbstractSignal> parse_aggregator_signal(pugi::xml_node aggregator_node)
{
	pugi::xml_node signal_node = xml_find_node(aggregator_node, "signal");
	std::vector<std::unique_ptr<PhysiCell::AbstractSignal>> signals;
	while (signal_node)
	{
		signals.push_back(parse_abstract_signal(signal_node));
		signal_node = signal_node.next_sibling("signal");
	}
	return std::unique_ptr<AggregatorSignal>(new AggregatorSignal(std::move(signals)));
}

std::unique_ptr<AbstractSignal> parse_elementary_signal(pugi::xml_node elementary_node)
{
	std::string name = elementary_node.attribute("name").value();
	std::string type = elementary_node.attribute("type").value();
	std::unique_ptr<AbstractSignal> pAS;
	if (type == "Hill")
	{
		pAS = parse_hill_signal(name, elementary_node);
	}
	else if (type == "PartialHill")
	{
		pAS = parse_partial_hill_signal(name, elementary_node);
	}
	else if (type == "Linear")
	{
		pAS = parse_linear_signal(name, elementary_node);
	}
	else if (type == "Heaviside" || type == "Step" || type == "Switch")
	{
		pAS = parse_heaviside_signal(name, elementary_node);
	}
	else
	{
		std::cerr << "ERROR: Elementary signal type not recognized: " << type << std::endl;
		exit(-1);
	}
	pugi::xml_node reference_node = xml_find_node(elementary_node, "reference");
	if (reference_node)
	{
		// Use dynamic_cast to check if pAS is actually an ElementarySignal
		ElementarySignal *pES = dynamic_cast<ElementarySignal *>(pAS.get());
		if (pES)
		{
			parse_reference(reference_node, pES);
		}
		else
		{
			std::cerr << "ERROR: Signal is not an ElementarySignal." << std::endl;
			exit(-1);
		}
	}
	return pAS;
}

std::unique_ptr<AbstractSignal> parse_hill_signal(std::string name, pugi::xml_node hill_node)
{
	double half_max = xml_get_double_value(hill_node, "half_max");
	double hill_power = xml_get_double_value(hill_node, "hill_power");
	bool applies_to_dead_cells = xml_get_bool_value(hill_node, "applies_to_dead_cells");
	return std::unique_ptr<HillSignal>(new HillSignal(name, applies_to_dead_cells, half_max, hill_power));
}

std::unique_ptr<AbstractSignal> parse_partial_hill_signal(std::string name, pugi::xml_node partial_hill_node)
{
	double half_max = xml_get_double_value(partial_hill_node, "half_max");
	double hill_power = xml_get_double_value(partial_hill_node, "hill_power");
	bool applies_to_dead_cells = xml_get_bool_value(partial_hill_node, "applies_to_dead_cells");
	return std::unique_ptr<PartialHillSignal>(new PartialHillSignal(name, applies_to_dead_cells, half_max, hill_power));
}

std::unique_ptr<AbstractSignal> parse_linear_signal(std::string name, pugi::xml_node linear_node)
{
	double signal_min = xml_get_double_value(linear_node, "signal_min");
	double signal_max = xml_get_double_value(linear_node, "signal_max");
	bool applies_to_dead_cells = xml_get_bool_value(linear_node, "applies_to_dead_cells");
	return std::unique_ptr<LinearSignal>(new LinearSignal(name, applies_to_dead_cells, signal_min, signal_max));
}

std::unique_ptr<AbstractSignal> parse_heaviside_signal(std::string name, pugi::xml_node heaviside_node)
{
	double threshold = xml_get_double_value(heaviside_node, "threshold");
	bool applies_to_dead_cells = xml_get_bool_value(heaviside_node, "applies_to_dead_cells");
	return std::unique_ptr<HeavisideSignal>(new HeavisideSignal(name, applies_to_dead_cells, threshold));
}

void parse_reference(pugi::xml_node reference_node, ElementarySignal *pES)
{
	std::string type = xml_get_string_value(reference_node, "type");
	double reference_value = xml_get_double_value(reference_node, "value");
	std::unique_ptr<SignalReference> pSR;
	if (type == "increasing")
	{
		pSR = std::unique_ptr<SignalReference>(new IncreasingSignalReference(reference_value));
	}
	else if (type == "decreasing")
	{
		pSR = std::unique_ptr<SignalReference>(new DecreasingSignalReference(reference_value));
	}
	else
	{
		std::cerr << "ERROR: Reference type not recognized: " << type << std::endl;
		exit(-1);
	}
    pES->add_reference(std::move(pSR));
	return;
}

void parse_behavior_rules_from_pugixml( void )
{
	pugi::xml_node node = physicell_config_root.child( "cell_rules" ); 
	if( !node )
	{ 
		std::cout << "Warning: Could not find <cell_rules> section of XML config file." << std::endl 
				 <<  "       Cannot parse cell rules, so disabling." << std::endl; 

		PhysiCell_settings.rules_enabled = false; 
		return; 
	}

	// find the set of rulesets 
	node = node.child( "rulesets" ); 
	if( !node )
	{ 
		std::cout << "Warning: Could not find <rulesets> in the <cell_rules> section of XML config file." << std::endl 
				 <<  "       Cannot parse cell rules, so disabling." << std::endl; 

		PhysiCell_settings.rules_enabled = false; 
		return; 
	}
	// find the first ruleset 
	node = node.child( "ruleset");
	if( !node )
	{ 
		std::cout << "Warning: Could not find any <ruleset> in the <rulesets> section of XML config file." << std::endl 
				 <<  "       Cannot parse cell rules, so disabling." << std::endl; 

		PhysiCell_settings.rules_enabled = false; 
		return; 
	}

	while( node )
	{
		if (node.attribute("enabled").as_bool() == true)
		{
			std::string folder = xml_get_string_value(node, "folder");
			std::string filename = xml_get_string_value(node, "filename");
			std::string input_filename = folder + "/" + filename;

			std::string format = node.attribute("format").as_string();
			std::string protocol = node.attribute("protocol").as_string();
			double version = node.attribute("version").as_double();

			parse_behavior_rules_from_file(input_filename, format, protocol, version);
		}
		else
		{
			std::cout << "\tRuleset disabled ... " << std::endl;
		}
		node = node.next_sibling("ruleset");
	}
	return;
}

void parse_behavior_rules_from_file(std::string path_to_file, std::string format, std::string protocol, double version) // see PhysiCell_rules.h for default values of format, protocol, and version
{
	std::cout << "\tProcessing ruleset in " << path_to_file << " ... " << std::endl;

	// get the file  extension of path_to_file
	if (format == "")
	{
		format = path_to_file.substr(path_to_file.find_last_of(".") + 1);
	}
	if (protocol == "")
	{
		protocol = "CBHG"; // default protocol (at least for CSVs)
	}
	if (version == -1.0)
	{
		version = 3.0; // default version (at least for CSVs)
	}

	if (format == "CSV" || format == "csv")
	{
		std::cerr << "Error: CSV format not supported for behavior rules. Quitting!" << std::endl;
		throw std::runtime_error("CSV format not supported for behavior rules.");
	}
	else if (format == "XML" || format == "xml")
	{
		std::cout << "\tFormat: XML" << std::endl;
		parse_xml_behavior_rules(path_to_file);
		PhysiCell_settings.rules_enabled = true;
	}
	else
	{
		std::cerr << "\tError: Unknown format (" << format << ") for ruleset " << path_to_file << ". Quitting!" << std::endl;
		exit(-1);
	}
	PhysiCell_settings.rules_enabled = true;
	return; 
}

};