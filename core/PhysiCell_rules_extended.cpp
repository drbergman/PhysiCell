#include "PhysiCell_rules_extended.h"
#include <algorithm> // For std::min_element, std::max_element, std::nth_element

namespace PhysiCell{

#ifndef __PhysiCell_rules_extended_cpp__
#define __PhysiCell_rules_extended_cpp__
#endif

double first_aggregator(std::vector<double> signals_in)
{
	return signals_in[0];
}

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

double product_aggregator(std::vector<double> signals_in)
{
	double signal_product = 1;
	for (auto &signal : signals_in)
	{
		if (signal == 0)
		{
			return 0;
		}
		signal_product *= signal;
	}
	return signal_product;
}

double mean_aggregator(std::vector<double> signals_in)
{
	return sum_aggregator(signals_in) / signals_in.size();
}

double max_aggregator(std::vector<double> signals_in)
{
    return *std::max_element(signals_in.begin(), signals_in.end());
}

double min_aggregator(std::vector<double> signals_in)
{
	return *std::min_element(signals_in.begin(), signals_in.end());
}

double median_aggregator(std::vector<double> signals_in)
{
    size_t n = signals_in.size();
    size_t mid = n / 2;

    std::nth_element(signals_in.begin(), signals_in.begin() + mid, signals_in.end());

    double signal_median;
    if (n % 2 == 0)
    {
        double mid1 = signals_in[mid];
        std::nth_element(signals_in.begin(), signals_in.begin() + mid - 1, signals_in.end());
        double mid2 = signals_in[mid - 1];
        signal_median = (mid1 + mid2) / 2.0;
    }
    else
    {
        signal_median = signals_in[mid];
    }
    return signal_median;
}

double geometric_mean_aggregator(std::vector<double> signals_in)
{
	double product = product_aggregator(signals_in);
	return pow(product, 1.0 / signals_in.size());
}

void AggregatorSignal::set_aggregator(std::string aggregator_name)
{
	if (aggregator_name == "sum")
	{
		aggregator = sum_aggregator;
	}
	else if (aggregator_name == "product")
	{
		aggregator = product_aggregator;
	}
	else if (aggregator_name == "multivariate hill" || aggregator_name == "multivariate_hill")
	{
		aggregator = multivariate_hill_aggregator;
	}
	else if (aggregator_name == "mean")
	{
		aggregator = mean_aggregator;
	}
	else if (aggregator_name == "max")
	{
		aggregator = max_aggregator;
	}
	else if (aggregator_name == "min")
	{
		aggregator = min_aggregator;
	}
	else if (aggregator_name == "median")
	{
		aggregator = median_aggregator;
	}
	else if (aggregator_name == "geometric mean" || aggregator_name == "geometric_mean")
	{
		aggregator = geometric_mean_aggregator;
	}
	else if (aggregator_name == "first")
	{
		aggregator = first_aggregator;
	}
	else if (aggregator_name == "custom")
	{
		aggregator = [](std::vector<double> signals_in)
		{
			std::cerr << "ERROR: Custom aggregator not set! Make sure to set one in custom.cpp if using a custom aggregator." << std::endl;
			exit(-1);
			return -1;
		};
	}
	else
	{
		std::cerr << "ERROR: Aggregator not recognized: " << aggregator_name << std::endl;
		std::cerr << "Must be one of: 'sum', 'product', 'multivariate hill', 'mean', 'max', 'min', 'median', 'geometric mean', 'custom'" << std::endl;
		exit(-1);
	}
	type = aggregator_name;
	return;
}

void AggregatorSignal::display(std::ostream &os, RuleLine rule_line, int indent, std::string additional_info)
{
	if (!has_signals())
	{
		return;
	}

	os << "// " << std::string(indent*2, ' ') << "└─" << additional_info << " using " << type << " aggregator" << std::endl;
	for (auto &signal : signals)
	{
		signal->display(os, rule_line, indent + 1);
	}
	return;
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

void MediatorSignal::set_mediator(std::string mediator_name)
{
	if (mediator_name == "decreasing dominant" || mediator_name == "decreasing_dominant")
	{
		aggregator = [this](std::vector<double> signals_in)
		{
			return this->decreasing_dominant_mediator(signals_in);
		};
	}
	else if (mediator_name == "increasing dominant" || mediator_name == "increasing_dominant")
	{
		aggregator = [this](std::vector<double> signals_in)
		{
			return this->increasing_dominant_mediator(signals_in);
		};
	}
	else if (mediator_name == "neutral" || mediator_name == "neutral_mediator" || mediator_name == "neutral mediator")
	{
		aggregator = [this](std::vector<double> signals_in)
		{
			return this->neutral_mediator(signals_in);
		};
	}
	else if (mediator_name == "custom")
	{
		aggregator = [](std::vector<double> signals_in)
		{
			std::cerr << "ERROR: Custom mediator not set! Make sure to set one in custom.cpp if using a custom mediator." << std::endl;
			exit(-1);
			return -1;
		};
	}
	else
	{
		std::cerr << "ERROR: Mediator not recognized: " << mediator_name << std::endl
				  << "Must be one of: 'decreasing dominant', 'increasing dominant', 'neutral', 'custom'" << std::endl;
		exit(-1);
	}
	type = mediator_name;
	return;
}

void MediatorSignal::display(std::ostream &os, RuleLine rule_line, int indent, std::string additional_info)
{
	os << "// " << std::string(indent*2, ' ') << "└─mediating between " << min_value << " ≤ " << base_value << " ≤ " << max_value << " using a " << type << " mediator" << std::endl;
	rule_line.response = "decreases";
	rule_line.max_response = min_value;
	decreasing_signal->display(os, rule_line, indent + 1, "decreasing");
	rule_line.response = "increases";
	rule_line.max_response = max_value;
	increasing_signal->display(os, rule_line, indent + 1, "increasing");
	return;
}

double RelativeSignal::evaluate(Cell *pCell)
{
	if (!applies_to_dead && pCell->phenotype.death.dead)
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

std::string RelativeSignal::construct_relative_signal_string( void )
{
	std::string signal;
	if (has_reference())
	{
		signal = "(" + signal_reference->get_type() + ") " + signal_name + " (from " + std::to_string(signal_reference->get_reference_value()) + ")";
	}
	else
	{
		signal = signal_name;
	}
	return signal;
}

void PartialHillSignal::display(std::ostream &os, RuleLine rule_line, int indent, std::string additional_info)
{
	os << rule_line.cell_type << "," << construct_relative_signal_string() << "," << rule_line.response << "," << rule_line.behavior << "," << rule_line.max_response << "," << half_max << "," << hill_power << "," << applies_to_dead << std::endl;
	return;
}

void HillSignal::display(std::ostream &os, RuleLine rule_line, int indent, std::string additional_info)
{
	os << rule_line.cell_type << "," << construct_relative_signal_string() << "," << rule_line.response << " (hill)," << rule_line.behavior << "," << rule_line.max_response << "," << half_max << "," << hill_power << "," << applies_to_dead << std::endl;
	return;
}

void IdentitySignal::display(std::ostream &os, RuleLine rule_line, int indent, std::string additional_info)
{
	// Note: the two commas after max_response are intentional to make sure this has 8 columns
	os << rule_line.cell_type << "," << construct_relative_signal_string() << "," << rule_line.response << " (identity)," << rule_line.behavior << "," << rule_line.max_response << ",,," << applies_to_dead << std::endl;
	return;
}

double AbsoluteSignal::evaluate(Cell *pCell)
{
	if (!applies_to_dead && pCell->phenotype.death.dead)
	{
		return 0;
	}
	return transformer(get_single_signal(pCell, signal_name));
}

std::string AbsoluteSignal::construct_absolute_signal_string( void )
{
	return "(" + type + ") " + signal_name;
}

void LinearSignal::display(std::ostream &os, RuleLine rule_line, int indent, std::string additional_info)
{
	os << rule_line.cell_type << "," << construct_absolute_signal_string() << "," << rule_line.response << " (linear)," << rule_line.behavior << "," << rule_line.max_response << "," << signal_min << "," << signal_max << "," << applies_to_dead << std::endl;
	return;
}

void HeavisideSignal::display(std::ostream &os, RuleLine rule_line, int indent, std::string additional_info)
{
	// Note: the two commas after threshold are intentional to make sure this has 8 columns
	os << rule_line.cell_type << "," << construct_absolute_signal_string() << "," << rule_line.response << " (heaviside)," << rule_line.behavior << "," << rule_line.max_response << "," << threshold << ",," << applies_to_dead << std::endl;
	return;
}

void BehaviorSetter::apply(Cell *pCell)
{
	double param = signal->evaluate(pCell);
	set_single_behavior(pCell, behavior, param);
}

double euler_direct_solve(double current, double rate, double target)
{
	return current + phenotype_dt * rate * (target - current);
}

double exponential_solve(double current, double rate, double target)
{
	return target + (current - target) * exp(-rate * phenotype_dt);
}

void BehaviorSetter::display(std::ostream &os, RuleLine rule_line)
{
	os << "//   └─set " << behavior << std::endl;
	rule_line.behavior = behavior;
	signal->display(os, rule_line, 2);
	return;
}

void BehaviorAccumulator::display(std::ostream &os, RuleLine rule_line)
{
	os << "//   └─accumulate " << behavior << " from " << behavior_base << " to " << behavior_saturation << std::endl;
	rule_line.behavior = behavior;
	signal->display(os, rule_line, 2);
	return;
}

void BehaviorAttenuator::display(std::ostream &os, RuleLine rule_line)
{
	os << "//   └─attenuate " << behavior << " from " << behavior_saturation << " to " << behavior_base << std::endl;
	rule_line.behavior = behavior;
	signal->display(os, rule_line, 2);
	return;
}

void BehaviorAccumulator::apply(Cell *pCell)
{
	double rate = signal->evaluate(pCell);
	double current_value = get_single_behavior(pCell, behavior);
	if (rate < 0)
	{
		// current_value = std::max(euler_direct_solve(current_value, -rate, behavior_base), behavior_base); // decay towards base value (note rate < 0); do not let it go below base value
		current_value = exponential_solve(current_value, -rate, behavior_base);
	}
	else if (rate > 0)
	{
		// current_value = std::min(euler_direct_solve(current_value, rate, behavior_saturation), behavior_saturation); // grow towards saturation value (note rate > 0); do not let it go above saturation value
		current_value = exponential_solve(current_value, rate, behavior_saturation);
	}
	set_single_behavior(pCell, behavior, current_value);
}

void BehaviorAttenuator::apply(Cell *pCell)
{
	double rate = signal->evaluate(pCell);
	double current_value = get_single_behavior(pCell, behavior);
	if (rate < 0)
	{
		// current_value = std::min(euler_direct_solve(current_value, -rate, behavior_base), behavior_base); // decay towards base value (note rate < 0); do not let it go below saturation value
		current_value = exponential_solve(current_value, -rate, behavior_base);
	}
	else if (rate > 0)
	{
		// current_value = std::max(euler_direct_solve(current_value, rate, behavior_saturation), behavior_saturation); // grow towards saturation value (note rate > 0); do not let it go above base value
		current_value = exponential_solve(current_value, rate, behavior_saturation);
	}
	set_single_behavior(pCell, behavior, current_value);
}

void BehaviorRuleset::apply(Cell *pCell)
{
	for (auto &rule : rules)
	{
		rule->apply(pCell);
	}
	return;
}

void BehaviorRuleset::display(std::ostream &os, RuleLine rule_line)
{
	for (auto &rule : rules)
	{
		rule->display(os, rule_line);
	}
	return;
}

std::unordered_map<Cell_Definition *, std::unique_ptr<BehaviorRuleset>> behavior_rulesets;

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
	for (auto &pCD : cell_definitions_by_index)
	{
		add_behavior_ruleset(pCD); 
	}
	return; 
}

void setup_behavior_rules( void )
{
	intialize_behavior_rulesets();

	// load rules
	if (argument_parser.path_to_rules_file != "")
	{parse_behavior_rules_from_file(argument_parser.path_to_rules_file);} // see function signature in PhysiCell_rules_extended.h for default values for rest of arguments
	else
	{parse_behavior_rules_from_pugixml();}

	display_behavior_rulesets( std::cout );

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
	std::string rules_file = PhysiCell_settings.folder + "/cell_rules_parsed.csv"; 
	export_behavior_rules( rules_file ); 
	return; 
}

BehaviorRuleset* find_behavior_ruleset( Cell_Definition* pCD )
{ return behavior_rulesets[pCD].get(); }

void apply_behavior_ruleset( Cell* pCell )
{
	Cell_Definition* pCD = find_cell_definition( pCell->type_name ); 
	behavior_rulesets[pCD]->apply( pCell );
	return; 
}

void parse_xml_behavior_rules(const std::string filename)
{
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

	pugi::xml_node ruleset_node;
	ruleset_node = xml_find_node(physicell_rules_root, "behavior_ruleset");
	std::vector<std::string> cell_definitions_ruled;

	while (ruleset_node)
	{
		std::string cell_type = ruleset_node.attribute("name").value();
		// check if cell_type has already had a ruleset defined for it
		if (std::find(cell_definitions_ruled.begin(), cell_definitions_ruled.end(), cell_type) != cell_definitions_ruled.end())
		{
			std::cerr << "XML Rules ERROR: Two rulesets for " << cell_type << " found." << std::endl
					  << "\tCombine them into a single ruleset please :)" << std::endl
					  << "\tSupport for rules across multiple files for the same cell type not yet supported." << std::endl;
			exit(-1);
		}
		cell_definitions_ruled.push_back(cell_type);

		std::vector<std::string> behaviors_ruled;
		Cell_Definition *pCD = find_cell_definition(cell_type);

		pugi::xml_node behavior_node = ruleset_node.child("behavior");
		while (behavior_node)
		{
			std::string behavior = behavior_node.attribute("name").value();
			if (std::find(behaviors_ruled.begin(), behaviors_ruled.end(), behavior) != behaviors_ruled.end())
			{
				std::cerr << "XML Rules ERROR: The behavior " << behavior << " is being set again for " << cell_type << "." << std::endl
						  << "\tSelect only one to keep. :)" << std::endl;
				exit(-1);
			}
			behaviors_ruled.push_back(behavior);

			std::unique_ptr<BehaviorRule> pBR = parse_behavior(cell_type, behavior, behavior_node);
			behavior_rulesets[pCD]->add_behavior_rule(std::move(pBR));

			behavior_node = behavior_node.next_sibling("behavior");
		}
		ruleset_node = ruleset_node.next_sibling("behavior_ruleset");
	}
	return;
}

std::unique_ptr<BehaviorRule> parse_behavior(std::string cell_type, std::string behavior, pugi::xml_node node)
{
	std::string type = "setter"; // default to setter
	if (xml_find_node(node, "type"))
	{
		type = xml_get_string_value(node, "type");
	}
	std::unique_ptr<AbstractSignal> pAS = parse_abstract_signal(node);
	if (type=="setter")
	{
		return std::unique_ptr<BehaviorRule>(new BehaviorSetter(behavior, std::move(pAS)));
	}
	double behavior_base;
	if (xml_find_node(node, "behavior_base"))
	{
		behavior_base = xml_get_double_value(node, "behavior_base");
	}
	else
	{
		behavior_base = get_single_base_behavior(find_cell_definition(cell_type), behavior); // base value for the behavior (not the rate) that negative rates will return the cell to
	}
	double behavior_saturation = xml_get_double_value(node, "behavior_saturation"); // saturation value for the behavior (not the rate) that positive rates will move the cell to
	if (type=="accumulator")
	{
		return std::unique_ptr<BehaviorRule>(new BehaviorAccumulator(behavior, std::move(pAS), behavior_base, behavior_saturation));
	}
	else if (type=="attenuator")
	{
		return std::unique_ptr<BehaviorRule>(new BehaviorAttenuator(behavior, std::move(pAS), behavior_base, behavior_saturation));
	}
	std::cerr << "ERROR: Behavior type not recognized: " << type << std::endl;
	std::cerr << "Must be one of: setter (or omitted), accumulator, attenuator" << std::endl;
	exit(-1);
}

std::unique_ptr<AbstractSignal> parse_abstract_signal(pugi::xml_node node)
{
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
	is_mediator |= (node.attribute("type")) && std::strcmp(node.attribute("type").value(), "mediator")==0;
	is_mediator |= xml_find_node(node, "decreasing_signals") != nullptr; // if the type was not declared but there are decreasing signals, then it is a mediator
	is_mediator |= xml_find_node(node, "increasing_signals") != nullptr; // if the type was not declared but there are increasing signals, then it is a mediator
	return is_mediator;
}

bool signal_is_aggregator(pugi::xml_node node)
{
	bool is_aggregator = (node.attribute("type")) && std::strcmp(node.attribute("type").value(), "aggregator")==0;
	is_aggregator |= std::string(node.name()) == "decreasing_signals" || std::string(node.name()) == "increasing_signals";
	return is_aggregator;
}

bool signal_is_elementary(pugi::xml_node node)
{
	return std::string(node.name()) == "signal" && node.attribute("name"); // make sure it fits the template for an elementary signal. the type attribute can be omitted, defaults to PartialHill
}

std::unique_ptr<AbstractSignal> parse_mediator_signal(pugi::xml_node mediator_node)
{
	double base_value = 1.0;
	if (mediator_node.child("base_value"))
	{
		base_value = xml_get_double_value(mediator_node, "base_value");
	}
	else if (std::string(mediator_node.name()) == "behavior")
	{
		std::string behavior = mediator_node.attribute("name").value();
		std::string cell_type = mediator_node.parent().attribute("name").value();
		base_value = get_single_base_behavior(find_cell_definition(cell_type), behavior);
	}
	else
	{
		std::cerr << "ERROR: Mediator signal must have a base value." << std::endl;
		exit(-1);
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
		if (decreasing_signals_node.child("max_response")) // note that even if the min value = base value, a decreasing signal can counteract an increasing signal
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
		if (increasing_signals_node.child("max_response")) // note that even if the max value = base value, an increasing signal can counteract a decreasing signal
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
		pMS->set_mediator(mediator_fn);
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
	std::unique_ptr<AggregatorSignal> pAS = std::unique_ptr<AggregatorSignal>(new AggregatorSignal(std::move(signals)));

	pugi::xml_node aggregator_fn_node = xml_find_node(aggregator_node, "aggregator");
	if (aggregator_fn_node)
	{
		std::string aggregator_fn = xml_get_my_string_value(aggregator_fn_node);
		pAS->set_aggregator(aggregator_fn);
	}
	return pAS;
}

std::unique_ptr<AbstractSignal> parse_elementary_signal(pugi::xml_node elementary_node)
{
	std::string name = elementary_node.attribute("name").value();
	std::string type = "partial_hill"; // default to partial hill
	if (elementary_node.attribute("type"))
	{ type = elementary_node.attribute("type").value(); }
	std::unique_ptr<AbstractSignal> pAS;
	if (type == "Hill" || type == "hill")
	{
		pAS = parse_hill_signal(name, elementary_node);
	}
	else if (type == "PartialHill" || type == "partialhill" || type == "partial_hill" || type == "partial hill")
	{
		pAS = parse_partial_hill_signal(name, elementary_node);
	}
	else if (type == "Identity" || type == "identity")
	{
		pAS = parse_identity_signal(name, elementary_node);
	}
	else if (type == "Linear" || type == "linear")
	{
		pAS = parse_linear_signal(name, elementary_node);
	}
	else if (type == "Heaviside" || type == "Step" || type == "Switch" || type == "heaviside" || type == "step" || type == "switch")
	{
		pAS = parse_heaviside_signal(name, elementary_node);
	}
	else
	{
		std::cerr << "ERROR: Elementary signal type not recognized: " << type << std::endl
		          << "Must be one of: 'Hill', 'PartialHill', 'Linear', 'Heaviside'."
				  << "Note: each of these options has alternative capitalizations and/or synonyms.";
		exit(-1);
	}
	pugi::xml_node reference_node = xml_find_node(elementary_node, "reference");
	if (reference_node)
	{
		// Use dynamic_cast to check if pAS is actually a RelativeSignal
		RelativeSignal *pES = dynamic_cast<RelativeSignal *>(pAS.get());
		if (pES)
		{
			parse_reference(reference_node, pES);
		}
		else
		{
			std::cerr << "ERROR: Reference defined for a Signal that is not a relative signal." << std::endl
					  << "Relative signals include: 'PartialHill' and 'Hill'.";
			exit(-1);
		}
	}
	return pAS;
}

std::unique_ptr<AbstractSignal> parse_hill_signal(std::string name, pugi::xml_node hill_node)
{
	double half_max = xml_get_double_value(hill_node, "half_max");
	double hill_power = xml_get_double_value(hill_node, "hill_power");
	bool applies_to_dead = xml_get_bool_value(hill_node, "applies_to_dead");
	return std::unique_ptr<HillSignal>(new HillSignal(name, applies_to_dead, half_max, hill_power));
}

std::unique_ptr<AbstractSignal> parse_partial_hill_signal(std::string name, pugi::xml_node partial_hill_node)
{
	double half_max = xml_get_double_value(partial_hill_node, "half_max");
	double hill_power = xml_get_double_value(partial_hill_node, "hill_power");
	bool applies_to_dead = xml_get_bool_value(partial_hill_node, "applies_to_dead");
	return std::unique_ptr<PartialHillSignal>(new PartialHillSignal(name, applies_to_dead, half_max, hill_power));
}

std::unique_ptr<AbstractSignal> parse_identity_signal(std::string name, pugi::xml_node identity_node)
{
	bool applies_to_dead = xml_get_bool_value(identity_node, "applies_to_dead");
	return std::unique_ptr<IdentitySignal>(new IdentitySignal(name, applies_to_dead));
}

std::unique_ptr<AbstractSignal> parse_linear_signal(std::string name, pugi::xml_node linear_node)
{
	bool is_decreasing = xml_find_node(linear_node, "type") && xml_get_string_value(linear_node, "type") == "decreasing"; // default is increasing, so only switch if this is given and is set to "decreasing"
	double signal_min = xml_get_double_value(linear_node, "signal_min");
	double signal_max = xml_get_double_value(linear_node, "signal_max");
	bool applies_to_dead = xml_get_bool_value(linear_node, "applies_to_dead");
	std::unique_ptr<LinearSignal> pLS;
	if (is_decreasing)
	{
		pLS = std::unique_ptr<LinearSignal>(new DecreasingLinearSignal(name, applies_to_dead, signal_min, signal_max));
	}
	else
	{
		pLS = std::unique_ptr<LinearSignal>(new IncreasingLinearSignal(name, applies_to_dead, signal_min, signal_max));
	}
	return pLS;
}

std::unique_ptr<AbstractSignal> parse_heaviside_signal(std::string name, pugi::xml_node heaviside_node)
{
	double threshold = xml_get_double_value(heaviside_node, "threshold");
	bool applies_to_dead = xml_get_bool_value(heaviside_node, "applies_to_dead");
	bool is_decreasing = xml_find_node(heaviside_node, "type") && xml_get_string_value(heaviside_node, "type") == "decreasing"; // default is increasing, so only switch if this is given and is set to "decreasing"
	std::unique_ptr<HeavisideSignal> pHS;
	if (is_decreasing)
	{
		pHS = std::unique_ptr<HeavisideSignal>(new DecreasingHeavisideSignal(name, applies_to_dead, threshold));
	}
	else
	{
		pHS = std::unique_ptr<HeavisideSignal>(new IncreasingHeavisideSignal(name, applies_to_dead, threshold));
	}
	return pHS;
}

void parse_reference(pugi::xml_node reference_node, RelativeSignal *pRelSig)
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
    pRelSig->add_reference(std::move(pSR));
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
			std::cout << "\tA ruleset is disabled ... " << std::endl;
		}
		node = node.next_sibling("ruleset");
	}
	return;
}

void parse_behavior_rules_from_file(std::string path_to_file, std::string format, std::string protocol, double version) // see PhysiCell_rules_extended.h for default values of format, protocol, and version
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
		std::cout << "\tFormat: CSV" << std::endl;
		parse_csv_behavior_rules(path_to_file, protocol, version);
	}
	else if (format == "XML" || format == "xml")
	{
		std::cout << "\tFormat: XML" << std::endl;
		parse_xml_behavior_rules(path_to_file);
	}
	else
	{
		std::cerr << "\tError: Unknown format (" << format << ") for ruleset " << path_to_file << ". Quitting!" << std::endl;
		exit(-1);
	}
	PhysiCell_settings.rules_enabled = true;
	copy_file_to_output(path_to_file);
	return; 
}

// loading rules from CSV files
void parse_csv_behavior_rules(std::string path_to_file, std::string protocol, double version)
{
	std::fstream fs( path_to_file, std::ios::in );
	if( !fs )
	{
		std::cout << "Warning: Rules file " << path_to_file << " failed to open." << std::endl; 
		return; 
	}

	std::cout << "Processing rules in file " << path_to_file << " ... " << std::endl; 

	while( fs.eof() == false )
	{
		std::string line; 	
		std::getline( fs , line, '\n'); 
		if( line.size() > 0 )
		{ parse_csv_behavior_rule(line); }
	}

	fs.close(); 

	std::cout << "Done!" << std::endl << std::endl; 

	return; 
}

void parse_csv_behavior_rule(std::string line)
{
	std::vector<std::string> tokenized_string; 
	split_csv( line , tokenized_string , ','); 

	// Make sure it was truly comma-separated. 
	// If not, try tab.
	if( tokenized_string.size() != 8 )
	{ split_csv( line , tokenized_string , '\t');  }

	// check for comment 
	if(tokenized_string[0][0] == '/' && tokenized_string[0][1] == '/' )
	{ std::cout << "Skipping commented rule (" << line << ")" << std::endl; return; }	

	if (is_csv_rule_misformed(tokenized_string))
	{
		std::cout << "Warning: Misformed rule (likely from an empty rules file). Skipping." << std::endl; 

		for( int n=0 ; n < tokenized_string.size(); n++ )
		{
			std::cout << n << " : " << tokenized_string[n] << std::endl; 
		}

		return; 
	}
	std::string temp = csv_strings_to_English_v3(tokenized_string, false); // need a v1 version of this
	// string portions of the rule
	std::string cell_type = tokenized_string[0];
	std::string signal = tokenized_string[1];
	std::string response = tokenized_string[2];
	std::string behavior = tokenized_string[3];

	double max_response = std::atof( tokenized_string[4].c_str() ); 
	double half_max  = std::atof( tokenized_string[5].c_str() );
	double hill_power = std::atof( tokenized_string[6].c_str() );
	bool use_for_dead = (bool) std::atof( tokenized_string[7].c_str() );

	std::cout << "Adding rule for " << cell_type << " cells:\n\t";
	std::cout << temp << std::endl;

	Cell_Definition* pCD = find_cell_definition(cell_type);
	if (pCD == NULL)
	{
		std::cout << "Error: Cell type " << cell_type << " not found." << std::endl;
		exit(-1);
	}
	double ref_base_value = get_single_base_behavior(pCD, behavior);

	BehaviorRule *pBR = behavior_rulesets[pCD]->find_behavior(behavior);

	if (pBR == nullptr)
	{
		std::unique_ptr<BehaviorSetter> newRule(new BehaviorSetter(behavior));
		pBR = newRule.get();
		behavior_rulesets[pCD]->add_behavior_rule(std::move(newRule));
	}
	
	MediatorSignal* pMS = dynamic_cast<MediatorSignal*>(pBR->signal.get());
	pMS->set_base_value(ref_base_value);
	AggregatorSignal* pAS = nullptr;
	if (response == "decreases")
	{
		pAS = dynamic_cast<AggregatorSignal*>(pMS->get_decreasing_signal());
		pMS->set_min_value(max_response); // set the min value for the decreasing signal
	}
	else
	{
		pAS = dynamic_cast<AggregatorSignal*>(pMS->get_increasing_signal());
		pMS->set_max_value(max_response); // set the min value for the decreasing signal
	}
	std::unique_ptr<PartialHillSignal> pPHS(new PartialHillSignal(signal, use_for_dead, half_max, hill_power));
	pAS->add_signal(std::move(pPHS));
	return;
}

bool is_csv_rule_misformed(std::vector<std::string> input)
{
	if (input.size() != 8)
	{
		return true;
	}
	// if any empty strings, skip
	for (int n = 0; n < input.size(); n++)
	{
		if (input[n].empty() == true)
		{
			return true;
		}
	}
	return false;
}

void split_csv( std::string input , std::vector<std::string>& output , char delim )
{
	output.resize(0); 

	std::istringstream is(input);
	std::string part;
	while( getline(is, part, delim ) )
	{ output.push_back(part); }

	return; 
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

BehaviorRule *BehaviorRuleset::find_behavior(std::string behavior)
{
	for (auto &rule : rules)
	{
		if (rule->behavior == behavior)
		{
			return rule.get();
		}
	}
	return nullptr;
}

void set_custom_mediator(const std::string &cell_definition_name, const std::string &behavior_name, double (*mediator_function)(std::vector<double>))
{
	MediatorSignal *pMS = get_top_level_mediator(cell_definition_name, behavior_name);
	pMS->aggregator = mediator_function;
	return;
}

void set_custom_mediator(const std::string &cell_definition_name, const std::string &behavior_name, double (*mediator_function)(MediatorSignal*, std::vector<double>))
{
	MediatorSignal *pMS = get_top_level_mediator(cell_definition_name, behavior_name);
	pMS->aggregator = [pMS, mediator_function](std::vector<double> signals_in)
	{
		return mediator_function(pMS, signals_in);
	};
	return;
}

void set_custom_aggregator(const std::string &cell_definition_name, const std::string &behavior_name, const std::string &response, double (*aggregator_function)(std::vector<double>))
{
	MediatorSignal *pMS = get_top_level_mediator(cell_definition_name, behavior_name);
	AggregatorSignal *pAS;
	if (response == "decreases")
	{
		pAS = dynamic_cast<AggregatorSignal *>(pMS->get_decreasing_signal());
	}
	else if (response == "increases")
	{
		pAS = dynamic_cast<AggregatorSignal *>(pMS->get_increasing_signal());
	}
	else
	{
		std::cerr << "Error: Response type not recognized: " << response << std::endl;
		exit(-1);
	}
	pAS->aggregator = aggregator_function;
	return;
}

MediatorSignal* get_top_level_mediator(const std::string &cell_definition_name, const std::string &behavior_name)
{
	Cell_Definition *pCD = find_cell_definition(cell_definition_name);
	if (pCD == NULL)
	{
		std::cerr << "Error: Cell type " << cell_definition_name << " not found." << std::endl;
		exit(-1);
	}
	BehaviorRuleset *pBR = find_behavior_ruleset(pCD);
	if (pBR == NULL)
	{
		std::cerr << "Error: Behavior ruleset for cell type " << cell_definition_name << " not found." << std::endl;
		exit(-1);
	}
	BehaviorRule* pBRule = pBR->find_behavior(behavior_name);
	if (pBRule == nullptr)
	{
		std::cerr << "Error: Behavior " << behavior_name << " not found in ruleset for cell type " << cell_definition_name << "." << std::endl;
		exit(-1);
	}
	auto pMS = dynamic_cast<MediatorSignal *>(pBRule->signal.get());
	if (pMS == nullptr)
	{
		std::cerr << "Error: Signal is not of type MediatorSignal for cell definition '"
				  << cell_definition_name << "' and behavior '" << behavior_name << "'." << std::endl;
		exit(-1);
	}
	return pMS;
}

void display_behavior_rulesets(std::ostream &os)
{

	RuleLine rule_line;

	for (int n = 0; n < cell_definitions_by_index.size(); n++)
	{
		Cell_Definition *pCD = cell_definitions_by_index[n];
		BehaviorRuleset *pBR = find_behavior_ruleset(pCD);

		std::string cell_type = pCD->name;

		os << "// " << cell_type << std::endl;

		rule_line.cell_type = cell_type;

		pBR->display(os, rule_line);

		os << std::endl;
	}

	return;
}

/*
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

void save_annotated_English_rules( void )
{
	std::string filename = PhysiCell_settings.folder + "/rules.txt";
	std::ofstream of( filename , std::ios::out );
	stream_annotated_English_rules( of ); 
	of.close(); 
}

void save_annotated_English_rules_HTML( void )
{
	std::string filename = PhysiCell_settings.folder + "/rules.html";
	std::ofstream of( filename , std::ios::out );
	stream_annotated_English_rules_HTML( of ); 
	of.close(); 
}
*/

void export_behavior_rules( std::string filename )
{
	std::fstream fs( filename, std::ios::out );
	if( !fs )
	{
		std::cerr << "ERROR: Rules export file " << filename << " failed to open." << std::endl; 
		exit(-1); 
	}

	std::cout << "Exporting rules to file " << filename << std::endl; 
	display_behavior_rulesets(fs);
	
	fs.close(); 

	return; 
}
};