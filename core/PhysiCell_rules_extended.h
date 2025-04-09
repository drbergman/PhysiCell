#include <vector>
#include <string>
#include <functional>
#include <memory>


#ifndef __PhysiCell_rules_extended__
#define __PhysiCell_rules_extended__

#include "./PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 

namespace PhysiCell{

// aggregator functions
double sum_aggregator(std::vector<double> signals_in);
double multivariate_hill_aggregator(std::vector<double> partial_hill_signals);
double product_aggregator(std::vector<double> signals_in);
double mean_aggregator(std::vector<double> signals_in);
double max_aggretator(std::vector<double> signals_in);
double min_aggregator(std::vector<double> signals_in);
double median_aggregator(std::vector<double> signals_in);
double geometric_mean_aggregator(std::vector<double> signals_in);

class AbstractSignal
{
public:
    virtual double evaluate(Cell *pCell) = 0;

    virtual ~AbstractSignal() {}
};

class AbstractAggregatorSignal : public AbstractSignal
{
private:
    virtual std::vector<double> get_signals(Cell *pCell) = 0;
    
public:
    std::function<double(std::vector<double>)> aggregator;

    double evaluate(Cell *pCell)
    {
        return aggregator(get_signals(pCell));
    }

    virtual ~AbstractAggregatorSignal() {}
};

class AggregatorSignal : public AbstractAggregatorSignal
{
private:
    std::vector<std::unique_ptr<PhysiCell::AbstractSignal>> signals;

public:
    std::vector<double> get_signals(Cell *pCell)
    {
        std::vector<double> signal_values(signals.size());
        for (size_t i = 0; i < signals.size(); ++i)
        {
            signal_values[i] = signals[i]->evaluate(pCell);
        }
        return signal_values;
    }

    void add_signal(std::unique_ptr<PhysiCell::AbstractSignal> pSignal)
    {
        signals.push_back(std::move(pSignal));
    }

    AggregatorSignal()
    {
        aggregator = multivariate_hill_aggregator;
    }

    AggregatorSignal(std::vector<std::unique_ptr<PhysiCell::AbstractSignal>> pSignals) : AggregatorSignal()
    {
        signals = std::move(pSignals);
    }
};

class MediatorSignal : public AbstractAggregatorSignal
{
private:
    std::unique_ptr<PhysiCell::AbstractSignal>decreasing_signal;
    std::unique_ptr<PhysiCell::AbstractSignal>increasing_signal;

    // values that the aggregator can use to calculate the output
    double min_value = 0.1;
    double base_value = 1;
    double max_value = 10;

public:
    std::vector<double> get_signals(Cell *pCell)
    {
        return {decreasing_signal->evaluate(pCell), increasing_signal->evaluate(pCell)};
    }

    double decreasing_dominant_mediator(std::vector<double> signals_in);
    double increasing_dominant_mediator(std::vector<double> signals_in);
    double neutral_mediator(std::vector<double> signals_in);

    MediatorSignal()
    {
        decreasing_signal = std::unique_ptr<PhysiCell::AbstractSignal>(new AggregatorSignal());
        increasing_signal = std::unique_ptr<PhysiCell::AbstractSignal>(new AggregatorSignal());
        aggregator = [this](std::vector<double> signals_in)
        {
            return this->decreasing_dominant_mediator(signals_in);
        };
    }

    MediatorSignal(double val) : MediatorSignal()
    {
        min_value = val;
        base_value = val;
        max_value = val;
    }

    MediatorSignal(std::unique_ptr<PhysiCell::AbstractSignal> pDecreasingSignal, std::unique_ptr<PhysiCell::AbstractSignal> pIncreasingSignal, double min = 0.1, double base = 1, double max = 10)
        : decreasing_signal(std::move(pDecreasingSignal)), increasing_signal(std::move(pIncreasingSignal)), min_value(min), base_value(base), max_value(max)
    {
        if (decreasing_signal == nullptr)
        {
            throw std::invalid_argument("Null pointer passed to MediatorSignal constructor for decreasing signal");
        }
        if (increasing_signal == nullptr)
        {
            throw std::invalid_argument("Null pointer passed to MediatorSignal constructor for increasing signal");
        }
        aggregator = [this](std::vector<double> signals_in)
        { return this->decreasing_dominant_mediator(signals_in); };
    }

    AbstractSignal* get_decreasing_signal() { return decreasing_signal.get(); }
    AbstractSignal* get_increasing_signal() { return increasing_signal.get(); }

    void set_decreasing_signal(std::unique_ptr<PhysiCell::AbstractSignal> pSignal) { decreasing_signal = std::move(pSignal); }
    void set_increasing_signal(std::unique_ptr<PhysiCell::AbstractSignal> pSignal) { increasing_signal = std::move(pSignal); }
    
    void set_min_value(double val) { min_value = val; }
    void set_base_value(double val) { base_value = val; }
    void set_max_value(double val) { max_value = val; }

    void use_increasing_dominant_mediator();
    void use_neutral_mediator();
};

class ElementarySignal : public AbstractSignal
{
private:

protected:
    std::string signal_name;
    bool applies_to_dead_cells;

public:

    virtual double transformer( double signal )
    {
        return signal;
    }

    ElementarySignal(std::string signal_name, bool applies_to_dead_cells) : signal_name(signal_name), applies_to_dead_cells(applies_to_dead_cells) {}

    virtual ~ElementarySignal() {}
};

class SignalReference
{
private:

protected:
    double reference_value = 0;

public:
    virtual double coordinate_transform(double signal) = 0;

    SignalReference(double reference_value_) : reference_value(reference_value_) {}

    virtual ~SignalReference() {}
};

class IncreasingSignalReference : public SignalReference
{
public:
    double coordinate_transform(double signal)
    {
        if (signal <= reference_value)
        {
            return 0;
        }
        return signal - reference_value;
    }

    using SignalReference::SignalReference; // Inherit constructors from SignalReference
    
    IncreasingSignalReference();

    IncreasingSignalReference(double reference_value_) : SignalReference(reference_value_) {}
};

class DecreasingSignalReference : public SignalReference
{
public:
    double coordinate_transform(double signal)
    {
        if (signal >= reference_value)
        {
            return 0;
        }
        return reference_value - signal;
    }

    using SignalReference::SignalReference; // Inherit constructors from SignalReference
    
    DecreasingSignalReference();

    DecreasingSignalReference(double reference_value_) : SignalReference(reference_value_) {}
};

class RelativeSignal : public ElementarySignal
{
protected:
    std::unique_ptr<SignalReference> signal_reference = nullptr;

public:
    virtual void add_reference(std::unique_ptr<SignalReference> pSR) = 0;

    double evaluate(Cell *pCell);

    using ElementarySignal::ElementarySignal; // Inherit constructors from ElementarySignal
};

class AbstractHillSignal : public RelativeSignal
{
private:
    double half_max;
    double hill_power;

public:
    void add_reference(std::unique_ptr<SignalReference> pSR)
    {
        signal_reference = std::move(pSR);
        half_max = signal_reference->coordinate_transform(half_max);
        if (half_max <= 0)
        {
            throw std::invalid_argument("Half max must be in the range in which the signal has effect.");
        }
    }

    double rescale(double signal)
    {
        signal /= half_max;
        return pow(signal, hill_power);
    }

    AbstractHillSignal(std::string signal_name, bool applies_to_dead_cells, double half_max, double hill_power)
        : RelativeSignal(signal_name, applies_to_dead_cells), half_max(half_max), hill_power(hill_power) {}
};

class PartialHillSignal : public AbstractHillSignal
{
public:
    double transformer(double signal)
    {
        return rescale(signal);
    }

    PartialHillSignal(std::string signal_name, bool applies_to_dead_cells, double half_max, double hill_power)
        : AbstractHillSignal(signal_name, applies_to_dead_cells, half_max, hill_power) {}
};

class HillSignal : public AbstractHillSignal
{
public:
    double transformer(double signal)
    {
        signal = rescale(signal);
        signal /= 1+signal;
        return signal;
    }

    HillSignal(std::string signal_name, bool applies_to_dead_cells, double half_max, double hill_power)
        : AbstractHillSignal(signal_name, applies_to_dead_cells, half_max, hill_power) {}
};

// AbsoluteSignals do not have need of a reference value
class AbsoluteSignal : public ElementarySignal
{
protected:

public:
    double evaluate(Cell *pCell);
    using ElementarySignal::ElementarySignal; // Inherit constructors from ElementarySignal
};

class LinearSignal : public AbsoluteSignal
{
private:

protected:
    double signal_min;
    double signal_max;
    double signal_range;

public:
    LinearSignal(std::string signal_name, bool applies_to_dead_cells, double signal_min, double signal_max)
        : AbsoluteSignal(signal_name, applies_to_dead_cells), signal_min(signal_min), signal_max(signal_max)
    {
        signal_range = signal_max - signal_min;
    }
};

class IncreasingLinearSignal : public LinearSignal
{
    double transformer(double signal)
    {
        if (signal <= signal_min) {
            return 0;
        }
        else if (signal >= signal_max) {
            return 1;
        }
        else {
            return (signal - signal_min) / signal_range;
        }
    }

    using LinearSignal::LinearSignal;
};

class DecreasingLinearSignal : public LinearSignal
{
    double transformer(double signal)
    {
        if (signal <= signal_min) {
            return 1;
        }
        else if (signal >= signal_max) {
            return 0;
        }
        else {
            return (signal_max - signal) / signal_range;
        }
    }

    using LinearSignal::LinearSignal;
};

class HeavisideSignal : public AbsoluteSignal
{
private:

protected:
    double threshold;

public:

    HeavisideSignal(std::string signal_name, bool applies_to_dead_cells, double threshold)
        : AbsoluteSignal(signal_name, applies_to_dead_cells), threshold(threshold) {}
};

class IncreasingHeavisideSignal : public HeavisideSignal
{
public:
    double transformer(double signal)
    {
        if (signal < threshold)
        {
            return 0;
        }
        else
        {
            return 1;
        }
    }

    IncreasingHeavisideSignal(std::string signal_name, bool applies_to_dead_cells, double threshold)
        : HeavisideSignal(signal_name, applies_to_dead_cells, threshold) {}
};

class DecreasingHeavisideSignal : public HeavisideSignal
{
public:
    double transformer(double signal)
    {
        if (signal > threshold)
        {
            return 0;
        }
        else
        {
            return 1;
        }
    }

    DecreasingHeavisideSignal(std::string signal_name, bool applies_to_dead_cells, double threshold)
        : HeavisideSignal(signal_name, applies_to_dead_cells, threshold) {}
};

class BehaviorRule
{
private:

protected:

public:
    std::string behavior;
    std::unique_ptr<AbstractSignal> signal;

    virtual void apply(Cell* pCell) = 0;

    BehaviorRule(std::string behavior)
        : behavior(behavior)
    {
        signal = std::unique_ptr<AbstractSignal>(new MediatorSignal());
    }

    BehaviorRule(std::string behavior, std::unique_ptr<AbstractSignal> pSignal)
        : behavior(behavior), signal(std::move(pSignal)) {}

    virtual ~BehaviorRule() {}

    // std::string cell_type;
    // Cell_Definition *pCell_Definition;

    // BehaviorRule();
    // BehaviorRule(std::string cell_type, std::string behavior, double min_behavior, double max_behavior);


    // void sync_to_cell_definition(Cell_Definition *pCD);


    // void add_signal(std::string, std::string response);


    /* to do
    void reduced_display(std::ostream &os);  // done
    void display(std::ostream &os);          // done
    void detailed_display(std::ostream &os); // done

    void English_display(std::ostream &os);
    void English_display_HTML(std::ostream &os);
    void English_detailed_display(std::ostream &os);
    void English_detailed_display_HTML(std::ostream &os);
    */
};

class BehaviorSetter : public BehaviorRule
{
public:
    void apply(Cell* pCell);
    using BehaviorRule::BehaviorRule; // Inherit constructors from BehaviorRule
};

class BehaviorAccumulator : public BehaviorRule
{
private:
    double base_value;
    double saturation_value;

public:
    void apply(Cell* pCell);

    BehaviorAccumulator(std::string behavior, std::unique_ptr<AbstractSignal> pSignal, double base_value, double saturation_value)
        : BehaviorRule(behavior, std::move(pSignal)), base_value(base_value), saturation_value(saturation_value) 
    {
        if (saturation_value < base_value)
        {
            throw std::invalid_argument("Saturation value must be greater than or equal to base value for an accumulator.");
        }
    }
};

class BehaviorAttenuator : public BehaviorRule
{
private:
    double base_value;
    double saturation_value;

public:
    void apply(Cell* pCell);

    BehaviorAttenuator(std::string behavior, std::unique_ptr<AbstractSignal> pSignal, double base_value, double saturation_value)
        : BehaviorRule(behavior, std::move(pSignal)), base_value(base_value), saturation_value(saturation_value) 
    {
        if (saturation_value > base_value)
        {
            throw std::invalid_argument("Saturation value must be less than or equal to base value for an attenuator.");
        }
    }
};

class BehaviorRuleset
{
private:
    std::vector<std::unique_ptr<BehaviorRule>> rules;

public:
    void apply(Cell* pCell);
    void add_behavior_rule(std::unique_ptr<BehaviorRule> pRule) {
        rules.push_back(std::move(pRule));
    }

    BehaviorRuleset() {}

    BehaviorRule *find_behavior(std::string behavior);

    // std::string cell_type;
    // Cell_Definition *pCell_Definition;

    // BehaviorRuleset();

    // void add_behavior(BehaviorRule *pRule) {
    //     rules.push_back(pRule);
    // }

    // void add_behavior(std::string behavior, double min_behavior, double max_behavior);

    // void apply(Cell* pCell) {
    //     for (auto& rule : rules) {
    //         rule->apply(pCell);
    //     }
    // }

    // void sync_to_cell_definition(Cell_Definition *pCD) {
    //     pCell_Definition = pCD;
    //     sync_to_cell_definition_finish(pCD->name);
    // }

    // void sync_to_cell_definition(std::string cell_type);

    // void sync_to_cell_definition_finish(std::string name)
    // {
    //     cell_type = name;
    //     for (auto &rule : rules)
    //     {
    //         rule->sync_to_cell_definition(pCell_Definition);
    //     }
    // };
};

void apply_behavior_ruleset( Cell* pCell );

void parse_xml_behavior_rules(const std::string filename);
std::unique_ptr<BehaviorRule> parse_behavior(std::string cell_type, std::string behavior, pugi::xml_node node);

bool signal_is_mediator(pugi::xml_node parent_node);
bool signal_is_aggregator(pugi::xml_node parent_node);
bool signal_is_elementary(pugi::xml_node parent_node);

std::unique_ptr<PhysiCell::AbstractSignal>parse_abstract_signal(pugi::xml_node parent_node);
std::unique_ptr<PhysiCell::AbstractSignal>parse_mediator_signal(pugi::xml_node mediator_node);
std::unique_ptr<PhysiCell::AbstractSignal>parse_aggregator_signal(pugi::xml_node aggregator_node);
std::unique_ptr<PhysiCell::AbstractSignal>parse_elementary_signal(pugi::xml_node elementary_node);

std::unique_ptr<AbstractSignal> parse_hill_signal(std::string name, pugi::xml_node hill_node);
std::unique_ptr<AbstractSignal> parse_partial_hill_signal(std::string name, pugi::xml_node partial_hill_node);
std::unique_ptr<AbstractSignal> parse_linear_signal(std::string name, pugi::xml_node linear_node);
std::unique_ptr<AbstractSignal> parse_heaviside_signal(std::string name, pugi::xml_node heaviside_node);

void parse_reference(pugi::xml_node reference_node, RelativeSignal *pRelSig);

void setup_behavior_rules( void );
void parse_behavior_rules_from_pugixml( void );
void parse_behavior_rules_from_file(std::string path_to_file, std::string format, std::string protocol, double version); // see PhysiCell_rules.h for default values of format, protocol, and version

double euler_direct_solve(double current, double rate, double target);
double exponential_solve(double current, double rate, double target);

// loading rules from CSV files
void parse_csv_behavior_rules(std::string path_to_file, std::string protocol, double version);
void parse_csv_behavior_rule(std::string line);
bool is_csv_rule_misformed(std::vector<std::string> input);
void split_csv( std::string input , std::vector<std::string>& output , char delim );
std::string csv_strings_to_English_v3( std::vector<std::string> strings , bool include_cell_header );

};

#endif
