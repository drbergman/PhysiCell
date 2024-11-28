#include <vector>
#include <string>
#include <functional>
#include <memory>


#ifndef __PhysiCell_rules_extended__
#define __PhysiCell_rules_extended__

#include "./PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 

namespace PhysiCell{

double sum_aggregator(std::vector<double> signals_in);
double multivariate_hill_aggregator(std::vector<double> partial_hill_signals);

class AbstractSignal
{
public:
    virtual double evaluate(Cell *pCell) = 0;

    virtual ~AbstractSignal() {}
};

class AbstractAggregatorSignal : public AbstractSignal
{
public:
    virtual std::vector<double> get_signals(Cell *pCell) = 0;
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

    void append_signal(std::unique_ptr<PhysiCell::AbstractSignal> pSignal)
    {
        signals.push_back(std::move(pSignal));
    }

    AggregatorSignal()
    {
        aggregator = sum_aggregator;
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

public:
    // values that the aggregator can use to calculate the output
    double min_value = 0.1;
    double base_value = 1;
    double max_value = 10;

    std::vector<double> get_signals(Cell *pCell)
    {
        return {increasing_signal->evaluate(pCell), decreasing_signal->evaluate(pCell)};
    }

    double decreasing_dominant_mediator(std::vector<double> signals_in);
    double increasing_dominant_mediator(std::vector<double> signals_in);
    double neutral_mediator(std::vector<double> signals_in);

    MediatorSignal()
    {
        aggregator = [this](std::vector<double> signals_in)
        {
            return this->decreasing_dominant_mediator(signals_in);
        };
    }

    MediatorSignal(std::unique_ptr<PhysiCell::AbstractSignal> pDecreasingSignal, std::unique_ptr<PhysiCell::AbstractSignal> pIncreasingSignal, double min = 0.1, double base = 1, double max = 10)
        : decreasing_signal(std::move(pDecreasingSignal)), increasing_signal(std::move(pIncreasingSignal)), min_value(min), base_value(base), max_value(max)
    {
        if (decreasing_signal == nullptr || increasing_signal == nullptr)
        {
            if (decreasing_signal == nullptr)
            {
                throw std::invalid_argument("Null pointer passed to MediatorSignal constructor for decreasing signal");
            }
            if (increasing_signal == nullptr)
            {
                throw std::invalid_argument("Null pointer passed to MediatorSignal constructor for increasing signal");
            }
        }
        aggregator = [this](std::vector<double> signals_in)
        { return this->decreasing_dominant_mediator(signals_in); };
    }

    void use_increasing_dominant_mediator();
    void use_neutral_mediator();
};

class SignalReference
{
public:
    double reference_value = 0;
    virtual double coordinate_transform(double signal) = 0;
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
    IncreasingSignalReference();

    IncreasingSignalReference(double reference_value_)
    {
        reference_value = reference_value_;
    }
};
class DecreasingSignalReference : public SignalReference
{
public:
    double coordinate_transform(double signal)
    {
        std::cout << "Inside coord transform for DecreasingSignalReference with signal = " << signal << std::endl;
        std::cout << "\treference_value = " << reference_value << std::endl;
        if (signal >= reference_value)
        {
            return 0;
        }
        std::cout << "returning " << reference_value - signal << std::endl;
        return reference_value - signal;
    }
    DecreasingSignalReference();

    DecreasingSignalReference(double reference_value_)
    {
        reference_value = reference_value_;
    }
};

class ElementarySignal : public AbstractSignal
{
private:
    std::unique_ptr<SignalReference> signal_reference = nullptr;

public:
    std::string signal_name;
    bool applies_to_dead_cells;

    virtual double transformer( double signal )
    {
        return signal;
    }

    std::function<void(SignalReference*)> add_signal_reference;

    void add_reference(std::unique_ptr<SignalReference> pSR);
    double evaluate(Cell *pCell);

    ElementarySignal(std::string signal_name, bool applies_to_dead_cells) : signal_name(signal_name), applies_to_dead_cells(applies_to_dead_cells) {}

    virtual ~ElementarySignal() {}
};

class AbstractHillSignal : public ElementarySignal
{
protected:
    void abstract_hill_signal_initializer()
    {
        add_signal_reference = [this](SignalReference *pSR)
        {
            half_max = pSR->coordinate_transform(half_max);
            if (half_max <= 0)
            {
                throw std::invalid_argument("Half max must be in the range in which the signal has effect.");
            }
        };
    }

public:
    double half_max;
    double hill_power;

    AbstractHillSignal(std::string signal_name, bool applies_to_dead_cells, double half_max, double hill_power)
        : ElementarySignal(signal_name, applies_to_dead_cells), half_max(half_max), hill_power(hill_power)
    {
        abstract_hill_signal_initializer();
    }

    double rescale(double signal)
    {
        signal /= half_max;
        return pow(signal, hill_power);
    }
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

class LinearSignal : public ElementarySignal
{
public:
    double signal_min;
    double signal_max;

    double transformer(double signal)
    {
        if (signal <= signal_min) {
            return 0;
        }
        else if (signal >= signal_max) {
            return 1;
        }
        else {
            return (signal - signal_min) / (signal_max - signal_min);
        }
    }

    LinearSignal(std::string signal_name, bool applies_to_dead_cells, double signal_min, double signal_max)
        : ElementarySignal(signal_name, applies_to_dead_cells), signal_min(signal_min), signal_max(signal_max)
    {
        add_signal_reference = [this](SignalReference *pSR) {
            double bound_1 = pSR->coordinate_transform(this->signal_min);
            double bound_2 = pSR->coordinate_transform(this->signal_max);
            if (bound_1 <= bound_2)
            {
                this->signal_min = bound_1;
                this->signal_max = bound_2;
            }
            else // the signal reference swapped the order of the bounds
            {
                this->signal_min = bound_2;
                this->signal_max = bound_1;
            }
        };
    }
};

class HeavisideSignal : public ElementarySignal
{
public:
    double threshold;

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

    HeavisideSignal(std::string signal_name, bool applies_to_dead_cells, double threshold)
        : ElementarySignal(signal_name, applies_to_dead_cells), threshold(threshold)
    {
        add_signal_reference = [this](SignalReference *pSR)
        {
            this->threshold = pSR->coordinate_transform(this->threshold);
        };
    }
};

class BehaviorRule
{
public:
    std::string behavior;
    std::unique_ptr<AbstractSignal> signal;

    void apply(Cell* pCell);

    BehaviorRule(std::string behavior, std::unique_ptr<AbstractSignal> pSignal)
        : behavior(behavior), signal(std::move(pSignal)) {}

    // std::string cell_type;
    // Cell_Definition *pCell_Definition;

    // BehaviorRule();
    // BehaviorRule(std::string cell_type, std::string behavior, double min_behavior, double max_behavior);


    // void sync_to_cell_definition(Cell_Definition *pCD);


    // void append_signal(std::string, std::string response);


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

void parse_xml_rules_extended(const std::string filename);
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

void parse_reference(pugi::xml_node reference_node, ElementarySignal *pES);

void setup_behavior_rules( void );
void parse_behavior_rules_from_pugixml( void );
void parse_behavior_rules_from_file(std::string path_to_file, std::string format, std::string protocol, double version); // see PhysiCell_rules.h for default values of format, protocol, and version

}; 

#endif
