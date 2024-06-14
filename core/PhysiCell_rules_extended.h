#include <vector>
#include <string>
#include <functional>


#ifndef __PhysiCell_rules_extended__
#define __PhysiCell_rules_extended__

#include "./PhysiCell.h"

namespace PhysiCell{

double default_aggregator( std::vector<double> signals_in ) {
    double signal_sum = 0;
    for (auto& signal : signals_in) {
        signal_sum += signal;
    }
    return signal_sum;
};

class AbstractSignal
{
public:
    virtual double evaluate( Cell* pCell ) = 0;

    virtual ~AbstractSignal() {}
};

class AbstractAggregatorSignal : public AbstractSignal
{
public:
    virtual std::vector<double> get_signals( Cell* pCell ) = 0;
    std::function<double(std::vector<double>)> aggregator;

    double evaluate( Cell* pCell ) {
        return aggregator(get_signals(pCell));
    }

    virtual ~AbstractAggregatorSignal() {}
};

class AggregatorSignal : public AbstractAggregatorSignal
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

    AggregatorSignal() {
        aggregator = default_aggregator;
    }
};

class MediatorSignal : public AbstractAggregatorSignal
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

    double default_mediator(std::vector<double> signals_in)
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
    };

    double aggregate_signals( Cell* pCell ) {
        std::vector<double> signal_values = get_signals(pCell);
        return aggregator(signal_values);
    }

    MediatorSignal() {
        aggregator = [this](std::vector<double> signals_in) {
            return this->default_mediator(signals_in);
        };
    }
    MediatorSignal(AbstractSignal *pDecreasingSignal, AbstractSignal *pIncreasingSignal) : MediatorSignal()
    {
        if (pDecreasingSignal == nullptr || pIncreasingSignal == nullptr)
        {
            throw std::invalid_argument("Null pointer passed to MediatorSignal constructor");
        }
        decreasing_signal = pDecreasingSignal;
        increasing_signal = pIncreasingSignal;
    }
    MediatorSignal(AbstractSignal *pDecreasingSignal, AbstractSignal *pIncreasingSignal, double min, double base, double max) : MediatorSignal(pDecreasingSignal, pIncreasingSignal)
    {
        min_value = min;
        base_value = base;
        max_value = max;
    }
};

class ElementarySignal : public AbstractSignal
{
public:
    std::string signal_name;
    bool applies_to_dead_cells;

    virtual double transformer( double signal ) {
        return signal;
    }

    double evaluate( Cell* pCell );

    virtual ~ElementarySignal() {}
};

class AbstractHillSignal : public ElementarySignal
{
public:
    double half_max;
    double hill_power;

    double rescale(double signal)
    {
        signal /= half_max;
        return pow(signal, hill_power);
    }
};

class PartialHillSignal : public AbstractHillSignal
{
public:
    double half_max;
    double hill_power;

    double transformer(double signal)
    {
        std::cout << "Evaluating Partial Hill signal" << std::endl;
        signal = rescale(signal);
        std::cout << "\tReturning " << signal << "\n" << std::endl;
        return signal;
        // signal /= half_max;
        // signal = pow(signal, hill_power);
        // return signal;
    }
};

class HillSignal : public AbstractHillSignal
{
public:
    double half_max;
    double hill_power;

    double transformer(double signal)
    {
        std::cout << "Evaluating Hill signal" << std::endl;
        signal = rescale(signal);
        signal /= 1+signal;
        std::cout << "\tReturning " << signal << "\n" << std::endl;
        return signal;
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


class BehaviorRule
{
public:
    std::string cell_type;
    Cell_Definition *pCell_Definition;

    std::string behavior;
    double min_value = 0.1;
    double base_value = 1.0;
    double max_value = 10.0;

    BehaviorRule();
    BehaviorRule(std::string cell_type, std::string behavior, double min_behavior, double max_behavior);


    void sync_to_cell_definition(Cell_Definition *pCD);

    AbstractSignal* signal;

    void add_signal(std::string, std::string response);

    double evaluate(Cell* pCell) {
        double return_value = signal->evaluate(pCell);
        std::cout << "Cell behavior set to " << return_value << std::endl;
        return return_value;
    }

    void apply(Cell* pCell);

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

public:
    std::vector<BehaviorRule *> rules;
    std::string cell_type;
    Cell_Definition *pCell_Definition;

    BehaviorRuleset();

    void add_behavior(BehaviorRule *pRule) {
        rules.push_back(pRule);
    }

    void add_behavior(std::string behavior, double min_behavior, double max_behavior);

    void apply(Cell* pCell) {
        for (auto& rule : rules) {
            rule->apply(pCell);
        }
    }

    void sync_to_cell_definition(Cell_Definition *pCD) {
        pCell_Definition = pCD;
        sync_to_cell_definition_finish(pCD->name);
    }

    void sync_to_cell_definition(std::string cell_type);

    void sync_to_cell_definition_finish(std::string name)
    {
        cell_type = name;
        for (auto &rule : rules)
        {
            rule->sync_to_cell_definition(pCell_Definition);
        }
    };
};

void parse_xml_rules_extended(std::string filename);
BehaviorRule *parse_abstract_signal(const pugi::xml_node parent_node);
BehaviorRule *parse_mediator_signal(const pugi::xml_node parent_node);
BehaviorRule *parse_aggregator_signal(const pugi::xml_node parent_node);
BehaviorRule *parse_elementary_signal(const pugi::xml_node parent_node);
bool signal_is_mediator(pugi::xml_node parent_node);
bool signal_is_aggregator(pugi::xml_node parent_node);
bool signal_is_elementary(pugi::xml_node parent_node);
}; 

#endif 
