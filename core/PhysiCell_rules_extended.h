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
/** This function returns the first signal from the input vector. */
double first_aggregator(std::vector<double> signals_in);
/** This function returns the sum of all signals in the input vector. */
double sum_aggregator(std::vector<double> signals_in);
/** This function returns the multivariate hill aggregation of the input signals. */
double multivariate_hill_aggregator(std::vector<double> partial_hill_signals);
/** This function returns the product of all signals in the input vector. */
double product_aggregator(std::vector<double> signals_in);
/** This function returns the mean of all signals in the input vector. */
double mean_aggregator(std::vector<double> signals_in);
/** This function returns the maximum signal from the input vector. */
double max_aggregator(std::vector<double> signals_in);
/** This function returns the minimum signal from the input vector. */
double min_aggregator(std::vector<double> signals_in);
/** This function returns the median signal from the input vector. */
double median_aggregator(std::vector<double> signals_in);
/** This function returns the geometric mean of all signals in the input vector. */
double geometric_mean_aggregator(std::vector<double> signals_in);

/** This class serves as the base for all signal types. */
class AbstractSignal
{
public:
    virtual double evaluate(Cell *pCell) = 0;

    virtual ~AbstractSignal() {}
};

/** 
 * @brief Abstract Class for all signals that are aggregators of other signals.
 * 
 * Combines 1+ signals into a single signal using a specified aggregation function.
 */
class AbstractAggregatorSignal : public AbstractSignal
{
private:
    virtual std::vector<double> evaluate_signals(Cell *pCell) = 0;
    
public:
    std::function<double(std::vector<double>)> aggregator;

    double evaluate(Cell *pCell) override
    {
        return aggregator(evaluate_signals(pCell));
    }

    virtual ~AbstractAggregatorSignal() {}
};

/** 
 * @brief Concrete class that (typically) handles signals in the same direction (e.g., all increasing).
 */
class AggregatorSignal : public AbstractAggregatorSignal
{
private:
    std::vector<std::unique_ptr<PhysiCell::AbstractSignal>> signals;

public:
    std::vector<double> evaluate_signals(Cell *pCell) override
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

    void use_custom_aggregator( void )
    {
        aggregator = [](std::vector<double> signals_in)
        {
            std::cerr << "ERROR: Custom aggregator not set! Make sure to set one in custom.cpp if using a custom aggregator." << std::endl;
            exit(-1);
            return -1;
        };
    }
};

/** 
 * @brief Concrete class that specifically mediates between decreasing and increasing signals.
 * 
 * This class is used to combine two signals, one decreasing and one increasing, into a single signal.
 * The two input signals can be of any type, but they are typically AggregatorSignals.
 * 
 * The class also includes slots for the minimum, base, and maximum values of the output signal.
 * These are automatically used in the class methods to calculate the output signal.
 */
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
    std::vector<double> evaluate_signals(Cell *pCell) override
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

    double get_min_value() const { return min_value; }
    double get_base_value() const { return base_value; }
    double get_max_value() const { return max_value; }

    void use_increasing_dominant_mediator();
    void use_neutral_mediator();
    void use_custom_mediator()
    {
        aggregator = [](std::vector<double> signals_in)
        {
            std::cerr << "ERROR: Custom mediator not set! Make sure to set one in custom.cpp if using a custom mediator." << std::endl;
            exit(-1);
            return -1;
        };
    }
};

/** 
 * @brief Abstract Class for all signals that come from a single, measured signal in PhysiCell.
 * 
 * Instances of this class must have a single signal that will determine how they read the PhysiCell simulation.
 * A transformer function is used to transform the signal before outputting it to the next layer up the signal chain.
 */
class ElementarySignal : public AbstractSignal
{
private:

protected:
    std::string signal_name;
    bool applies_to_dead;

public:
    virtual double evaluate(Cell *pCell) override = 0;

    virtual double transformer(double signal) = 0;

    ElementarySignal(std::string signal_name, bool applies_to_dead) : signal_name(signal_name), applies_to_dead(applies_to_dead) {}

    virtual ~ElementarySignal() {}
};

/**
 * @brief Concrete class that sets a reference for a RelativeSignal.
 * 
 * The reference is applied to the signal before passing it to the transformer function.
 */
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

/**
 * @brief Concrete class for a reference that sets the minimum value for a signal to have effect.
 * 
 * If a RelativeSignal reads a signal that is less than the reference value, it will be set to 0.
 * Otherwise, the amount above the reference value is passed to the transformer function.
 */
class IncreasingSignalReference : public SignalReference
{
public:
    double coordinate_transform(double signal) override
    {
        if (signal <= reference_value)
        {
            return 0;
        }
        return signal - reference_value;
    }

    using SignalReference::SignalReference; // Inherit constructors from SignalReference
};

/**
 * @brief Concrete class for a reference that sets the maximum value for a signal to have effect.
 * 
 * If a RelativeSignal reads a signal that is greater than the reference value, it will be set to 0.
 * Otherwise, the amount below the reference value is passed to the transformer function.
 */
class DecreasingSignalReference : public SignalReference
{
public:
    double coordinate_transform(double signal) override
    {
        if (signal >= reference_value)
        {
            return 0;
        }
        return reference_value - signal;
    }

    using SignalReference::SignalReference; // Inherit constructors from SignalReference
};

/**
 * @brief Abstract Class for all signals that are relative to a reference signal.
 * 
 * This class is used to create signals that are relative to a reference signal.
 * The reference signal is set using a SignalReference object, which is passed to the add_reference function.
 * The add_reference method will update internals of the signal to use the reference value.
 */
class RelativeSignal : public ElementarySignal
{
protected:
    std::unique_ptr<SignalReference> signal_reference = nullptr;

public:
    virtual double transformer(double signal) override = 0;
    virtual void add_reference(std::unique_ptr<SignalReference> pSR) = 0;

    double evaluate(Cell *pCell) override;

    using ElementarySignal::ElementarySignal; // Inherit constructors from ElementarySignal

    virtual ~RelativeSignal() {}
};

/**
 * @brief Abstract Class for all signals that are based on a Hill function.
 * 
 * Holds the two parameters of the Hill function: half_max and hill_power.
 * If a non-trivial reference is set, the half_max value is updated to be relative to the reference.
 * This avoids repeating calculations in the transformer function.
 */
class AbstractHillSignal : public RelativeSignal
{
private:
    double half_max;
    double hill_power;

public:
    virtual double transformer(double signal) override = 0;
    void add_reference(std::unique_ptr<SignalReference> pSR) override
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

    AbstractHillSignal(std::string signal_name, bool applies_to_dead, double half_max, double hill_power)
        : RelativeSignal(signal_name, applies_to_dead), half_max(half_max), hill_power(hill_power) {}
};

/**
 * @brief Concrete class that is a partial Hill function signal.
 * 
 * The partial Hill function is defined as (signal / half_max)^hill_power.
 * This is typically used to compute the multivariate Hill function.
 * Importantly, the output of this signal is on [0, ∞).
 */
class PartialHillSignal : public AbstractHillSignal
{
public:
    double transformer(double signal) override
    {
        return rescale(signal);
    }

    PartialHillSignal(std::string signal_name, bool applies_to_dead, double half_max, double hill_power)
        : AbstractHillSignal(signal_name, applies_to_dead, half_max, hill_power) {}
};

/**
 * @brief Concrete class that is a Hill function signal.
 * 
 * The Hill function is defined as x / (1 + x) where x is (signal / half_max)^hill_power.
 * The output of this signal is on [0, 1).
 */
class HillSignal : public AbstractHillSignal
{
public:
    double transformer(double signal) override
    {
        signal = rescale(signal);
        signal /= 1 + signal;
        return signal;
    }

    HillSignal(std::string signal_name, bool applies_to_dead, double half_max, double hill_power)
        : AbstractHillSignal(signal_name, applies_to_dead, half_max, hill_power) {}
};

/**
 * @brief Concrete class that is an identity signal.
 * 
 * This signal does not transform the input signal.
 * An alternative would be a PartialHillSignal with a half_max of 1 and a hill_power of 1.
 * A reference value can be set.
 */
class IdentitySignal : public RelativeSignal
{
public:
    double transformer(double signal) override
    {
        return signal;
    }

    void add_reference(std::unique_ptr<SignalReference> pSR) override
    {
        signal_reference = std::move(pSR);
    }

    IdentitySignal(std::string signal_name, bool applies_to_dead)
        : RelativeSignal(signal_name, applies_to_dead) {}
};

// AbsoluteSignals do not have need of a reference value
/**
 * @brief Abstract Class for all signals that are absolute.
 * 
 * Absolute signals cannot have a reference value.
 */
class AbsoluteSignal : public ElementarySignal
{
protected:

public:
    double evaluate(Cell *pCell) override;
    virtual double transformer(double signal) override = 0;
    using ElementarySignal::ElementarySignal; // Inherit constructors from ElementarySignal

    virtual ~AbsoluteSignal() {}
};

/**
 * @brief Abstract Class for all signals that are linear.
 * 
 * Linear signals are defined by a minimum and maximum signal.
 * Between these two values, the signal changes linearly between 0 and 1.
 * Outside of this range, the signal is set to 0 or 1.
 */
class LinearSignal : public AbsoluteSignal
{
private:

protected:
    double signal_min;
    double signal_max;
    double signal_range;

public:
    virtual double transformer(double signal) override = 0;
    LinearSignal(std::string signal_name, bool applies_to_dead, double signal_min, double signal_max)
        : AbsoluteSignal(signal_name, applies_to_dead), signal_min(signal_min), signal_max(signal_max)
    {
        signal_range = signal_max - signal_min;
    }

    virtual ~LinearSignal() {}
};

/**
 * @brief Concrete class that is an increasing linear signal.
 * 
 * Increases from 0 (at signal_min) to 1 (at signal_max).
 * The increasing linear signal is defined as (signal - signal_min) / (signal_max - signal_min).
 */
class IncreasingLinearSignal : public LinearSignal
{
    double transformer(double signal) override
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

/**
 * @brief Concrete class that is a decreasing linear signal.
 * 
 * Decreases from 1 (at signal_min) to 0 (at signal_max).
 * The decreasing linear signal is defined as (signal_max - signal) / (signal_max - signal_min).
 */
class DecreasingLinearSignal : public LinearSignal
{
    double transformer(double signal) override
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

/**
 * @brief Abstract Class for all signals that are Heaviside functions.
 * 
 * Heaviside functions are defined by a threshold value.
 * If the signal is at or beyond the threshold, the output is 1.
 * Otherwise, the output is 0.
 */
class HeavisideSignal : public AbsoluteSignal
{
private:

protected:
    double threshold;

public:
    virtual double transformer(double signal) override = 0;

    HeavisideSignal(std::string signal_name, bool applies_to_dead, double threshold)
        : AbsoluteSignal(signal_name, applies_to_dead), threshold(threshold) {}

    virtual ~HeavisideSignal() {}
};

/**
 * @brief Concrete class that is an increasing Heaviside signal.
 * 
 * The increasing Heaviside signal is 0 below the threshold and 1 at or above the threshold.
 */
class IncreasingHeavisideSignal : public HeavisideSignal
{
public:
    double transformer(double signal) override
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

    IncreasingHeavisideSignal(std::string signal_name, bool applies_to_dead, double threshold)
        : HeavisideSignal(signal_name, applies_to_dead, threshold) {}
};

/**
 * @brief Concrete class that is a decreasing Heaviside signal.
 * 
 * The decreasing Heaviside signal is 0 above the threshold and 1 at or below the threshold.
 */
class DecreasingHeavisideSignal : public HeavisideSignal
{
public:
    double transformer(double signal) override
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

    DecreasingHeavisideSignal(std::string signal_name, bool applies_to_dead, double threshold)
        : HeavisideSignal(signal_name, applies_to_dead, threshold) {}
};

/**
 * @brief Abstract Class for all behavior rules.
 * 
 * This class sets the behavior targeted and the signal that is used to set it.
 * Typically, this signal is a MediatorSignal, but it can be any signal type.
 */
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

/**
 * @brief Concrete class that sets a behavior for a cell.
 * 
 * This class is the standard way to set a behavior for a cell.
 * After evaluating the signals, the behavior is set to the value of the signal.
 */
class BehaviorSetter : public BehaviorRule
{
public:
    void apply(Cell* pCell) override;
    using BehaviorRule::BehaviorRule; // Inherit constructors from BehaviorRule
};

/**
 * @brief Abstract class for behavior rules that gradually change behaviors toward saturation values.
 * 
 * This class serves as the base for both BehaviorAccumulator and BehaviorAttenuator classes,
 * which implement behaviors that change according to ODE dynamics rather than being set directly.
 * 
 * The two derived classes are BehaviorAccumulator and BehaviorAttenuator.
 * Both classes allow a behavior to increase or decrease over time, respectively, so long as the input signal is sufficiently strong.
 * If the input signal is not strong enough, the behavior will relax back to the base value.
 * 
 * The base_value field is the value that the behavior will return to if the input signal is not strong enough.
 * The saturation_limit field is the value that the behavior will approach if the input signal is strong enough.
 * 
 * The output of the signal is the rate of change of the behavior.
 * If the rate is positive, the behavior will approach the saturation value.
 * If the rate is negative, the behavior will approach the base value.
 * Thus, users should make sure that the output of the signal can be positive or negative.
 */
class BehaviorRateSetter : public BehaviorRule
{
protected:
    double base_value;
    double saturation_limit;

public:
    virtual void apply(Cell* pCell) override = 0;
    
    BehaviorRateSetter(std::string behavior, std::unique_ptr<AbstractSignal> pSignal, double base_value, double saturation_limit)
        : BehaviorRule(behavior, std::move(pSignal)), base_value(base_value), saturation_limit(saturation_limit) {}
    
    virtual ~BehaviorRateSetter() {}
};

/**
 * @brief Concrete class that accumulates a behavior over time.
 * 
 * This class allows for a behavior to increase over time if the input signal is sufficiently strong.
 * Otherwise, it will relax back to the base value.
 * 
 * See BehaviorRateSetter for more details.
 */
class BehaviorAccumulator : public BehaviorRateSetter
{
public:
    void apply(Cell* pCell) override;

    BehaviorAccumulator(std::string behavior, std::unique_ptr<AbstractSignal> pSignal, double base_value, double saturation_limit)
        : BehaviorRateSetter(behavior, std::move(pSignal), base_value, saturation_limit) 
    {
        if (saturation_limit < base_value)
        {
            throw std::invalid_argument("Saturation value must be greater than or equal to base value for an accumulator.");
        }
    }
};

/**
 * @brief Concrete class that attenuates a behavior over time.
 * 
 * This class allows for a behavior to decrease over time if the input signal is sufficiently strong.
 * Otherwise, it will relax back to the saturation value.
 * 
 * See BehaviorRateSetter for more details.
 */
class BehaviorAttenuator : public BehaviorRateSetter
{
public:
    void apply(Cell* pCell) override;

    BehaviorAttenuator(std::string behavior, std::unique_ptr<AbstractSignal> pSignal, double base_value, double saturation_limit)
        : BehaviorRateSetter(behavior, std::move(pSignal), base_value, saturation_limit) 
    {
        if (saturation_limit > base_value)
        {
            throw std::invalid_argument("Saturation value must be less than or equal to base value for an attenuator.");
        }
    }
};

/**
 * @brief Class that holds a set of behavior rules for a cell type.
 * 
 * This class is used to apply behavior rules to a cell.
 * It holds a vector of BehaviorRule objects and applies them to the cell.
 */
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
std::unique_ptr<AbstractSignal> parse_identity_signal(std::string name, pugi::xml_node identity_node);
std::unique_ptr<AbstractSignal> parse_linear_signal(std::string name, pugi::xml_node linear_node);
std::unique_ptr<AbstractSignal> parse_heaviside_signal(std::string name, pugi::xml_node heaviside_node);

void parse_reference(pugi::xml_node reference_node, RelativeSignal *pRelSig);

void setup_behavior_rules( void );
BehaviorRuleset* find_behavior_ruleset( Cell_Definition* pCD );
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

void set_custom_mediator(const std::string &cell_definition_name, const std::string &behavior_name, double (*mediator_function)(std::vector<double>));
void set_custom_mediator(const std::string &cell_definition_name, const std::string &behavior_name, double (*mediator_function)(MediatorSignal*, std::vector<double>));
void set_custom_aggregator(const std::string &cell_definition_name, const std::string &behavior_name, const std::string &response, double (*aggregator_function)(std::vector<double>));
MediatorSignal* get_top_level_mediator(const std::string &cell_definition_name, const std::string &behavior_name);
};

#endif
