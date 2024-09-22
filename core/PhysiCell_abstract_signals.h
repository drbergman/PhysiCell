#include <vector>
#include <string>
#include <iostream>
#include <functional>


double signal_sum( std::vector<double> signals_in );
double multivariate_hill_aggregator( std::vector<double> partial_hill_signals );
double signal_difference( std::vector<double> signals_in );
// double decreasing_dominant_mediator(std::vector<double> signals_in, double min_value, double base_value, double max_value);
class Cell
{
public:
    double pressure = 0.5;
    double oxygen = 7.0;
    bool is_dead = false;
};
double get_single_signal(Cell* pC, std::string signal_name);

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
        std::cout << "Inside evaluate..." << std::endl;
        std::vector<double> signal_values = get_signals(pCell);
        std::cout << "Inside evaluate with " << signal_values.size() << " signals." << std::endl;
        return aggregator(signal_values);
    }
};

class AggregatorSignal : public AbstractAggregatorSignal
{
private:
    std::vector<AbstractSignal *> signals;

public:
    std::vector<double> get_signals( Cell* pCell ) {
        std::vector<double> signal_values;
        signal_values.resize(signals.size());
        for (size_t i = 0; i < signals.size(); ++i) {
            signal_values[i] = signals[i]->evaluate(pCell);
            std::cout << "signal_values[i]: " << signal_values[i] << std::endl;
        }
        std::cout << "Returning from get_signals inside AggregatorSignal with " << signal_values.size() << " signals." << std::endl;
        return signal_values;
    }

    void append_signal(AbstractSignal *pSignal) {
        signals.push_back(pSignal);
    }

    AggregatorSignal() {
        aggregator = signal_sum;
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

    double decreasing_dominant_mediator(std::vector<double> signals_in)
    {
        std::cout << "Default mediator aggregator" << std::endl;
        std::cout << "\tSignal 1: " << signals_in[0] << ", Signal 2: " << signals_in[1] << ", Min: " << min_value << ", Base: " << base_value << ", Max: " << max_value << std::endl;
        double out = base_value;
        std::cout << "\tBase value: " << out << std::endl;
        out += (max_value - base_value) * signals_in[1];
        std::cout << "\tAfter adding signal 2: " << out << std::endl;
        out *= 1 - signals_in[0];
        std::cout << "\tAfter multiplying by signal 1: " << out << std::endl;
        out += min_value * signals_in[0];
        std::cout << "\tAfter adding signal 1: " << out << std::endl;
        return out;
        // return min_value * signals_in[0] + (base_value + (max_value - base_value) * signals_in[1]) * (1 - signals_in[0]);
    };

    MediatorSignal() {
        aggregator = [this](std::vector<double> signals_in) {
            return this->decreasing_dominant_mediator(signals_in);
        };
    }

    MediatorSignal(AbstractSignal *pIncreasingSignal, AbstractSignal *pDecreasingSignal, double min = 0.1, double base = 1, double max = 10)
        : increasing_signal(pIncreasingSignal), decreasing_signal(pDecreasingSignal), min_value(min), base_value(base), max_value(max)
    {
        if (pIncreasingSignal == nullptr || pDecreasingSignal == nullptr)
        {
            throw std::invalid_argument("Null pointer passed to MediatorSignal constructor");
        }
        std::cout << "Creating mediator" << std::endl;
        aggregator = [this](std::vector<double> signals_in)
        { return this->decreasing_dominant_mediator(signals_in); };
    }
};

class SignalReference
{
public:
    double reference_value = 0;
    virtual double coordinate_transform(double signal) = 0;
    virtual ~SignalReference() {}
};

class IncreasingSignalReference : SignalReference
{
public: 
    double coordinate_transform(double signal) {
        if (signal <= reference_value)
        { return 0;}
        return signal - reference_value;
    }
    IncreasingSignalReference();

    IncreasingSignalReference(double reference_value_) {
        reference_value = reference_value_;
    }
};
class DecreasingSignalReference : public SignalReference
{
public: 
    double coordinate_transform(double signal) {
        std::cout << "Inside coord transform for DecreasingSignalReference with signal = " << signal << std::endl;
        std::cout << "\treference_value = " << reference_value << std::endl;
        if (signal >= reference_value)
        { return 0;}
        std::cout << "returning " << reference_value - signal << std::endl;
        return reference_value - signal;
    }
    DecreasingSignalReference();

    DecreasingSignalReference(double reference_value_) {
        reference_value = reference_value_;
    }
};

class ElementarySignal : public AbstractSignal
{
public:
    std::string signal_name;
    bool applies_to_dead_cells;
    SignalReference* signal_reference = nullptr;

    virtual double transformer( double signal ) {
        return signal;
    }

    std::function<void(SignalReference*)> add_signal_reference;

    void add_reference(SignalReference* pSR) {
        signal_reference = pSR;
        add_signal_reference(pSR);
    }

    double evaluate( Cell* pCell ) {
        if (!applies_to_dead_cells && pCell->is_dead)
        {
            return 0;
        }
        std::cout << "applying to this cell" << std::endl;
        double signal = get_single_signal(pCell, signal_name);
        std::cout << "Signal: " << signal << std::endl;
        std::cout << "signal_reference for " << signal_name << ": " << signal_reference << std::endl;
        if (signal_reference!=nullptr) {
            signal = signal_reference->coordinate_transform(signal);
        }
        return transformer(signal);
    }

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

    HillSignal()
    {
        add_signal_reference = [this](SignalReference *pSR) {
            this->half_max = pSR->coordinate_transform(this->half_max);
        };
    }
};

class PartialHillSignal : public ElementarySignal
{
public:
    double half_max;
    double hill_power;

    double transformer(double signal)
    {
        std::cout << "Evaluating Partial Hill signal" << std::endl;
        signal /= half_max;
        signal = pow(signal, hill_power);
        std::cout << "\tReturning " << signal << "\n" << std::endl;
        return signal;
    }

    PartialHillSignal()
    {
        add_signal_reference = [this](SignalReference *pSR) {
            double new_half_max = pSR->coordinate_transform(this->half_max);
            std::cout << "New half max: " << new_half_max << std::endl;
            this->half_max = new_half_max;
        };
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
        if (signal <= signal_min) {
            std::cout << "\tReturning 0\n" << std::endl;
            return 0;
        }
        else if (signal >= signal_max) {
            std::cout << "\tReturning 1\n" << std::endl;
            return 1;
        }
        else {
            std::cout << "\tReturning " << (signal - signal_min) / (signal_max - signal_min) << "\n" << std::endl;
            return (signal - signal_min) / (signal_max - signal_min);
        }
    }

    LinearSignal()
    {
        add_signal_reference = [this](SignalReference *pSR) {
            double bound_1 = pSR->coordinate_transform(this->signal_min);
            double bound_2 = pSR->coordinate_transform(this->signal_max);
            if (bound_1 <= bound_2)
            {
                this->signal_min = bound_1;
                this->signal_max = bound_2;
            }
            else
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
        if (signal < threshold) {
            return 0;
        }
        else {
            return 1;
        }
    }

    HeavisideSignal()
    {
        add_signal_reference = [this](SignalReference *pSR) {
            this->threshold = pSR->coordinate_transform(this->threshold);
        };
    }
};

class BehaviorRule
{
public:
    AbstractSignal* signal;

    void apply(Cell* pCell) {
        double return_value = signal->evaluate(pCell);
        std::cout << "Cell behavior set to " << return_value << std::endl;
        return;
    }
};