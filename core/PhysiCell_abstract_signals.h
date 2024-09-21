#include <vector>
#include <string>
#include <iostream>
#include <functional>


double signal_sum( std::vector<double> signals_in );
double signal_difference( std::vector<double> signals_in );
// double decreasing_dominant_mediator(std::vector<double> signals_in, double min_value, double base_value, double max_value);

class Cell
{
public:
    double signal = 0.6;
    bool is_dead = false;
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

    double aggregate_signals( Cell* pCell ) {
        std::vector<double> signal_values = get_signals(pCell);
        return aggregator(signal_values);
    }

    MediatorSignal() {
        aggregator = [this](std::vector<double> signals_in) {
            return this->decreasing_dominant_mediator(signals_in);
        };
    }
    MediatorSignal(AbstractSignal *pIncreasingSignal, AbstractSignal *pDecreasingSignal)
    {
        if (pIncreasingSignal == nullptr || pDecreasingSignal == nullptr)
        {
            throw std::invalid_argument("Null pointer passed to MediatorSignal constructor");
        }
        std::cout << "Creating mediator" << std::endl;
        increasing_signal = pIncreasingSignal;
        decreasing_signal = pDecreasingSignal;
        aggregator = [this](std::vector<double> signals_in)
        { return this->decreasing_dominant_mediator(signals_in); };
    }
    MediatorSignal(AbstractSignal *pIncreasingSignal, AbstractSignal *pDecreasingSignal, double min, double base, double max)
    {
        if (pIncreasingSignal == nullptr || pDecreasingSignal == nullptr)
        {
            throw std::invalid_argument("Null pointer passed to MediatorSignal constructor");
        }
        std::cout << "Creating mediator" << std::endl;
        increasing_signal = pIncreasingSignal;
        decreasing_signal = pDecreasingSignal;
        min_value = min;
        base_value = base;
        max_value = max;
        aggregator = [this](std::vector<double> signals_in)
        { return this->decreasing_dominant_mediator(signals_in); };
    }
};

class ElementarySignal : public AbstractSignal
{
public:
    std::string signal;
    bool applies_to_dead_cells;

    virtual double transformer( double signal ) {
        return signal;
    }

    double evaluate( Cell* pCell ) {
        if (!applies_to_dead_cells && pCell->is_dead)
        {
            return 0;
        }
        return transformer(pCell->signal);
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