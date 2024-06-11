#include <vector>
#include <string>
#include <iostream>

double sum_fn( std::vector<double> signals_in );
double diff_fn( std::vector<double> signals_in );
double default_mediator_aggregator(std::vector<double> signals_in, double min_value, double base_value, double max_value);

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

class AbstractSignalAggregator : public AbstractSignal
{
public:
    virtual std::vector<double> get_signals( Cell* pCell ) = 0;
    double (*aggregator)(std::vector<double> signals_in) = 0;

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
    std::vector<double> get_signals( Cell* pCell ) {
        return {increasing_signal->evaluate(pCell), decreasing_signal->evaluate(pCell)};
    }

    double aggregate_signals( Cell* pCell ) {
        std::vector<double> signal_values = get_signals(pCell);
        return aggregator(signal_values);
    }

    SignalMediator() {
        aggregator = diff_fn;
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
        aggregator = diff_fn;
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
        double signal_value = signal->evaluate(pCell);
        double return_value;
        if (signal_value < 0) {
            return_value = min_value * signal_value + base_value * (1 - signal_value);
        }
        else if (signal_value > 0) {
            return_value = max_value * signal_value + base_value * (1 - signal_value);
        }
        else {
            return_value = base_value;
        }
        std::cout << "Cell behavior set to " << return_value << std::endl;
    }
};