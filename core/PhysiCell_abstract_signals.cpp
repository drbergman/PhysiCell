#include "./PhysiCell_abstract_signals.h"

double sum_fn( std::vector<double> signals_in )
{
    std::cout << "Summing signals" << std::endl;
    double signal_sum = 0;
    for (auto& signal : signals_in) {
        std::cout << "\tSignal: " << signal << std::endl;
        signal_sum += signal;
    }
    return signal_sum;
}

double diff_fn( std::vector<double> signals_in )
{
    std::cout << "Mediating signals" << std::endl;
    std::cout << "\tSignal 1: " << signals_in[0] << std::endl;
    std::cout << "\tSignal 2: " << signals_in[1] << std::endl << std::endl;
    return signals_in[1] - signals_in[0];
}

double default_mediator_aggregator(std::vector<double> signals_in, double min_value, double base_value, double max_value)
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
    // std::cout << "Mediating signals" << std::endl;
    // std::cout << "\tSignal 1: " << signals_in[0] << std::endl;
    // std::cout << "\tSignal 2: " << signals_in[1] << std::endl << std::endl;
    // return min_value * signals_in[0] + (base_value + (max_value - base_value) * signals_in[1]) * (1 - signals_in[0]);
}

int main (int argc, char* argv[])
{
    BehaviorRule rule;
    rule.min_value = 0.3;
    rule.base_value = 0.5;
    rule.max_value = 0.9;

    // SignalAggregator aggregator;
    // aggregator.signal_aggregator = sum_fn;

    LinearSignal signal;
    signal.signal_min = 0.0;
    signal.signal_max = 1.0;

    // aggregator.signals.push_back(&signal);

    LinearSignal signal2;
    signal2.signal_min = 0.5;
    signal2.signal_max = 1.0;

    HillSignal signal3;
    signal3.half_max = 0.5;
    signal3.hill_power = 2.0;

    SignalMediator mediator(&signal, &signal2);
    // use lambda function to set the aggregator
    mediator.aggregator = [](std::vector<double> signals_in) -> double {
        return default_mediator_aggregator(signals_in, 0.3, 0.5, 0.9);
    };
    // aggregator.signals.push_back(&signal2);

    SignalAggregator aggregator;
    aggregator.add_signal(&mediator);
    aggregator.add_signal(&signal3); 

    rule.signal = &mediator;

    Cell cell;
    rule.apply(&cell);
}