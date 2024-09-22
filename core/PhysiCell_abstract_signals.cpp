#include "./PhysiCell_abstract_signals.h"

double signal_sum( std::vector<double> signals_in )
{
    std::cout << "Summing signals" << std::endl;
    double signal_sum = 0;
    for (auto& signal : signals_in) {
        std::cout << "\tSignal: " << signal << std::endl;
        signal_sum += signal;
    }
    return signal_sum;
}

double multivariate_hill_aggregator( std::vector<double> partial_hill_signals )
{
    std::cout << "Inside multivariate_hill_aggregator with " << partial_hill_signals.size() << " signals" << std::endl;
    for (size_t i = 0; i < partial_hill_signals.size(); ++i )
    {
        std::cout << "partial_hill_signals[i]: " << partial_hill_signals[i] << std::endl;
    }
    double out = signal_sum(partial_hill_signals);
    return out / (1 + out);
}

double signal_difference( std::vector<double> signals_in )
{
    std::cout << "Mediating signals" << std::endl;
    std::cout << "\tSignal 1: " << signals_in[0] << std::endl;
    std::cout << "\tSignal 2: " << signals_in[1] << std::endl << std::endl;
    return signals_in[1] - signals_in[0];
}

double get_single_signal(Cell* pC, std::string signal_name)
{
    std::cout << "inside get_single_signal" << std::endl;
    if (signal_name == "pressure")
    {return pC->pressure;}
    else if (signal_name == "oxygen" )
    { return pC->oxygen;}
    else
    { return 0; }
}


int main (int argc, char* argv[])
{
    BehaviorRule rule;
    double min_value = 0.3;
    double base_value = 0.5;
    double max_value = 0.9;

    PartialHillSignal pressure_signal;
    pressure_signal.signal_name = "pressure";
    pressure_signal.half_max = 0.5;
    pressure_signal.hill_power = 1;

    PartialHillSignal low_oxygen_signal;
    low_oxygen_signal.signal_name = "oxygen";
    low_oxygen_signal.half_max = 3;
    low_oxygen_signal.hill_power = 1;


    DecreasingSignalReference oxygen_signal_reference = DecreasingSignalReference(10);
    low_oxygen_signal.add_reference(&oxygen_signal_reference);

    // AggregatorSignal aggregator;
    // aggregator.signal_aggregator = signal_sum;

    // LinearSignal signal;
    // signal.signal_min = 0.0;
    // signal.signal_max = 1.0;

    // aggregator.signals.push_back(&signal);

    // LinearSignal signal2;
    // signal2.signal_min = 0.5;
    // signal2.signal_max = 1.0;

    // HillSignal signal3;
    // signal3.half_max = 0.5;
    // signal3.hill_power = 2.0;

    AggregatorSignal decreasing_signals;
    decreasing_signals.append_signal(&pressure_signal);
    decreasing_signals.append_signal(&low_oxygen_signal);
    decreasing_signals.aggregator = multivariate_hill_aggregator;

    AggregatorSignal increasing_signals;
    increasing_signals.aggregator = multivariate_hill_aggregator;
    MediatorSignal mediator(&decreasing_signals, &increasing_signals, min_value, base_value, max_value);

    // use lambda function to set the aggregator
    // mediator.aggregator = [](std::vector<double> signals_in) -> double {
    //     return decreasing_dominant_mediator(signals_in, 0.3, 0.5, 0.9);
    // };
    // aggregator.signals.push_back(&signal2);
    // mediator

    // AggregatorSignal aggregator;
    // aggregator.append_signal(&mediator);
    // aggregator.append_signal(&signal3); 

    rule.signal = &mediator;

    std::cout << "APPLYING" << std::endl;

    Cell cell;
    rule.apply(&cell);
}