# Introduction
This PR implements an extension of the rules grammar using XMLs. In this more robust format, the rules are extended in many ways (see below). **Importantly, this PR still allows for CSVs to be read in.** In fact, all sample projects in this PR are updated to use this new framework while still using CSVs for their rules. See the `template_xml_rules_extended` sample project for examples exploring the breadth of this implementation.

A core design principle of this implementation is to make the original rules the default settings. That is, if `type`s, `aggregator`s, and `mediator`s are omitted, the original rules are used.

# Table of Contents
- [Introduction](#introduction)
- [Table of Contents](#table-of-contents)
- [Key Features](#key-features)
  - [New signal `transformer`s](#new-signal-transformers)
    - [XML Examples](#xml-examples)
  - [Signal reference values](#signal-reference-values)
    - [XML Examples](#xml-examples-1)
  - [New signal `aggregator`s](#new-signal-aggregators)
    - [XML Examples](#xml-examples-2)
  - [New signal `mediator`s](#new-signal-mediators)
    - [XML Examples](#xml-examples-3)
  - [Behavior `setter`s, `accumulator`s, and `attenuator`s](#behavior-setters-accumulators-and-attenuators)
    - [Accumulators and attenuators](#accumulators-and-attenuators)
    - [XML Examples](#xml-examples-4)
  - [Parsed output](#parsed-output)
- [Useful examples](#useful-examples)
  - [Example 1: Switch behaviors](#example-1-switch-behaviors)
  - [Example 2: Signal X AND Signal Y](#example-2-signal-x-and-signal-y)
  - [Example 3: LOW Signal X (from $x\_0$)](#example-3-low-signal-x-from-x_0)
  - [Example 4: Hysteresis](#example-4-hysteresis)
- [Advanced features](#advanced-features)
- [Automatic conversion of formats](#automatic-conversion-of-formats)
  - [From CSV to XML](#from-csv-to-xml)
  - [From XML to CSV](#from-xml-to-csv)
- [To-dos](#to-dos)
- [Reviewer notes](#reviewer-notes)

# Key Features
## New signal `transformer`s

The elementary signals of PhysiCell can first be transformed before being passed on to the next layer. By default, the signals are passed as partial Hill functions (see Table below). This expands on this by introducing the `ElementarySignal` class and the `transformer` method. In a rules XML, these would go into a `<signal>` element.
See below for examples. In the Table below, the `Name` indicates the attribute value.

| Name | Description | Parameters | Formula |
|-------|-------------|-----------|---------|
| `partial_hill` | default | `half_max` ($\gamma$), `hill_power` ($n$) | $\left(\frac{s}{\gamma}\right)^n$ |
| `hill` | Hill function | `half_max` ($\gamma$), `hill_power` ($n$) | $\frac{x}{1+x},\ x=\left(\frac{s}{\gamma}\right)^n$ |
| `linear` | linearly vary on $[0,1]$ on an range of signal values, extending constantly and continuously outisde the interval | `signal_min` ($s_1$), `signal_max` ($s_2$) | $\frac{s-s_1}{s_2-s_1}$ or $\frac{s_2-s}{s_2-s_1}$ on $[s_1,s_2]$ |
| `heaviside` | 0 or 1 depending on the value of the signal relative to the threshold | `threshold` ($T$) | $H(s-T)$ or $H(T-s)$|
| `identity` | identity function | NA | $s$ |

### XML Examples
If the `type` attribute is omitted below (the one specifying `partial_hill`), the default is `partial_hill`.
```xml
<signal name="pressure" type="partial_hill">
    <half_max>0.5</half_max>
    <hill_power>2</hill_power>
    <applies_to_dead>0</applies_to_dead>
</signal>
```
If the `type` element is omitted below (the one specifying `increasing`), the default is `increasing`.
```xml
<signal name="oxygen" type="linear">
    <type>increasing</type>
    <signal_min>0</signal_min>
    <signal_max>10</signal_max>
    <applies_to_dead>0</applies_to_dead>
</signal>
```

## Signal reference values
For cases in which the signal is in reference to a fixed value, a `<reference>` element can be added to the signal element. References can be either increasing or decreasing. The difference with the reference value is passed on to the transformer. If increasing, this differences is $s-s_0$. If decreasing, this difference is $s_0-s$. Whenever the difference is negative, 0 is passed on to the transformer.

Only a subset of transformers accept reference values: `partial_hill`, `hill`, and `identity`. In the former two cases, the `half_max` must be in the support of the reference, i.e. above an increasing reference and below a decreasing reference. The `half_max` is adjusted accordingly in this case (see below).

### XML Examples
In the first example, the signal is in reference to a fixed value of 0.25. Thus, the signal is passed on to the transformer as $s-0.25$. **Importantly**, the `half_max` is adjusted accordingly as well so that the transformation is

$$
\left(\frac{s - 0.25}{0.5 - 0.25}\right)^2
$$

In this way, the transformer still behaves as expected: at the `half_max`, the transformed signal is 0.5. This transformation to the `half_max` is done one time when the reference is set.

```xml
<signal name="contact with dead" type="partial_hill">
    <half_max>0.5</half_max>
    <hill_power>2</hill_power>
    <reference>
        <type>increasing</type>
        <value>0.25</value>
    </reference>
    <applies_to_dead>1</applies_to_dead>
</signal>
```

In the second example, *decreasing* `pressure` from 1.5 is passed on to the transformer.
```xml
<signal name="pressure" type="hill">
    <half_max>0.5</half_max>
    <hill_power>2</hill_power>
    <reference>
        <type>decreasing</type>
        <value>1.5</value>
    </reference>
    <applies_to_dead>0</applies_to_dead>
</signal>
```

If the `pressure` signal is 1.0, the Hill (not partial Hill) transformer outputs

$$
\left(\frac{1.5-1.0}{1.5-0.5}\right)^2 = 0.25 \to \frac{0.25}{1+0.25} = 0.2
$$

## New signal `aggregator`s
An `AggregatorSignal` is used to combine multiple signals into one. These perform the function of the increasing signals (or decreasing signals) in the original rules. The default for these is the `multivariate_hill` aggregator (see below). Every `AggregatorSignal` has a `aggregator` method that takes the vector of doubles from its constituent signals and returns a double.

In the Table below, $x_i$ is one of the constituent (post-transformation) signals into the aggregator.

| Name | Description | Formula |
|-------|-------------|---------|
| `multivariate_hill` | default | $\frac{\sum_i x_i}{1+\sum_i x_i}$ |
| `sum` | sum of signals | $\sum_i x_i$ |
| `product` | product of signals | $\prod_i x_i$ |
| `min` | minimum of signals | $\min_i x_i$ |
| `max` | maximum of signals | $\max_i x_i$ |
| `mean` | mean of signals | $\frac{1}{n}\sum_i x_i$ |
| `median` | median of signals | $\text{median}(x_1,x_2,\ldots,x_n)$ |
| `geometric_mean` | geometric mean of signals | $\sqrt[n]{\prod_i x_i}$ |
| `first` | first signal (avoids overhead if only one signal) | $x_1$ |
| `custom` | custom aggregator (set in `custom.cpp`) | `double f(std::vector<double> signals_in);` |

Note: if using the `custom` aggregator, see the example in the `template_xml_rules_extended` project.

### XML Examples
The most common place to find aggregators is in the `<increasing_signals>` and `<decreasing_signals>` elements. See [Advanced features](#advanced-features) for other ways to use aggregators.

```xml
<increasing_signals>
    <aggregator>sum</aggregator>
    <signal name="oxygen" type="linear">
        ...
    </signal>
    ...more signals...
</increasing_signals>
```

## New signal `mediator`s
A `MediatorSignal` is a signal that is used to mediate between two or more signals. These perform the function of combining the increasing and decreasing signals in the original rules. The default for these is the `decreasing_dominant` mediator (see below). Every `MediatorSignal` has an `aggregator` method (they are a subclass of `AbstractAggregator`) that takes in the decreasing signal and the increasing signal (in that order) and returns a double.

In the Table below, the min behavior is $b_m$, base is $b_0$, and max is $b_M$. $D$ is the decreasing signal and $U$ is the increasing signal ($U$ for up). If omitted, the behavior values default to $0.1$, $1.0$, and $10.0$, respectively.

| Name | Description | Formula |
|-------|-------------|---------|
| `decreasing_dominant` | default | $Db_m + (1-D)((1-U)b_0 + Ub_M)$ |
| `increasing_dominant` | increasing signals can override the decreasing signals | $Ub_M + (1-U)((1-D)b_0 + Db_m)$ |
| `neutral` | neutral mediator | $b_0 + D(b_m-b_0) + U(b_M-b_0)$ |
| `custom` | custom mediator (set in `custom.cpp`) | `double f(MediatorSignal* pMS, std::vector<double> signals_in);` or `double f(std::vector<double> signals_in);` |

Note: the first custom mediator signature includes a `MediatorSignal*` pointer to allow the custom mediator to also make use of the `min`, `base`, and `max` values in the `MediatorSignal` class. If using the `custom` mediator, see the example in the `template_xml_rules_extended` project.

### XML Examples
The most common place to find mediators is implicitly in the `<behavior>` element.
```xml
<behavior name="migration speed">
    <mediator>increasing_dominant</mediator>
    <increasing_signals>
        ...signals...
    </increasing_signals>
    <decreasing_signals>
        ...signals...
    </decreasing_signals>
</behavior>
```

## Behavior `setter`s, `accumulator`s, and `attenuator`s
Behaviors can be a `setter` (default) in which they set the behavior value of the cell. They can now also be `accumulator`s or `attenuator`s that cause the behavior to increase or decrease, respectively, over time.

### Accumulators and attenuators
When a behavior is an `accumulator` or `attenuator`, the behavior value is not set directly. Instead, the behavior value is updated according to the rate of change of the behavior. The rate of change is determined by the signals. The output of the top level signal (typically a `MediatorSignal`) is used as the rate of change of the behavior. If the rate is positive, the behavior moves away from the base value. If the rate is negative, the behavior moves towards the base value. **It is therefore very important to make sure that the rate values include positive and negative values.**

To explain, let $r_-<0$ be the minimum rate (what the decreasing signals drive the rate towards) and let $r_+>0$ be the maximum rate (what the increasing signals drive the rate towards). Also, let $r_0$ be the base rate (the rate in the absence of any signals). $r_0$ will often be negative indicating that without stimulus, the behavior relaxes to the base value. The user must further set a saturation limit for the behavior, $b_s$. Note: none of these values are part of PhysiCell prior to this PR and none are directly observable in the output data since they are only part of the rule, not the cell. Finally, the behavior does have a base value of $b_0$ that is part of PhysiCell, e.g. the base migration speed. This can be read from the configuration file or set in the `<behavior_base>` element.

Then, when the behavior is calculated, the top level mediator will balance the decreasing and increasing signals to settle on a rate, $r$, in the range $r_- \le r_0 \le r_+$. The behavior is then updated as follows:

$$
b' = \begin{cases}
|r|(b_0-b), & r<0 \\
0, & r=0 \\
|r|(b_s-b), & r>0
\end{cases}
$$

In other words, the behavior is exponentially "decaying" towards the equilibrium value, which is determined by the sign of the rate.

### XML Examples
```xml
<behavior name="migration speed">
    <type>accumulator</type>
    <base_value note="r₀">-0.01</base_value>
    <behavior_base note="b₀">0.0</behavior_base>
    <behavior_saturation note="bₛ">5.0</behavior_saturation>
    <mediator>increasing_dominant</mediator>
    <decreasing_signals>
        <max_response note="r₋">-0.1</max_response>
        <aggregator>first</aggregator>
        <signal name="pressure" type="linear">
            <signal_min>0.2</signal_min>
            <signal_max>1.0</signal_max>
            <applies_to_dead>0</applies_to_dead>
        </signal>
    </decreasing_signals>
    <increasing_signals>
        <max_response note="r₊">0.2</max_response>
        <aggregator>first</aggregator>
        <signal name="oxygen" type="heaviside">
            <threshold>12.0</threshold>
            <applies_to_dead>0</applies_to_dead>
        </signal>
    </increasing_signals>
</behavior>
```

## Parsed output
PhysiCell will create a CSV file from the parsed XML rules file and output it to `<output_folder>/cell_rules_parsed.csv`. See the [`template_xml_rules_extended`](https://github.com/drbergman/PhysiCell/blob/1.14.2-drbergman-2.1.2/sample_projects/template_xml_rules_extended/config/cell_rules_parsed.csv) sample project for an example of the output.

# Useful examples
## Example 1: Switch behaviors
This is literally the purpose of including the `heaviside` transformer, but it is still worth noting that we can now do this very elegantly.
```xml
<signal name="pressure" type="heaviside">
    <threshold>0.5</threshold>
    <applies_to_dead>0</applies_to_dead>
</signal>
```

## Example 2: Signal X AND Signal Y
A limitation of the previous framework was that two increasing signals were combined through an OR operator, i.e., only one of X or Y needed to be high for the behavior to be increased. Now, an AND conjunction is straightforward to implement. A future goal is to make convert this into a syntax that can easily be read by the user. For now, the following is a valid XML snippet:
```xml
<behavior name="migration speed">
    <base_value>0.0</base_value>
    <increasing_signals>
        <max_response>2.0</max_response>
        <aggregator>product</aggregator>
        <signal name="oxygen" type="linear">
            <signal_min>5</signal_min>
            <signal_max>15</signal_max>
            <applies_to_dead>0</applies_to_dead>
        </signal>
        <signal name="pressure" type="linear">
            <type>decreasing</type>
            <signal_min>0</signal_min>
            <signal_max>2.0</signal_max>
            <applies_to_dead>0</applies_to_dead>
        </signal>
    </increasing_signals>
</behavior>
```

The `migration speed` can be increased up to the `max_response` of 2.0 if `oxygen` is at or above 15 AND pressure is at 0. If either oxygen is below 5 or pressure is above 2.0, the `migration speed` is set to 0.0. Otherwise, the product of the two linearly-transformed signals will be on $[0,1]$; this value is then multiplied by the difference `max_response - base_value` and added to the `base_value` to get the final `migration speed`.

## Example 3: LOW Signal X (from $x_0$)
The ability to have a signal decrease result in a stronger effect is now possible, opening up rules statements like "LOW oxygen increases transform to mesenchymal". In fact, the above example showed such an example, but this is worth highlighting on its own.
```xml
<behavior name="transform to mesenchymal">
    <base_value>0.0</base_value>
    <increasing_signals>
        <aggregator>first</aggregator>
        <max_response>0.1</max_response>
        <signal name="oxygen" type="hill">
            <half_max>8.0</half_max>
            <hill_power>4.0</hill_power>
            <applies_to_dead>0</applies_to_dead>
            <reference>
                <type>decreasing</type>
                <value>10.0</value>
            </reference>
        </signal>
    </increasing_signals>
</behavior>
```

Note that in this example, the `transform to mesenchymal` behavior is set to 0.0 when `oxygen` is above 10.0. The `half_max` of 8.0 means that the behavior will be 50% towards its max value when `oxygen` is at 8.0, i.e., it will be at 0.05. When `oxygen` is 0.0, then the behavior will be $0.1\cdot(5^4/(1+5^4)) \approx 0.09984$.

## Example 4: Hysteresis
Oftentimes, we want a signal to initiate a process that eventaully leads to change in the cell behavior. This lag between observed signal and subsequent behavior can be modeled using the `accumulator` and `attenuator` behavior types. For example, the following XML snippet shows how to implement a hysteresis effect:
```xml
<behavior name="damage">
    <type>accumulator</type>
    <base_value>-0.01</base_value>
    <behavior_base>0.0</behavior_base>
    <behavior_saturation>1.0</behavior_saturation>
    <mediator>neutral</mediator>
    <increasing_signals>
        <max_response>0.2</max_response>
        <aggregator>first</aggregator>
        <signal name="chemotherapy" type="heaviside">
            <threshold>0.5</threshold>
            <applies_to_dead>0</applies_to_dead>
        </signal>
    </increasing_signals>
    <decreasing_signals>
        <max_response>-0.1</max_response>
        <aggregator>first</aggregator>
        <signal name="custom:repair_factor" type="hill">
            <half_max>1.5</half_max>
            <hill_power>2.0</hill_power>
            <applies_to_dead>0</applies_to_dead>
        </signal>
    </decreasing_signals>
</behavior>
```

In this example, the substrate `chemotherapy` induces the accumulation of damage to the cell.
So long as the extracellular `chemotherapy` concentration is at or above 0.5, this signal will increase the damage accumulation rate by 0.2 - (-0.01) = 0.21. The damage will accumulate towards the saturation limit of 1.0. If the `chemotherapy` concentration drops below 0.5, the damage accumulation rate will revert back to the `base_value` defined here (-0.01). Since it is negative, that means the cell damage will exponentially "decay" towards the cell-type-defined base damage value, which in this case is set explicitly in the rules to 0.0 by the `<behavior_base>` element.

If the cell has `custom:repair_factor` (perhaps defined by an intracellular model), this will decrease the damage accumulation rate towards -0.1. At the `half_max` value, 50% of the decrease in rate will be achieved so that the rate will be decreased by $0.5\cdot(-0.1-(-0.01))=-0.045$.

Since a `neutral` mediator is used, the increasing and decreasing deltas will be computed separately and simply added together. Thus, at `chemotherapy` $\ge0.5$ and `custom:repair_factor` $\gg1.5$, the damage accumulation rate will be $-0.01 + (-0.045) + 0.2 = 0.145$.

Note: PhysiCell has a mechanism for this already for repairing `damage`. However, this is unique to `damage` and no other behavior has this. This new framework allows for any behavior to have this hysteresis effect.

# Advanced features
The expected way for users to take advantage of this new framework is to use the previous rules structure:
- mediate between increasing and decreasing elementary signals using the `decreasing_dominant` mediator
- aggregate increasing signals using the `multivariate_hill` aggregator, similarly for decreasing signals
- use the `partial_hill` transformer for all (elementary) signals

Users can simply swap out the mediators, aggregators, and transformers to the new ones presented here.

Some users may also find value in the reference values and in the accumulator/attenuator features. These are conceptually new and, in time, the community will develop terminology to describe them.

This framework does go deeper. Users can make any hierarchy of signals they like. In the custom code, the user can craft any behavior they want.

In the XML, the top two levels for a given rule must be the mediator and then the increasing/decreasing aggregators. After that, anything goes. For example, the following is a valid XML snippet:
```xml
<behavior name="migration speed">
    <increasing_signals>
        <signal type="aggregator">
            <aggregator>product</aggregator>
            <signal name="oxygen" type="linear">
                <signal_min>0</signal_min>
                <signal_max>10</signal_max>
                <applies_to_dead>0</applies_to_dead>
            </signal>
            <signal name="glucose" type="linear">
                <signal_min>0</signal_min>
                <signal_max>10</signal_max>
                <applies_to_dead>0</applies_to_dead>
            </signal>
            <signal type="mediator">
                <mediator>neutral</mediator>
                <base_value>0.1</base_value>
                <decreasing_signals>
                    <max_response>0.01</max_response>
                    <signal name="pressure" type="partial_hill">
                        <half_max>0.5</half_max>
                        <hill_power>2</hill_power>
                        <applies_to_dead>0</applies_to_dead>
                    </signal>
                </decreasing_signals>
                <increasing_signals>
                    <max_response>1.0</max_response>
                    <signal name="contact with dead" type="heaviside">
                        <threshold>1</threshold>
                        <applies_to_dead>0</applies_to_dead>
                    </signal>
                </increasing_signals>
            </signal>
        </signal>
    </increasing_signals>
    <decreasing_signals />
</behavior>
```

# Automatic conversion of formats
The Julia package [PhysiCellXMLRules.jl](https://github.com/drbergman/PhysiCellXMLRules.jl/) can convert a CSV rules file to the new XML format and convert it back. The conversion to CSV is currently lossy and not (yet?) readable by PhysiCell.

See the following for examples of the conversion from this unreleased branch.

## From CSV to XML
For the following content of `rules.csv`:
```csv
cell_type_1,pressure,decreases,cycle entry,0.0,0.5,4.0,0
```
Running
```julia
using PhysiCellXMLRules
path_to_csv = "rules.csv"
path_to_xml = "rules.xml"
writeXMLRules(path_to_xml, path_to_csv)
```
will produce the following XML:
```xml
<?xml version="1.0" encoding="utf-8"?>
<behavior_rulesets>
  <behavior_ruleset name="cell_type_1">
    <behavior name="cycle entry">
      <type>setter</type>
      <mediator>decreasing_dominant</mediator>
      <decreasing_signals>
        <aggregator>multivariate_hill</aggregator>
        <max_response>0.0</max_response>
        <signal name="pressure" type="partial_hill">
          <applies_to_dead>0</applies_to_dead>
          <half_max>0.5</half_max>
          <hill_power>4.0</hill_power>
        </signal>
      </decreasing_signals>
    </behavior>
  </behavior_ruleset>
</behavior_rulesets>
```

## From XML to CSV
With the above XML file, running
```julia
using PhysiCellXMLRules
path_to_xml = "rules.xml"
path_to_csv = "rules_export.csv"
exportCSVRules(path_to_csv, path_to_xml)
```
will produce the following CSV:
```csv
// XML Rules Export

// This file was generated by PhysiCellXMLRules.jl

// README:
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// cell_type,(s₁) signal (s₂),response (s₃),behavior,max_response,p₁,p₂,applies_to_dead
// s₁: if 'decreasing', then as the signal decreases, the behavior <response>; (s₁) omitted if 'increasing'
// s₂: if 'from <value>', then the signal affects the behavior from <value> onwards (depending on s₁); omitted implies 'from 0.0'
// s₃: the type of function (one of 'partial_hill', 'hill', 'linear', 'heaviside'); (s₃) omitted if 'partial_hill'
// pᵢ: the parameters for the function:
//   - half_max and hill_power for partial_hill and hill
//   - signal_min and signal_max for linear
//   - threshold for heaviside; leave p₂ empty, i.e. '...,<threshold>,,...'
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// cell_type_1
// └─cycle entry
//   └─setter using a decreasing_dominant mediator
//     └─decreasing to 0.0 using a multivariate_hill aggregator
cell_type_1,pressure,decreases,cycle entry,0.0,0.5,4.0,0
```

# To-dos
- [X] Export parsed rules to CSV
- [ ] Export to human-readable formats
- [ ] Read CSVs into PhysiCell with extensions
- [ ] Integrate into studio
- [ ] Finish adding docstrings

# Reviewer notes
The core changes are in `core/PhysiCell_rules_extended.cpp` and `core/PhysiCell_rules_extrended.h`. After just testing the sample projects, I recommend looking through those.