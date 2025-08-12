# ANOVA

ANOVA implementation in Elixir.

ANOVA, or Analysis of Variance is a widely used statistical technique for comparing means
among multiple groups, providing insights into whether observed differences are due to
actual effects or simply random variability. It's an extension of the t-test, allowing for
comparisons among multiple groups simultaneously. One-Way ANOVA compares means of three or
more independent groups for one factor.

## Installation

If [available in Hex](https://hex.pm/docs/publish), the package can be installed by adding
`anova` to your list of dependencies in `mix.exs`:

```elixir
def deps do
  [
    {:anova, "~> 0.1.0"}
  ]
end
```

## Key Features

* Takes a list of groups (each group is a list of observations)
* Calculates F-statistic, p-value, and all intermediate statistics
* Returns a detailed map with all ANOVA results

Complete Statistical Calculations:

* Sum of Squares Between (SSB), Within (SSW), and Total (SST)
* Mean Squares and degrees of freedom
* F-statistic and p-value approximation
* Significance testing

## Usage

```elixir
groups = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
alpha = 0.05

results = ANOVA.one_way(groups, alpha)
```