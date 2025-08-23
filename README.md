# ANOVA

ANOVA implementation in Elixir.

ANOVA, or Analysis of Variance is a widely used statistical technique for comparing means
among multiple groups, providing insights into whether observed differences are due to
actual effects or simply random variability. It's an extension of the t-test, allowing for
comparisons among multiple groups simultaneously. One-Way ANOVA compares means of three or
more independent groups for one factor.

## Installation

The package can be installed by adding `anova` to your list of dependencies in `mix.exs`:

```elixir
def deps do
  [
    {:anova, "~> 0.5.1"}
  ]
end
```

## Usage

```elixir
groups = [[1,2,3,4,5], [6,7,8,9,10], [11,12,13,14,15]]

alpha = 0.05

anova_results = ANOVA.one_way(groups)

tukey_results = TukeyHSD.test(anova_results, alpha)
```

## Limitations

* Using rough approximations to estimate Tukey's p-value.

## Theoretical Background

One-way ANOVA (Analysis of Variance) is used to test the null hypothesis that the means of
several groups are equal against the alternative hypothesis that at least one group mean
is different.

One-way ANOVA partitions the total variability in the data into two sources: variability
between groups (due to the factor) and variability within groups (due to random error). It
then compares these two sources of variability using the F-statistic. A large F-statistic
suggests that the between-group variability is much larger than the within-group
variability, indicating a significant difference between at least two group means.
Post-hoc tests: If the one-way ANOVA result is significant (meaning at least two groups
have different means), post-hoc tests are often used to determine which specific groups
differ significantly from each other.

### Common Post-hoc Tests

* Tukey's HSD (Honestly Significant Difference):
  * Most commonly used post-hoc test
  * Controls family-wise error rate
  * Good balance of power and Type I error control

* Bonferroni Correction
  * When you have few comparisons planned in advance
  * Divides alpha by number of comparisons
  * Simple but often too conservative

* Holm-Bonferroni (Step-down)
  * More powerful than standard Bonferroni, good general choice
  * Sequential testing procedure
  * Still controls family-wise error rate

* Scheffé Test
  * Most conservative test
  * Allows for any contrast (not just pairwise)
  * Best when you have many planned contrasts

* Fisher's LSD (Least Significant Difference)
  * Least conservative (should only be used when ANOVA F is significant)
  * No correction for multiple comparisons
  * Highest power but higher Type I error risk

### Key Concepts

* Independent Variable (Factor): The variable that defines the groups being compared.
* Dependent Variable: The variable that is measured and compared between groups.
* Groups: The different levels or categories within the independent variable.
* Null Hypothesis: The assumption that all group means are equal (no difference).
* Alternative Hypothesis: The statement that at least two group means are different.
* F-statistic: A calculated value that summarizes the variability between groups compared
  to the variability within groups.
* p-value: The probability of observing the results (or more extreme results) if the null
  hypothesis is true.
