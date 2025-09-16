defmodule TukeyHSD do
  @moduledoc """
  Tukey's Honestly Significant Difference (HSD) post-hoc test
  for use with ANOVA results in the specified format.
  """

  @type int_pair :: {pos_integer(), pos_integer()}
  @type float_pair :: {float(), float()}

  @type test_result :: %{
          anova: ANOVA.one_way_result(),
          post_hoc_test: %{
            summary: %{
              test: binary(),
              total_comparisons: pos_integer(),
              significant_comparisons: non_neg_integer(),
              non_significant_comparisons: non_neg_integer(),
              significant_pairs: [int_pair()],
              difference_stats: %{
                max: float(),
                min: float(),
                median: float(),
                mean: float()
              },
              effect_size_stats: %{median: float(), mean: float()}
            },
            pairwise_comparisons: [
              %{
                standard_error: float(),
                difference: float(),
                groups: int_pair(),
                effect_size: float(),
                significant?: boolean(),
                confidence_interval: %{
                  level: float(),
                  upper: float(),
                  lower: float()
                },
                means: float_pair(),
                p_value: float(),
                q_statistic: float()
              }
            ]
          }
        }

  @doc """
  Performs Tukey's HSD test using ANOVA results.
  """
  @spec test(ANOVA.one_way_result(), float()) :: test_result()
  def test(
        %{
          summary: %{
            group_sizes: group_sizes,
            group_means: group_means
          },
          anova_table: %{
            within: %{df: df_within, ms: ms_within}
          }
        } = anova_result,
        alpha
      ) do
    comparisons =
      perform_pairwise_comparisons(group_means, group_sizes, df_within, ms_within, alpha)

    %{
      anova: anova_result,
      post_hoc_test: %{
        summary: summarize_results(comparisons),
        pairwise_comparisons: comparisons
      }
    }
  end

  defp perform_pairwise_comparisons(group_means, group_sizes, df_within, ms_within, alpha) do
    indexed_means = Enum.with_index(group_means, 1)
    n_groups = length(group_means)

    for {mean_i, i} <- indexed_means,
        {mean_j, j} <- indexed_means,
        i < j do
      n_i = Enum.at(group_sizes, i - 1)
      n_j = Enum.at(group_sizes, j - 1)

      standard_error = :math.sqrt(ms_within * (1 / n_i + 1 / n_j) / 2)

      q_critical = StudentizedRange.qtukey(1 - alpha, n_groups, df_within)

      hsd = q_critical * standard_error

      difference = abs(mean_i - mean_j)

      q_statistic = difference / standard_error

      %{
        groups: {i, j},
        means: {mean_i, mean_j},
        difference: difference,
        standard_error: standard_error,
        q_statistic: q_statistic,
        p_value: 1 - StudentizedRange.ptukey(q_statistic, n_groups, df_within),
        significant?: difference > hsd,
        confidence_interval: calculate_confidence_interval(mean_i - mean_j, hsd),
        effect_size: effect_size(mean_i, mean_j, ms_within)
      }
    end
  end

  defp calculate_confidence_interval(raw_difference, hsd) do
    %{
      lower: raw_difference - hsd,
      upper: raw_difference + hsd,
      level: 95.0
    }
  end

  # Approximation of Cohen's d
  def effect_size(group1_mean, group2_mean, ms_within) when ms_within > 0.0 do
    diff = abs(group1_mean - group2_mean)
    diff / :math.sqrt(ms_within)
  end

  defp summarize_results(comparisons) do
    total_comparisons = length(comparisons)
    significant_comparisons = Enum.count(comparisons, & &1.significant?)

    differences = Enum.map(comparisons, & &1.difference)
    effect_sizes = Enum.map(comparisons, & &1.effect_size)

    %{
      test: "Tukey's HSD",
      total_comparisons: total_comparisons,
      significant_comparisons: significant_comparisons,
      non_significant_comparisons: total_comparisons - significant_comparisons,
      significant_pairs:
        comparisons
        |> Enum.filter(& &1.significant?)
        |> Enum.map(& &1.groups),
      difference_stats: %{
        mean: Statistex.average(differences),
        median: Statistex.median(differences),
        min: Enum.min(differences),
        max: Enum.max(differences)
      },
      effect_size_stats: %{
        mean: Statistex.average(effect_sizes),
        median: Statistex.median(effect_sizes)
      }
    }
  end
end
