defmodule TukeyHSD do
  @moduledoc """
  Tukey's Honestly Significant Difference (HSD) post-hoc test
  for use with ANOVA results in the specified format.
  """

  @doc """
  Performs Tukey's HSD test using ANOVA results.
  """
  def test(
        %{
          summary: %{
            groups: n_groups,
            group_sizes: group_sizes,
            group_means: group_means
          },
          anova_table: %{
            within: %{df: df_within, ms: ms_within}
          }
        } = anova_result,
        alpha
      ) do
    q_critical = get_q_critical(n_groups, df_within, alpha)

    comparisons = perform_pairwise_comparisons(group_means, group_sizes, ms_within, q_critical)

    anova_result
    |> Map.put(:posthoc_summary, summarize_results(comparisons))
    |> Map.put(:pairwise_comparisons, comparisons)
  end

  defp perform_pairwise_comparisons(group_means, group_sizes, ms_within, q_critical) do
    indexed_means = Enum.with_index(group_means, 1)

    for {mean_i, i} <- indexed_means,
        {mean_j, j} <- indexed_means,
        i < j do
      # Get sample sizes for these specific groups
      n_i = Enum.at(group_sizes, i - 1)
      n_j = Enum.at(group_sizes, j - 1)

      # SE = sqrt(MS_within * (1/n_i + 1/n_j))
      standard_error = :math.sqrt(ms_within * (1 / n_i + 1 / n_j))

      # Calculate HSD critical difference for this pair
      hsd = q_critical * standard_error

      difference = abs(mean_i - mean_j)
      significant = difference > hsd

      q_statistic = difference / standard_error

      p_value = estimate_tukey_p_value(q_statistic, length(group_means))

      ci = calculate_confidence_interval(mean_i - mean_j, hsd)

      %{
        groups: "Group #{i} vs Group #{j}",
        group_i: i,
        group_j: j,
        mean_i: mean_i,
        mean_j: mean_j,
        difference: difference,
        raw_difference: mean_i - mean_j,
        standard_error: standard_error,
        q_statistic: q_statistic,
        hsd_critical: hsd,
        p_value: p_value,
        significant?: significant,
        confidence_interval: ci,
        effect_size: calculate_cohens_d(difference, standard_error)
      }
    end
  end

  # Approximates p-value for Tukey's HSD q-statistic using normal approximation.
  # Works best when df is large.
  defp estimate_tukey_p_value(q_statistic, n_groups) do
    phi = fn x ->
      0.5 * (1.0 + :math.erf(x / :math.sqrt(2.0)))
    end

    cdf_val = phi.(q_statistic / :math.sqrt(2.0))
    1.0 - :math.pow(cdf_val, n_groups - 1)
  end

  defp calculate_confidence_interval(raw_difference, hsd) do
    %{
      lower: raw_difference - hsd,
      upper: raw_difference + hsd,
      level: 95.0
    }
  end

  defp calculate_cohens_d(difference, standard_error) do
    # Rough approximation of Cohen's d
    # This is simplified - true Cohen's d would need pooled standard deviation
    difference / (standard_error * :math.sqrt(2))
  end

  defp get_q_critical(k, df, alpha) do
    get_q_critical_lookup(k, df, alpha)
  end

  # Comprehensive lookup table for Studentized Range critical values
  defp get_q_critical_lookup(k, df, 0.05 = alpha) do
    lookup_table = %{
      {2, 5} => 3.64,
      {2, 6} => 3.46,
      {2, 7} => 3.34,
      {2, 8} => 3.26,
      {2, 9} => 3.20,
      {2, 10} => 3.15,
      {2, 12} => 3.08,
      {2, 15} => 3.01,
      {2, 20} => 2.95,
      {2, 30} => 2.89,
      {2, 60} => 2.83,
      {3, 5} => 4.60,
      {3, 6} => 4.34,
      {3, 7} => 4.16,
      {3, 8} => 4.04,
      {3, 9} => 3.95,
      {3, 10} => 3.88,
      {3, 12} => 3.77,
      {3, 15} => 3.67,
      {3, 20} => 3.58,
      {3, 30} => 3.49,
      {3, 60} => 3.40,
      {4, 5} => 5.22,
      {4, 6} => 4.90,
      {4, 7} => 4.68,
      {4, 8} => 4.53,
      {4, 9} => 4.41,
      {4, 10} => 4.33,
      {4, 12} => 4.20,
      {4, 15} => 4.08,
      {4, 20} => 3.96,
      {4, 30} => 3.85,
      {4, 60} => 3.74,
      {5, 5} => 5.67,
      {5, 6} => 5.30,
      {5, 7} => 5.06,
      {5, 8} => 4.89,
      {5, 9} => 4.76,
      {5, 10} => 4.65,
      {5, 12} => 4.51,
      {5, 15} => 4.37,
      {5, 20} => 4.23,
      {5, 30} => 4.10,
      {5, 60} => 3.98
    }

    Map.get(lookup_table, {k, df}) || interpolate_q_critical(k, df, alpha, lookup_table)
  end

  defp get_q_critical_lookup(k, df, 0.01 = alpha) do
    lookup_table = %{
      {3, 5} => 6.98,
      {3, 6} => 6.33,
      {3, 7} => 5.92,
      {3, 8} => 5.64,
      {3, 9} => 5.43,
      {3, 10} => 5.27,
      {3, 12} => 5.05,
      {3, 15} => 4.84,
      {3, 20} => 4.64,
      {3, 30} => 4.45,
      {3, 60} => 4.25,
      {4, 5} => 8.12,
      {4, 6} => 7.33,
      {4, 7} => 6.85,
      {4, 8} => 6.51,
      {4, 9} => 6.25,
      {4, 10} => 6.04,
      {4, 12} => 5.74,
      {4, 15} => 5.49,
      {4, 20} => 5.24,
      {4, 30} => 4.99,
      {4, 60} => 4.75
    }

    Map.get(lookup_table, {k, df}) || interpolate_q_critical(k, df, alpha, lookup_table)
  end

  defp get_q_critical_lookup(k, df, alpha) do
    # For other alpha values, use 0.05 as base and adjust
    base_q = get_q_critical_lookup(k, df, 0.05)

    cond do
      alpha <= 0.01 -> base_q * 1.35
      alpha <= 0.025 -> base_q * 1.20
      alpha >= 0.10 -> base_q * 0.85
      true -> base_q
    end
  end

  defp interpolate_q_critical(k, df, _alpha, lookup_table) do
    # Find nearest values for interpolation
    df_values =
      lookup_table
      |> Map.keys()
      |> Enum.filter(fn {key_k, _} -> key_k == k end)
      |> Enum.map(fn {_, key_df} -> key_df end)
      |> Enum.sort()

    case df_values do
      [] ->
        # Conservative fallback
        3.5

      [single_df] ->
        Map.get(lookup_table, {k, single_df}, 3.5)

      _ ->
        # Linear interpolation between closest df values
        {lower_df, upper_df} = find_interpolation_bounds(df, df_values)
        lower_q = Map.get(lookup_table, {k, lower_df}, 3.5)
        upper_q = Map.get(lookup_table, {k, upper_df}, 3.5)

        if lower_df == upper_df do
          lower_q
        else
          # Linear interpolation
          weight = (df - lower_df) / (upper_df - lower_df)
          lower_q + weight * (upper_q - lower_q)
        end
    end
  end

  defp find_interpolation_bounds(df, df_values) do
    lower = df_values |> Enum.filter(&(&1 <= df)) |> Enum.max(fn -> hd(df_values) end)
    upper = df_values |> Enum.filter(&(&1 >= df)) |> Enum.min(fn -> List.last(df_values) end)
    {lower, upper}
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
        mean: Statistics.mean(differences),
        median: Statistics.median(differences),
        min: Enum.min(differences),
        max: Enum.max(differences)
      },
      effect_size_stats: %{
        mean: Statistics.mean(effect_sizes),
        median: Statistics.median(effect_sizes)
      }
    }
  end
end
