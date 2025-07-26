defmodule ANOVA do
  @moduledoc """
  ANOVA (Analysis of variance)
  """

  @doc """
  Computes one-way ANOVA for the given groups.
  groups: List of samples by group.
  alpha: Significance level.
  """
  def one_way(groups, alpha) do
    if length(groups) < 2 do
      raise ArgumentError, "At least 2 groups are required for ANOVA"
    end

    if Enum.any?(groups, fn group -> length(group) == 0 end) do
      raise ArgumentError, "All groups must contain at least one observation"
    end

    # Step 1: Calculate basic statistics
    all_observations = List.flatten(groups)
    n_total = length(all_observations)
    k = length(groups)

    # Group statistics
    group_stats = calculate_group_stats(groups)

    # Overall mean
    overall_mean = Enum.sum(all_observations) / n_total

    # Step 2: Calculate Sum of Squares
    {sst, ssr, sse} = calculate_sum_of_squares(groups, all_observations, overall_mean, group_stats)

    # Step 3: Degrees of freedom
    df_between = k - 1
    df_within = n_total - k
    df_total = n_total - 1

    # Step 4: Mean squares
    ms_between = ssr / df_between
    ms_within = sse / df_within

    # Step 5: F-statistic
    f_statistic = ms_between / ms_within

    # Step 6: P-value
    p_value = 1 - Statistics.Distributions.F.cdf(df_between, df_within).(f_statistic)

    %{
      summary: %{
        groups: k,
        total_observations: n_total,
        overall_mean: overall_mean,
        group_means: Enum.map(group_stats, fn %{mean: mean} -> mean end)
      },
      anova_table: %{
        between: %{ss: ssr, df: df_between, ms: ms_between},
        within: %{ss: sse, df: df_within, ms: ms_within},
        total: %{ss: sst, df: df_total}
      },
      test_results: %{
        f_statistic: f_statistic,
        p_value: p_value,
        significant?: p_value < alpha,
        effect_size: ssr / sst
      }
    }
  end

  defp calculate_group_stats(groups) do
    Enum.map(groups, fn group ->
      %{
        data: group,
        n: length(group),
        mean: Enum.sum(group) / length(group),
        sum: Enum.sum(group)
      }
    end)
  end

  defp calculate_sum_of_squares(groups, all_observations, overall_mean, group_stats) do
    sst =
      all_observations
      |> Enum.map(fn x -> :math.pow(x - overall_mean, 2) end)
      |> Enum.sum()

    # Sum of Squares Regression (Between groups)
    ssr =
      group_stats
      |> Enum.map(fn %{n: n, mean: group_mean} ->
        n * :math.pow(group_mean - overall_mean, 2)
      end)
      |> Enum.sum()

    # Sum of Squares Error (Within groups)
    sse =
      groups
      |> Enum.zip(group_stats)
      |> Enum.map(fn {group, %{mean: group_mean}} ->
        group
        |> Enum.map(fn x -> :math.pow(x - group_mean, 2) end)
        |> Enum.sum()
      end)
      |> Enum.sum()

    # Verify the fundamental ANOVA identity: SST = SSR + SSE
    if abs(sst - (ssr + sse)) > 0.001 do
      IO.puts("WARNING: SST != SSR + SSE. Check calculations!")
      IO.puts("SST: #{sst}, SSR + SSE: #{ssr + sse}")
    end

    {sst, ssr, sse}
  end
end
