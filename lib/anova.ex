defmodule Anova do
  @moduledoc """
  ANOVA (Analysis of variance)
  """

  @doc """
  Computes one-way ANOVA for the given groups.
  groups: List of samples by group.
  alpha: Significance level.
  """
  @spec one_way(list(), float()) :: map()
  def one_way(groups, alpha) do
    if length(groups) < 2 do
      raise ArgumentError, "At least 2 groups are required for ANOVA"
    end

    if Enum.any?(groups, fn group -> length(group) == 0 end) do
      raise ArgumentError, "All groups must contain at least one observation"
    end

    {ssr, sse} = ss(groups)
    df_r = length(groups) - 1
    df_e = length(List.flatten(groups)) - length(groups)

    ms_r = ssr / df_r
    ms_e = sse / df_e

    f_value = ms_r / ms_e

    f_crit = Statistics.Distributions.F.ppf(df_r, df_e).(1 - alpha)
    p_value = 1 - Statistics.Distributions.F.cdf(df_r, df_e).(f_value)

    # Effect size (eta-squared)
    effect_size = if ssr + sse > 0, do: ssr / (ssr + sse), else: 0.0

    %{
      between: %{ss: ssr, df: df_r, ms: ms_r},
      within: %{ss: sse, df: df_e, ms: ms_e},
      total: %{ss: ssr + sse, df: df_r + df_e},
      f_value: f_value,
      f_crit: f_crit,
      p: p_value,
      effect_size: effect_size,
      significant: p_value < alpha
    }
  end

  defp ss(groups) do
    group_means =
      Enum.map(groups, fn group ->
        Enum.sum(group) / length(group)
      end)

    all_observations = List.flatten(groups)
    overall_mean = Enum.sum(all_observations) / length(all_observations)

    # Calculate Sum of Squares Regression (Between groups) - SSR
    ssr =
      groups
      |> Enum.zip(group_means)
      |> Enum.map(fn {group, group_mean} ->
        length(group) * :math.pow(group_mean - overall_mean, 2)
      end)
      |> Enum.sum()

    # Calculate Sum of Squares Error (Within groups) - SSE
    sse =
      groups
      |> Enum.zip(group_means)
      |> Enum.map(fn {group, group_mean} ->
        group
        |> Enum.map(fn value -> :math.pow(value - group_mean, 2) end)
        |> Enum.sum()
      end)
      |> Enum.sum()

    {ssr, sse}
  end
end

defmodule AnovaAlternative do
  @moduledoc """
  Alternative ANOVA implementation with clearer mathematical relationships.
  """

  def one_way_verbose(groups, alpha) do
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
    group_stats =
      Enum.map(groups, fn group ->
        %{
          data: group,
          n: length(group),
          mean: Enum.sum(group) / length(group),
          sum: Enum.sum(group)
        }
      end)

    # Overall mean
    overall_mean = Enum.sum(all_observations) / n_total

    # Step 2: Calculate Sum of Squares
    # SST = SSR + SSE (fundamental ANOVA identity)

    # Total Sum of Squares
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

    # Verify the fundamental identity: SST = SSR + SSE
    if abs(sst - (ssr + sse)) > 0.001 do
      IO.puts("WARNING: SST != SSR + SSE. Check calculations!")
      IO.puts("SST: #{sst}, SSR + SSE: #{ssr + sse}")
    end

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
end
