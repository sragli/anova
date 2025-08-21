defmodule ANOVA do
  @moduledoc """
  ANOVA (Analysis of variance)
  """

  @type number_list :: [number()]
  @type groups :: [number_list()]

  @doc """
  One-way ANOVA for k independent groups (unequal n allowed).

  Call `ANOVA.one_way(groups)` where `groups` is a list of lists of numbers.

  Returns a map with the computed values.
  """
  @spec one_way(groups()) :: map()

  def one_way(groups) when length(groups) < 2 do
    raise(ArgumentError, "at least 2 groups are required for ANOVA")
  end

  def one_way([[] | _]) do
    raise(ArgumentError, "all groups must contain at least one observation")
  end

  def one_way(groups) do
    groups = Enum.map(groups, &only_numbers!/1)
    ns = Enum.map(groups, &length/1)
    k = length(groups)
    n_total = Enum.sum(ns)

    means = Enum.map(groups, &mean/1)
    overall_mean = weighted_mean(means, ns)

    # Sums of squares
    ss_between =
      Enum.zip(means, ns)
      |> Enum.reduce(0.0, fn {m, n}, acc -> acc + n * :math.pow(m - overall_mean, 2) end)

    ss_within =
      groups
      |> Enum.reduce(0.0, fn g, acc ->
        m = mean(g)
        acc + Enum.reduce(g, 0.0, fn x, a -> a + :math.pow(x - m, 2) end)
      end)

    ss_total = ss_between + ss_within

    # Degrees of freedom
    df_between = k - 1
    df_within = n_total - k
    df_total = n_total - 1

    # Mean squares
    ms_between = ss_between / df_between
    ms_within = ss_within / df_within

    # F and p-value
    f = ms_between / ms_within
    p_value = 1.0 - Statistics.Distributions.F.cdf(df_between, df_within).(f)

    # Effect sizes
    eta_squared = ss_between / ss_total
    omega_squared = (ss_between - df_between * ms_within) / (ss_total + ms_within)

    %{
      summary: %{
        groups: k,
        group_sizes: Enum.map(groups, fn g -> length(g) end),
        total_observations: n_total,
        overall_mean: overall_mean,
        group_means: means
      },
      anova_table: %{
        between: %{ss: ss_between, df: df_between, ms: ms_between},
        within: %{ss: ss_within, df: df_within, ms: ms_within},
        total: %{ss: ss_total, df: df_total}
      },
      test_results: %{
        f_statistic: f,
        p_value: p_value,
        eta_squared: eta_squared,
        omega_squared: omega_squared
      }
    }
  end

  defp only_numbers!(list) when is_list(list) do
    case Enum.all?(list, &is_number/1) and length(list) >= 2 do
      true -> list
      false -> raise(ArgumentError, "each group must be a list of at least 2 numbers")
    end
  end

  defp mean(list), do: Enum.sum(list) / length(list)

  defp weighted_mean(means, ns) do
    total = Enum.sum(ns)
    if total == 0, do: raise(ArgumentError, "sum of weights is zero")

    means
    |> Enum.zip(ns)
    |> Enum.map(fn {m, n} -> m * n end)
    |> Enum.sum()
    |> Kernel./(total)
  end
end
