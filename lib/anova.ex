defmodule Anova do
  @moduledoc """
  ANOVA (Analysis of variance).
  """

  @spec one_way(list(), float()) :: list()
  @doc """
  Computes one-way ANOVA for the given groups.
  """
  def one_way(groups, alpha) do
    group_means =
      for g <- groups do
        Enum.sum(g) / length(g)
      end

    overall_mean = Enum.sum(group_means) / length(group_means)

    ssr =
      for i <- 0..(length(groups) - 1) do
        length(Enum.at(groups, i)) * :math.pow(Enum.at(group_means, i) - overall_mean, 2)
      end
      |> Enum.sum()

    sse =
      for i <- 0..(length(groups) - 1) do
        for v <- Enum.at(groups, i) do
          :math.pow(v - Enum.at(group_means, i), 2)
        end
      end
      |> List.flatten()
      |> Enum.sum()

    df_r = length(groups) - 1
    df_e = length(List.flatten(groups)) - length(groups)

    ms_r = ssr / df_r
    ms_e = sse / df_e

    f_value = ms_r / ms_e

    f_crit = Statistics.Distributions.F.ppf(df_r, df_e).(1 - alpha)
    p_value = 1 - Statistics.Distributions.F.cdf(df_r, df_e).(f_value)

    [
      %{"id" => "Between groups", "ss" => ssr, "df" => df_r, "ms" => ms_r, "fvalue" => f_value, "p" => p_value, "fcrit" => f_crit},
      %{"id" => "Within groups", "ss" => sse, "df" => df_e, "ms" => ms_e},
      %{"id" => "Total", "ss" => ssr + sse, "df" => df_r + df_e}
    ]
  end
end
