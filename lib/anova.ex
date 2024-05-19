defmodule Anova do
  @moduledoc """
  ANOVA (Analysis of variance).
  """

  @spec one_way(list(), float()) :: map()
  @doc """
  Computes one-way ANOVA for the given groups.
  """
  def one_way(groups, alpha) do
    {ssr, sse} = ss(groups)

    df_r = length(groups) - 1
    df_e = length(List.flatten(groups)) - length(groups)

    ms_r = ssr / df_r
    ms_e = sse / df_e

    f_value = ms_r / ms_e

    f_crit = Statistics.Distributions.F.ppf(df_r, df_e).(1 - alpha)
    p_value = 1 - Statistics.Distributions.F.cdf(df_r, df_e).(f_value)

    effect_size = ssr / (ssr + sse)

    %{
      "between" => %{"ss" => ssr, "df" => df_r, "ms" => ms_r},
      "within" => %{"ss" => sse, "df" => df_e, "ms" => ms_e},
      "total" => %{"ss" => ssr + sse, "df" => df_r + df_e},
      "f_value" => f_value,
      "f_crit" => f_crit,
      "p" => p_value,
      "effect_size" => effect_size
    }
  end

  defp ss(groups) do
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

    {ssr, sse}
  end
end
