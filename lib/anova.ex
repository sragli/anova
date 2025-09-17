defmodule ANOVA do
  @moduledoc """
  ANOVA (Analysis of variance)
  """

  @type number_list :: [number()]
  @type groups :: [number_list()]

  @type one_way_result :: %{
          summary: %{
            groups: pos_integer(),
            group_sizes: [pos_integer()],
            total_observations: pos_integer(),
            overall_mean: float(),
            group_means: [float()]
          },
          anova_table: %{
            between: %{ss: float(), df: float(), ms: float()},
            within: %{ss: float(), df: float(), ms: float()},
            total: %{ss: float(), df: float()}
          },
          test_results: %{
            f_statistic: float(),
            p_value: float(),
            eta_squared: float(),
            omega_squared: float()
          }
        }

  @doc """
  One-way ANOVA for k independent groups (unequal n allowed).

  Call `ANOVA.one_way(groups)` where `groups` is a list of lists of numbers.

  Returns a map with the computed values.
  """
  @spec one_way(groups()) :: one_way_result()

  def one_way(groups) when length(groups) < 2 do
    raise(ArgumentError, "at least 2 groups are required for ANOVA")
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
      Enum.reduce(groups, 0.0, fn g, acc ->
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
    p_value = 1.0 - f_cdf(f, df_between, df_within)

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

    means
    |> Enum.zip(ns)
    |> Enum.map(fn {m, n} -> m * n end)
    |> Enum.sum()
    |> Kernel./(total)
  end

  # ---------- F distribution ----------
  # F CDF with df1 = d1, df2 = d2
  # CDF_F(f; d1, d2) = I_{x}(a,b) with
  #   a = d1/2, b = d2/2, x = (d1*f)/(d1*f + d2)
  @spec f_cdf(float(), pos_integer(), pos_integer()) :: float()
  defp f_cdf(f, d1, d2) when f >= 0.0 and d1 > 0 and d2 > 0 do
    a = d1 / 2.0
    b = d2 / 2.0
    x = d1 * f / (d1 * f + d2)
    betainc_regularized(a, b, x)
  end

  # ---------- Regularized incomplete beta I_x(a,b) ----------
  # We implement:
  #   I_x(a,b) = B_x(a,b) / B(a,b)
  # using Lentz’s algorithm for the continued fraction in betacf,
  # and a standard symmetry transform for x > (a+1)/(a+b+2) to improve accuracy.
  @spec betainc_regularized(float(), float(), float()) :: float()
  defp betainc_regularized(a, b, x) do
    cond do
      x <= 0.0 ->
        0.0

      x >= 1.0 ->
        1.0

      true ->
        bt =
          :math.exp(
            MathUtils.lgamma(a + b) - MathUtils.lgamma(a) - MathUtils.lgamma(b) +
              a * :math.log(x) + b * :math.log(1.0 - x)
          )

        # Use symmetry to improve convergence
        if x < (a + 1.0) / (a + b + 2.0) do
          bt * betacf(a, b, x) / a
        else
          1.0 - bt * betacf(b, a, 1.0 - x) / b
        end
    end
  end

  # Continued fraction for incomplete beta (Cephes / NR style)
  @spec betacf(float(), float(), float()) :: float()
  defp betacf(a, b, x) do
    max_iter = 200
    eps = 1.0e-12
    fpmin = 1.0e-300

    aa0 = 1.0
    bb0 = 1.0 - (a + b) * x / (a + 1.0)
    bb0 = if abs(bb0) < fpmin, do: fpmin, else: bb0
    c0 = 1.0
    d0 = 1.0 / bb0
    h0 = d0

    {h, _c, _d, _aa, _bb, _m} =
      Enum.reduce_while(1..max_iter, {h0, c0, d0, aa0, bb0, 1}, fn m, {h, c, d, _aa, _bb, _} ->
        m2 = 2 * m

        # Even step
        num = m * (b - m) * x
        den = (a + m2 - 1.0) * (a + m2)
        aa = num / den

        d = 1.0 + aa * d
        d = if abs(d) < fpmin, do: fpmin, else: d
        c = 1.0 + aa / c
        c = if abs(c) < fpmin, do: fpmin, else: c
        d = 1.0 / d
        h = h * d * c

        # Odd step
        num = -(a + m) * (a + b + m) * x
        den = (a + m2) * (a + m2 + 1.0)
        aa = num / den

        d = 1.0 + aa * d
        d = if abs(d) < fpmin, do: fpmin, else: d
        c = 1.0 + aa / c
        c = if abs(c) < fpmin, do: fpmin, else: c
        d = 1.0 / d

        h_new = h * d * c

        if abs(h_new - h) < eps * abs(h) do
          {:halt, {h_new, c, d, aa, den, m}}
        else
          {:cont, {h_new, c, d, aa, den, m}}
        end
      end)

    h
  end
end
