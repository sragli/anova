defmodule StudentizedRange do
  @moduledoc """
  Studentized range distribution for Tukey's HSD.

  Provides:
    * ptukey(q, k, df) -> CDF P(Q ≤ q) for k means and residual df
    * qtukey(p, k, df) -> upper-tail quantile q such that P(Q ≤ q) = p

  The implementation follows the standard double-integral representation:
    F(q; k, ν) = C * ∫_{s=0}^{∞} s^{ν-1} φ(√ν s) ∫_{z=-∞}^{∞} φ(z) [Φ(z + q s) - Φ(z)]^{k-1} dz ds

  with C = √(2π) * k * ν^{ν/2} / (Γ(ν/2) * 2^{(ν/2 - 1)}).

  Special case ν=∞ uses:
    F(q; k) = k ∫ φ(z) [Φ(z+q) - Φ(z)]^{k-1} dz

  References:
    - Copenhaver & Holland (1988), J. Stat. Comput. Simul. (distribution of the maximum studentized range).
    - R ‘stats’ :: ptukey/qtukey documentation (algorithmic details and references).
    - Wikipedia: “Studentized range distribution” (integral forms).
  """

  @type k :: pos_integer()
  @type df :: pos_integer() | :infinity

  # ===== Public API ============================================================

  @doc """
  CDF of the studentized range: P(Q ≤ q) for `k` means and residual `df`.

  * `q`  — nonnegative real
  * `k`  — number of group means (k ≥ 2)
  * `df` — residual degrees of freedom (integer ≥ 1) or :infinity

  Options:
    * :inner_nodes (default 16)   — Gauss–Legendre nodes for inner z-integral
    * :outer_nodes (default 16)   — Gauss–Legendre nodes for outer s-integral
    * :z_limit    (default 8.0)   — truncate z ∈ [−L, L]
    * :s_limit    (default 6.0)   — truncate s ∈ [0, L] (effective for moderate/large ν)
    * :panels     (default 32)    — panel count per dimension (adaptive-like)

  Returns a float in [0,1].
  """
  @spec ptukey(float(), pos_integer(), pos_integer() | :infinity) :: float()

  def ptukey(q, k, df, opts \\ [])

  def ptukey(q, _k, _df, _opts) when q <= 0.0, do: 0.0
  def ptukey(_q, k, _df, _opts) when k < 2, do: raise(ArgumentError, "k must be ≥ 2")

  def ptukey(q, k, :infinity, opts) do
    # F(q; k) = k ∫ φ(z) [Φ(z+q) - Φ(z)]^{k-1} dz over z ∈ ℝ
    inner_nodes = Keyword.get(opts, :inner_nodes, 16)
    zlim = Keyword.get(opts, :z_limit, 8.0)
    panels = Keyword.get(opts, :panels, 32)

    integrand = fn z ->
      phi(z) * :math.pow(cdf_norm(z + q) - cdf_norm(z), k - 1)
    end

    k * integrate_gl(-zlim, zlim, integrand, inner_nodes, panels)
  end

  def ptukey(q, k, df, opts) when is_integer(df) and df >= 1 do
    # General finite-df case (double integral in s and z)
    inner_nodes = Keyword.get(opts, :inner_nodes, 16)
    outer_nodes = Keyword.get(opts, :outer_nodes, 16)
    zlim = Keyword.get(opts, :z_limit, 8.0)
    panels = Keyword.get(opts, :panels, 32)

    c = c_prefactor(k, df)

    # outer integrand in s ∈ [0, ∞) → truncate to [0, s_limit]
    # For typical Tukey settings, s_limit ≈ 6 is ample; you can raise if needed.
    s_limit = Keyword.get(opts, :s_limit, 6.0)

    outer_integrand = fn s ->
      # weight from chi distribution piece: s^(ν-1) * φ(√ν s)
      w = :math.pow(s, df - 1) * phi(:math.sqrt(df) * s)

      # inner integral over z ∈ (−∞, ∞) → truncate to [−zlim, zlim]
      inner = fn z ->
        a = cdf_norm(z + q * s) - cdf_norm(z)
        phi(z) * :math.pow(max(a, 0.0), k - 1)
      end

      w * integrate_gl(-zlim, zlim, inner, inner_nodes, panels)
    end

    c * integrate_gl(0.0, s_limit, outer_integrand, outer_nodes, panels)
  end

  @doc """
  Quantile function for the studentized range: `q` such that P(Q ≤ q) = p.

  * `p`  — probability in (0,1)
  * `k`  — number of means
  * `df` — residual df (integer ≥ 1) or :infinity

  Options (passed to ptukey/4 and to the root finder):
    * :tol      (default 1.0e-8) — absolute tolerance for p
    * :max_iter (default 100)
    * as in ptukey/4: :inner_nodes, :outer_nodes, :z_limit, :s_limit, :panels

  Returns a positive float.
  """
  @spec qtukey(float(), pos_integer(), pos_integer() | :infinity) :: :infinity | float()

  def qtukey(p, k, df, opts \\ [])

  def qtukey(p, _k, _df, _opts) when p <= 0.0, do: 0.0
  def qtukey(p, _k, _df, _opts) when p >= 1.0, do: :infinity

  def qtukey(p, k, df, opts) do
    # Bracket q: for Tukey, q typically in [0, ~20]
    max_iter = Keyword.get(opts, :max_iter, 100)
    tol = Keyword.get(opts, :tol, 1.0e-8)

    # Start with a heuristic upper bound that grows with k and tighter p
    q_hi = initial_qhi(p, k, df, opts)
    q_lo = 0.0

    f = fn q -> ptukey(q, k, df, opts) - p end

    # Ensure the bracket
    {q_lo, q_hi} =
      if f.(q_hi) < 0.0 do
        expand_hi(q_hi, f)
      else
        {q_lo, q_hi}
      end

    bisection(q_lo, q_hi, f, tol, max_iter)
  end

  # ===== Numerics: Gauss–Legendre quadrature on [a,b] with paneling ============

  # Precompute 16-node Gauss–Legendre abscissae/weights (symmetric, positive side).
  # These are standard constants for n=16 on [-1,1].
  @gl16_x [
    0.09501250983763744,
    0.2816035507792589,
    0.4580167776572274,
    0.6178762444026438,
    0.7554044083550030,
    0.8656312023878318,
    0.9445750230732326,
    0.9894009349916499
  ]
  @gl16_w [
    0.1894506104550685,
    0.1826034150449236,
    0.1691565193950025,
    0.1495959888165767,
    0.1246289712555339,
    0.0951585116824928,
    0.0622535239386479,
    0.0271524594117541
  ]

  defp integrate_gl(a, b, f, nodes, panels) when a < b and nodes in [8, 16] do
    # We implement for nodes=16 (default). For nodes=8, we just reuse a subset.
    # Split [a,b] into `panels` equal subintervals to get adaptive-like behavior.
    h = (b - a) / panels

    Enum.reduce(0..(panels - 1), 0.0, fn i, acc ->
      aa = a + i * h
      bb = aa + h
      acc + gl_panel(aa, bb, f, nodes)
    end)
  end

  defp gl_panel(a, b, f, 16) do
    # Transform [-1,1] → [a,b], x = (b+a)/2 + (b-a)/2 * t; dx = (b-a)/2 dt
    m = 0.5 * (a + b)
    r = 0.5 * (b - a)
    # symmetric sum
    sum =
      Enum.zip(@gl16_x, @gl16_w)
      |> Enum.reduce(0.0, fn {x, w}, s ->
        s1 = f.(m + r * x)
        s2 = f.(m - r * x)
        s + w * (s1 + s2)
      end)

    r * sum
  end

  # ===== Special functions =====================================================

  # Standard normal pdf and cdf
  @sqrt2 :math.sqrt(2.0)
  @sqrt2pi :math.sqrt(2.0 * :math.pi())

  defp phi(x), do: :math.exp(-0.5 * x * x) / @sqrt2pi
  defp cdf_norm(x), do: 0.5 * (1.0 + :math.erf(x / @sqrt2))

  # Prefactor C(k, ν) = √(2π) k ν^{ν/2} / (Γ(ν/2) 2^{(ν/2 - 1)})
  defp c_prefactor(k, df) do
    num = @sqrt2pi * k * :math.pow(df, df / 2)
    den = :math.pow(2.0, df / 2 - 1.0) * :math.exp(MathUtils.lgamma(df / 2))
    num / den
  end

  # ===== Root finding for qtukey ===============================================

  defp initial_qhi(p, k, :infinity, _opts) do
    # rough heuristic: grows with k and tail tightness
    # ~ maps p to a z-ish starting point
    base = inv_normal(0.5 + 0.5 * p)
    max(2.0, base + 2.0 + 0.6 * :math.log(k + 0.5))
  end

  defp initial_qhi(_p, k, df, _opts) when is_integer(df) do
    max(2.5, 2.0 + 0.7 * :math.log(k + 0.5) + 1.0 / :math.sqrt(:erlang.float(df)))
  end

  defp expand_hi(q_hi, f, factor \\ 1.6, max \\ 1.0e3) do
    cond do
      q_hi > max -> {0.0, q_hi}
      f.(q_hi) >= 0.0 -> {0.0, q_hi}
      true -> expand_hi(q_hi * factor, f, factor, max)
    end
  end

  defp bisection(a, b, f, tol, max_iter) do
    fa = f.(a)
    fb = f.(b)
    if fa > 0.0 or fb < 0.0, do: raise("qtukey: invalid bracket")

    do_bisect(a, b, fa, fb, f, tol, max_iter, 0)
  end

  defp do_bisect(a, b, fa, fb, f, tol, max_iter, it) do
    mid = 0.5 * (a + b)
    fm = f.(mid)

    cond do
      it >= max_iter or :erlang.abs(b - a) < 1.0e-10 or :erlang.abs(fm) < tol -> mid
      fm < 0.0 -> do_bisect(mid, b, fm, fb, f, tol, max_iter, it + 1)
      true -> do_bisect(a, mid, fa, fm, f, tol, max_iter, it + 1)
    end
  end

  # Simple inverse-normal for heuristic only (Acklam-ish 3-region rational approx).
  defp inv_normal(p) when p > 0.0 and p < 1.0 do
    # coefficients from a common approximation; adequate for heuristics
    a1 = -3.969683028665376e+01
    a2 = 2.209460984245205e+02
    a3 = -2.759285104469687e+02
    a4 = 1.383577518672690e+02
    a5 = -3.066479806614716e+01
    a6 = 2.506628277459239e+00

    b1 = -5.447609879822406e+01
    b2 = 1.615858368580409e+02
    b3 = -1.556989798598866e+02
    b4 = 6.680131188771972e+01
    b5 = -1.328068155288572e+01

    c1 = -7.784894002430293e-03
    c2 = -3.223964580411365e-01
    c3 = -2.400758277161838e+00
    c4 = -2.549732539343734e+00
    c5 = 4.374664141464968e+00
    c6 = 2.938163982698783e+00

    d1 = 7.784695709041462e-03
    d2 = 3.224671290700398e-01
    d3 = 2.445134137142996e+00
    d4 = 3.754408661907416e+00

    pl = 0.02425
    pu = 1.0 - pl

    x =
      cond do
        p < pl ->
          q = :math.sqrt(-2.0 * :math.log(p))

          (((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) /
            ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0)

        p > pu ->
          q = :math.sqrt(-2.0 * :math.log(1.0 - p))

          -(((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) /
            ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0)

        true ->
          q = p - 0.5
          r = q * q

          (((((a1 * r + a2) * r + a3) * r + a4) * r + a5) * r + a6) * q /
            (((((b1 * r + b2) * r + b3) * r + b4) * r + b5) * r + 1.0)
      end

    x
  end
end
