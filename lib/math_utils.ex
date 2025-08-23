defmodule MathUtils do
  @moduledoc """
  Utility functions used in ANOVA and post-hoc test calculations.
  """

  @doc """
  Lanczos log-gamma function.

  Returns log Γ(z) for real z > 0. Accuracy is ~1e-14 relative.
  For very small z close to 0, the reflection formula ensures stability.
  """
  def lgamma(z) do
    # Lanczos approximation (g=7, n=9)
    p = [
      0.99999999999980993,
      676.5203681218851,
      -1259.1392167224028,
      771.32342877765313,
      -176.61502916214059,
      12.507343278686905,
      -0.13857109526572012,
      9.9843695780195716e-6,
      1.5056327351493116e-7
    ]

    cond do
      z < 0.5 ->
        # Reflection formula: Γ(z)Γ(1−z) = π / sin(πz)
        :math.log(:math.pi()) - :math.log(:math.sin(:math.pi() * z)) - lgamma(1.0 - z)

      true ->
        z1 = z - 1.0
        # sum over coefficients
        x =
          Enum.reduce(1..(length(p) - 1), hd(p), fn i, acc ->
            acc + Enum.at(p, i) / (z1 + i)
          end)

        t = z1 + 7.5

        0.5 * :math.log(2.0 * :math.pi()) +
          (z1 + 0.5) * :math.log(t) - t + :math.log(x)
    end
  end
end
