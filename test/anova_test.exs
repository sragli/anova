defmodule AnovaTest do
  @moduledoc false

  use ExUnit.Case
  doctest Anova

  test "group differences are significant" do
    groups = [
      [1, 2, 1, 0],
      [2, 3, 2, 1],
      [10, 7, 10, 8],
      [5, 4, 5, 6]
    ]

    alpha = 0.05

    result = Anova.one_way(groups, alpha)
    IO.inspect(result)

    assert result[:p] < alpha == result[:significant]
  end
end
