defmodule AnovaTest do
  @moduledoc false

  use ExUnit.Case
  doctest Anova

  test "greets the world" do
    groups = [
      [1, 2, 1, 0],
      [2, 3, 2, 1],
      [10, 7, 10, 8],
      [5, 4, 5, 6]
    ]

    alpha = 0.05

    result = Anova.one_way(groups, alpha)
    IO.inspect(result)

    assert Enum.at(result, 0)["p"] < alpha
  end
end
