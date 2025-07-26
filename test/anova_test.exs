defmodule ANOVATest do
  @moduledoc false

  use ExUnit.Case
  doctest ANOVA

  test "group differences are significant" do
    groups = [
      [1, 2, 1, 0],
      [2, 3, 2, 1],
      [10, 7, 10, 8],
      [5, 4, 5, 6]
    ]

    alpha = 0.05

    result = ANOVA.one_way(groups, alpha)
    IO.inspect(result)

    assert true == result.test_results.significant?
    assert result.test_results.p_value < alpha
  end
end
