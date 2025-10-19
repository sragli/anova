defmodule ANOVATest do
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

    result = ANOVA.one_way(groups)

    assert result.test_results.p_value < alpha
  end

  test "raises when fewer than 2 groups" do
    assert_raise ArgumentError, fn -> ANOVA.one_way([[1, 2]]) end
  end

  test "raises when a group has fewer than 2 numbers or non-numeric entries" do
    assert_raise ArgumentError, fn -> ANOVA.one_way([[1, 2], [3]]) end
    assert_raise ArgumentError, fn -> ANOVA.one_way([[1, 2], [3, "a"]]) end
  end
end
