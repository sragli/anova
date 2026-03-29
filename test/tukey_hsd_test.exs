defmodule TukeyHSDTest do
  use ExUnit.Case

  @groups [[1, 2, 1, 0], [2, 3, 2, 1], [10, 7, 10, 8], [5, 4, 5, 6]]

  describe "confidence interval level" do
    test "reflects alpha=0.05 as 95.0" do
      result = TukeyHSD.test(ANOVA.one_way(@groups), 0.05)
      levels = Enum.map(result.post_hoc_test.pairwise_comparisons, & &1.confidence_interval.level)
      assert Enum.all?(levels, &(&1 == 95.0))
    end

    test "reflects alpha=0.01 as 99.0" do
      result = TukeyHSD.test(ANOVA.one_way(@groups), 0.01)
      levels = Enum.map(result.post_hoc_test.pairwise_comparisons, & &1.confidence_interval.level)
      assert Enum.all?(levels, &(&1 == 99.0))
    end

    test "reflects alpha=0.10 as 90.0" do
      result = TukeyHSD.test(ANOVA.one_way(@groups), 0.10)
      levels = Enum.map(result.post_hoc_test.pairwise_comparisons, & &1.confidence_interval.level)
      assert Enum.all?(levels, &(&1 == 90.0))
    end
  end

  describe "zero within-group variance" do
    # Groups with no within-group spread produce ms_within = 0,
    # standard_error = 0, and undefined q-statistics.
    @zero_var_groups [[5, 5, 5], [6, 6, 6], [7, 7, 7]]

    test "does not crash" do
      anova = ANOVA.one_way(@zero_var_groups)
      assert anova.anova_table.within.ms == 0.0
      assert %{post_hoc_test: %{pairwise_comparisons: [_ | _]}} = TukeyHSD.test(anova, 0.05)
    end

    test "q_statistic is :infinity and p_value is 0.0" do
      result = TukeyHSD.test(ANOVA.one_way(@zero_var_groups), 0.05)

      for comp <- result.post_hoc_test.pairwise_comparisons do
        assert comp.q_statistic == :infinity
        assert comp.p_value == 0.0
      end
    end

    test "effect_size is 0.0" do
      result = TukeyHSD.test(ANOVA.one_way(@zero_var_groups), 0.05)

      for comp <- result.post_hoc_test.pairwise_comparisons do
        assert comp.effect_size == 0.0
      end
    end
  end
end
