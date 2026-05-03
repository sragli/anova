defmodule StudentizedRangeTest do
  use ExUnit.Case

  describe "basic guards" do
    test "ptukey returns 0 for nonpositive q" do
      assert StudentizedRange.ptukey(0.0, 3, 5) == 0.0
      assert StudentizedRange.ptukey(-1.0, 4, :infinity) == 0.0
    end

    test "qtukey endpoint behaviour" do
      assert StudentizedRange.qtukey(0.0, 3, 10) == 0.0
      assert StudentizedRange.qtukey(1.0, 3, 10) == :infinity
    end
  end

  describe "consistency of qtukey <-> ptukey (invertibility)" do
    test "finite df: qtukey then ptukey ≈ p" do
      p = 0.95
      k = 4
      df = 10

      q = StudentizedRange.qtukey(p, k, df)
      assert is_number(q) and q > 0.0

      p_est = StudentizedRange.ptukey(q, k, df)
      assert_in_delta(p_est, p, 1.0e-6)
    end

    test "infinite df: qtukey then ptukey ≈ p" do
      p = 0.975
      k = 3
      df = :infinity

      q = StudentizedRange.qtukey(p, k, df)
      assert is_number(q) and q > 0.0

      p_est = StudentizedRange.ptukey(q, k, df)
      assert_in_delta(p_est, p, 1.0e-6)
    end
  end

  describe "large df (overflow regression)" do
    # Regression: :math.pow(df, df/2) overflows for large df, causing ArithmeticError.
    # c_prefactor/2 must be computed in log-space.
    test "ptukey does not raise for large df" do
      assert is_number(StudentizedRange.ptukey(3.0, 4, 301))
    end

    test "qtukey does not raise for large df" do
      q = StudentizedRange.qtukey(0.95, 4, 301)
      assert is_number(q) and q > 0.0
    end
  end

  describe "unsupported options" do
    test "inner_nodes: 8 raises FunctionClauseError" do
      assert_raise FunctionClauseError, fn ->
        StudentizedRange.ptukey(3.0, 3, 10, inner_nodes: 8)
      end
    end

    test "outer_nodes: 8 raises FunctionClauseError" do
      assert_raise FunctionClauseError, fn ->
        StudentizedRange.ptukey(3.0, 3, 10, outer_nodes: 8)
      end
    end
  end
end
