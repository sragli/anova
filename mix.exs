defmodule ANOVA.MixProject do
  use Mix.Project

  def project do
    [
      app: :anova,
      version: "0.5.1",
      elixir: "~> 1.15",
      start_permanent: Mix.env() == :prod,
      description: description(),
      package: package(),
      deps: deps(),
      name: "ANOVA",
      source_url: "https://github.com/sragli/anova",
      docs: docs()
    ]
  end

  def application do
    [
      extra_applications: [:logger]
    ]
  end

  defp description() do
    "ANOVA implementation in Elixir."
  end

  defp package() do
    [
      files: ~w(lib .formatter.exs mix.exs README.md LICENSE CHANGELOG),
      licenses: ["Apache-2.0"],
      links: %{"GitHub" => "https://github.com/sragli/anova"}
    ]
  end

  defp docs() do
    [
      main: "ANOVA",
      extras: ["README.md", "LICENSE", "examples.livemd", "CHANGELOG"]
    ]
  end

  defp deps do
    [
      {:statistex, "~> 1.0"},
      {:ex_doc, ">= 0.0.0", only: :dev, runtime: false}
    ]
  end
end
