using Pkg
Pkg.activate(@__DIR__)
CI = get(ENV, "CI", nothing) == "true"
using Documenter, ROTIR

makedocs(;
    modules=[ROTIR],
    sitename = "ROTIR",
    checkdocs = :none,
    doctest = false,
    format = Documenter.HTML(;
          prettyurls = CI,
          collapselevel = 2,
          repolink = "https://github.com/fabienbaron/ROTIR.jl",
      ),
    authors = "Fabien Baron and contributors",
    pages = [
        "Home" => "index.md",
        "Installation" => "install.md",
        "Guides" => [
            "Overview"                => "guides/overview.md",
            "Conventions"             => "guides/conventions.md",
            "Tessellation"            => "guides/tessellation.md",
            "Surface Types"           => "guides/surfaces.md",
            "Image Reconstruction"    => "guides/reconstruction.md",
            "Multi-resolution"        => "guides/multires.md",
            "Plotting"                => "guides/plotting.md",
        ],
        "API Reference" => [
            "Tessellation"           => "api/tessellation.md",
            "Geometry & Surfaces"    => "api/geometry.md",
            "Chi-squared & Imaging"  => "api/chi2.md",
            "Fused Polyft"           => "api/fused_polyft.md",
            "Shape Gradients"        => "api/shape_gradient.md",
            "Orbits"                 => "api/orbits.md",
            "Plotting"               => "api/plotting.md",
        ],
    ],
)

if CI
    deploydocs(;
        repo   = "github.com/fabienbaron/ROTIR.jl",
        target = "build",
    )
end
