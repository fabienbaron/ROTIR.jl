using Pkg
Pkg.activate(@__DIR__)
CI = get(ENV, "CI", nothing) == "true"
using Documenter, ROTIR

makedocs(;
    modules=[ROTIR],
    sitename = "ROTIR",
    repo = "https://github.com/fabienbaron/ROTIR.jl/blob/{commit}{path}#{line}",
    # :exports, but REPORTING rather than blocking. `:none` hid the backlog entirely, which is
    # why adding an `@ref` to an undocumented export broke the build three times in one session:
    # the first thing that referenced it was also the first thing that checked it.
    #
    # MEASURED 2026-10-04: 394 exported bindings, 65 of them in a `@docs` block. 140 have NO
    # docstring at all (`provider_intensity`, `readoifits`, `create_star`, `gui`, the whole
    # `plot2d*`/`plot3d*` family) and 189 have one the manual never includes. A hard
    # `checkdocs = :exports` therefore fails with ~329 errors today, so `warnonly` keeps the
    # build green while making the list visible on every run. Remove `:missing_docs` from
    # `warnonly` once that list is short enough to hold the line.
    checkdocs = :exports,
    warnonly = [:missing_docs],
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
            "Stellar Atmosphere models" => [
                "Overview"             => "guides/atmospheres.md",
                "Korg and MARCS"       => "guides/atmospheres_korg.md",
                "Kurucz ATLAS9"        => "guides/atmospheres_kurucz.md",
            ],
            "Parametric Fitting"      => "guides/parametric_fitting.md",
            "Image Reconstruction"    => "guides/reconstruction.md",
            "Direct Imaging Methods"  => "guides/imaging_methods.md",
            "Binary Orbits"           => "guides/orbits.md",
            "Orbit Fitting"           => "guides/orbit_fitting.md",
            "Multi-resolution"        => "guides/multires.md",
            "Sky images for OITOOLS"  => "guides/sky_export.md",
            "Plotting"                => "guides/plotting.md",
        ],
        "API Reference" => [
            "Tessellation"           => "api/tessellation.md",
            "Geometry & Surfaces"    => "api/geometry.md",
            "Chi-squared & Imaging"  => "api/chi2.md",
            "Radial Regularizers"    => "api/radial_regularizers.md",
            "Binary Geometry"        => "api/binary_geometry.md",
            "Mutual Irradiation"     => "api/reflection.md",
            "Fused Polygon Fourier Transform" => "api/fused_polyft.md",
            "Rasterization & NFFT"   => "api/rasterize_nfft.md",
            "Shape Gradients"        => "api/shape_gradient.md",
            "Parametric Model & AD"  => "api/parametric_gradient.md",
            "Atmospheres & Spectra"  => "api/atmospheres.md",
            "Bootstrap Uncertainties" => "api/bootstrap.md",
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
