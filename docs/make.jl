using Documenter
using Beyn

makedocs(
    sitename = "Beyn.jl",
    modules = [Beyn],
    remotes = nothing,     # no git remote configured yet; set this once the repo is on GitHub
    checkdocs = :none,     # the source has no docstrings yet
    pages = [
        "Home" => "index.md",
        "Algorithm" => "algorithm.md",
        "API" => "api.md",
    ],
)

deploydocs(
    repo = "github.com/$(get(ENV, "GITHUB_REPOSITORY", "OWNER/Beyn.jl")).git",
    devbranch = "main",
    push_preview = true,
)
