# SPDX-FileCopyrightText: © 2025 TEC <contact@tecosaur.net>
# SPDX-License-Identifier: MPL-2.0

using Documenter
using BioIdentifiers

makedocs(;
    modules=[BioIdentifiers],
    pages=[
        "Index" => "index.md",
    ],
    format=Documenter.HTML(assets=["assets/favicon.ico"]),
    sitename="BioIdentifiers.jl",
    authors="tecosaur",
    warnonly=[:missing_docs],
)

deploydocs(repo="github.com/tecosaur/BioIdentifiers.jl")
