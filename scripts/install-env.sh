#!/usr/bin/env bash
# install-env.sh -- download and install the pre-built conda environment for
# rnaseq-star-tetranscripts, published as a GitHub Release asset alongside each
# version tag. The tarball is built for Linux x86_64 (see
# .github/workflows/release-env.yml); there is no need to install conda.
#
# Usage:
#   ./scripts/install-env.sh [-o PREFIX] [-r VERSION] [-f]
#
# Options:
#   -o PREFIX   install into PREFIX (default: $HOME/software/rnaseq-star-tetranscripts-env)
#   -r VERSION  release tag to fetch, e.g. v0.1.0 (default: latest release)
#   -f          overwrite an existing PREFIX and skip the platform check
#   -h          show this help
#
# Afterwards, either activate the environment:
#   source "$PREFIX/bin/activate"
# or just prepend its bin directory to your PATH (handy on SLURM compute nodes):
#   export PATH="$PREFIX/bin:$PATH"
set -euo pipefail

repo="altintasali/rnaseq-star-tetranscripts"
stem="rnaseq-star-tetranscripts"
default_prefix="$HOME/software/${stem}-env"

usage() {
    sed -n '2,19p' "$0"
    exit "${1:-0}"
}

prefix="$default_prefix"
version=""
force=0
while getopts ":o:r:fh" opt; do
    case $opt in
        o) prefix=$OPTARG ;;
        r) version=$OPTARG ;;
        f) force=1 ;;
        h) usage ;;
        :) echo "option -$OPTARG requires an argument" >&2; usage 1 ;;
        \?) echo "invalid option: -$OPTARG" >&2; usage 1 ;;
    esac
done

if [[ ! -e "$prefix" ]]; then
    mkdir -p "$prefix"
fi
if [[ -n "$(ls -A "$prefix")" && $force -eq 0 ]]; then
    echo "error: $prefix is not empty; pick another prefix with -o or re-run with -f" >&2
    exit 1
fi

if [[ $force -eq 0 ]]; then
    os=$(uname -s)
    arch=$(uname -m)
    if [[ "$os" != "Linux" || "$arch" != "x86_64" ]]; then
        echo "error: the pre-built environment is Linux x86_64 only (got $os/$arch);" >&2
        echo "       create the env locally instead with: conda env create -f environment.yaml" >&2
        exit 1
    fi
fi

# Resolve the release to fetch.
# RNASEQ_INSTALL_API_URL overrides the GitHub API endpoint (used by tests and
# mirrors); the response must look like a GitHub "get a release" payload.
if [[ -z "$version" ]]; then
    api="${RNASEQ_INSTALL_API_URL:-https://api.github.com/repos/$repo/releases/latest}"
else
    api="${RNASEQ_INSTALL_API_URL:-https://api.github.com/repos/$repo/releases/tags/$version}"
fi
json=$(curl -fsSL "$api" 2>/dev/null || {
    echo "error: could not query $api" >&2
    exit 1
})
version=$(printf '%s' "$json" | grep -oE '"tag_name": *"[^"]*"' | head -n1 | sed -E 's/.*"([^"]+)"[[:space:]]*$/\1/')
if [[ -z "$version" ]]; then
    echo "error: could not determine a release to download" >&2
    exit 1
fi
asset_stem="${stem}-${version}-env.tar.gz"
asset_urls=$(printf '%s' "$json" \
    | grep -oE '"browser_download_url": *"[^"]*"' \
    | sed -E 's/.*"([^"]+)"[[:space:]]*$/\1/' \
    | grep -E "${asset_stem}(\.part\..*)?$|SHA256SUMS$" \
    || true)
if [[ -z "$asset_urls" ]]; then
    echo "error: release $version has no $asset_stem assets" >&2
    exit 1
fi

work_dir=$(mktemp -d)
trap 'rm -rf "$work_dir"' EXIT

echo "Downloading $stem env assets from release $version ..."
(
    cd "$work_dir"
    while IFS= read -r url; do
        [[ -n "$url" ]] || continue
        echo "  $url"
        curl -fsSL --retry 3 -O "$url"
    done <<< "$asset_urls"
)

# Verify checksums (or at least download SHA256SUMS to confirm they line up).
if command -v sha256sum >/dev/null 2>&1; then
    cksum="sha256sum"
elif command -v shasum >/dev/null 2>&1; then
    cksum="shasum -a 256"
else
    echo "error: neither sha256sum nor shasum found on PATH" >&2
    exit 1
fi
(cd "$work_dir" && $cksum -c SHA256SUMS)

# Recombine split parts, if the release shipped them.
tarball="$work_dir/${asset_stem}"
if [[ ! -f "$tarball" ]]; then
    if compgen -G "$tarball.part.*" >/dev/null; then
        cat "$tarball".part.* > "$tarball"
    else
        echo "error: no env tarball (or split parts) found in the release" >&2
        exit 1
    fi
fi

echo "Extracting into $prefix ..."
tar -xzf "$tarball" -C "$prefix"

echo "Relocating hard-coded prefixes (conda-unpack) ..."
"$prefix/bin/conda-unpack"

echo
echo "Done. Activate it with:"
echo "  source \"$prefix/bin/activate\""
echo "or add it to your PATH (e.g. once in ~/.bashrc or on SLURM compute nodes):"
echo "  export PATH=\"$prefix/bin:\$PATH\""
echo "Then run the workflow as usual:"
echo "  snakemake --configfile config/test.yaml --cores 4   # smoke test"
echo
echo "For a shared cluster install, extract to shared storage once (e.g. -o /shared/software/${stem}-env)"
echo "so all nodes see it; no conda or container runtime is needed."
