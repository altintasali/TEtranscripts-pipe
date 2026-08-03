#!/usr/bin/env bash
# activate-env.sh -- activate the pre-built rnaseq-star-tetranscripts conda
# environment. SOURCE this file (it sets environment variables in your current
# shell and must run in it):
#
#   source scripts/activate-env.sh [PREFIX]
#
# PREFIX defaults to $HOME/software/rnaseq-star-tetranscripts-env, matching the
# default install location used by install-env.sh. Pass a different PREFIX when
# you installed elsewhere with `./scripts/install-env.sh -o PREFIX`.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    echo "error: source this file instead of running it:" >&2
    echo "  source scripts/activate-env.sh [PREFIX]" >&2
    exit 1
fi

prefix="${1:-$HOME/software/rnaseq-star-tetranscripts-env}"
if [[ ! -f "$prefix/bin/activate" ]]; then
    echo "error: no environment found at $prefix" >&2
    echo "       install it first with: ./scripts/install-env.sh [-o $prefix]" >&2
    return 1
fi
source "$prefix/bin/activate"
