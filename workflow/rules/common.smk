# Everything the rule files share: config, samples, references, runtime
# switches, target lists, and generated conda environments.
#
# This was one 1,217-line file doing eight unrelated jobs. It is now six, in
# dependency order -- each part uses names the ones above it defined, so the
# order here is load-bearing and the parts are not independently includable.
#
# Snakemake resolves an included path relative to the including file, so these
# resolve under workflow/rules/common/.


include: "common/config.smk"
include: "common/samples.smk"
include: "common/refs.smk"
include: "common/runtime.smk"
include: "common/targets.smk"
include: "common/envs.smk"
