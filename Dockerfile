FROM mambaorg/micromamba:2

ARG MAMBA_DOCKERFILE_ACTIVATE=1

COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/default-config/versions.yaml /tmp/versions.yaml

# Install Python in the base environment first
RUN micromamba install -y -c conda-forge python=3.13 pyyaml

COPY --chown=$MAMBA_USER:$MAMBA_USER <<"EOF" /tmp/install_tools.py
import yaml, subprocess, sys

with open("/tmp/versions.yaml") as f:
    v = yaml.safe_load(f)["versions"]

core  = [f"star={v['star']}", f"samtools={v['samtools']}",
         f"multiqc={v['multiqc']}",
         f"ucsc-gtftogenepred={v['ucsc_gtftogenepred']}",
         f"ucsc-genepredtobed={v['ucsc_genepredtobed']}"]
r_env = [f"rseqc={v['rseqc']}", f"tetranscripts={v['tetranscripts']}"]

def run(env_name, pkgs, channels):
    cmd = ["micromamba", "create", "-y", "-n", env_name] + channels + pkgs
    r = subprocess.run(cmd)
    if r.returncode:
        sys.exit(r.returncode)

run("core", core, ["-c", "bioconda", "-c", "conda-forge"])
run("r_env", r_env, ["-c", "bioconda", "-c", "conda-forge"])
EOF

RUN /opt/conda/bin/python /tmp/install_tools.py

ENV PATH=/opt/conda/envs/core/bin:/opt/conda/envs/r_env/bin:/opt/conda/envs/base/bin:$PATH

RUN micromamba clean --all --yes
