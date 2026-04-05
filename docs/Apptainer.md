# Apptainer

## Overview

The image is fully self-contained: the latest BuscoClade release is bundled
at `/opt/buscoclade/` during the build. Input data is mounted from the host
at runtime via `-B` (bind mount).

Two image variants are available:

| Variant | Tag | Contents | First run |
|---|---|---|---|
| **Minimal** | `:latest`, `:vX.Y.Z-minimal` | Snakemake + pipeline code | conda envs created on first run (~10–15 min) |
| **Full** | `:vX.Y.Z-full` | + all pre-built conda envs | starts immediately |

---

## 1. Pull pre-built image from GHCR

The easiest way to get started — no build required:

```bash
# Latest minimal image
apptainer pull buscoclade.sif oras://ghcr.io/tomarovsky/buscoclade:latest

# Specific release, minimal
apptainer pull buscoclade.sif oras://ghcr.io/tomarovsky/buscoclade:vX.Y.Z-minimal

# Specific release, full (conda envs pre-built, no wait on first run)
apptainer pull buscoclade.sif oras://ghcr.io/tomarovsky/buscoclade:vX.Y.Z-full
```

---

## 2. Build locally

If you prefer to build the image yourself:

### Minimal image (default)

```bash
apptainer build --fakeroot apptainer/buscoclade.sif apptainer/buscoclade.def
```

### Full image (all conda environments pre-built)

```bash
apptainer build --fakeroot --build-arg BUILD_ENVS=true \
    apptainer/buscoclade.sif apptainer/buscoclade.def
```

> If `--fakeroot` is not available on your system, build on a machine with
> root access and transfer the `.sif` file — it is fully portable.

---

## 3. Run the pipeline

The `%runscript` handles `--use-conda`, `--conda-frontend`, `--conda-prefix`,
and `--directory` automatically — pass any Snakemake arguments directly after
the image name.

> **Note:** Use `--unsquash` to avoid squashfuse timeout issues on network
> filesystems and clusters.

### Basic run

```bash
mkdir -p /host/path/to/.snakemake
mkdir -p /host/path/to/results

apptainer run --unsquash \
    -B /host/path/to/.snakemake:/opt/buscoclade/.snakemake \
    -B /host/path/to/input/:/opt/buscoclade/input \
    -B /host/path/to/results/:/opt/buscoclade/results \
    -B /host/path/to/busco_datasets/ortho_odb12/:/opt/ortho_odb12/ \
    buscoclade.sif \
    --config busco_dataset_path="/opt/ortho_odb12/" \
    --cores 10
```

> **Note:** When using the full image, add `use_existing_envs=True` to the `--config` flag.

### Dry run

```bash
apptainer run --unsquash \
    ... \
    buscoclade.sif --dry-run
```

### Using a custom config file from the host

```bash
apptainer run --unsquash \
    ... \
    -B /host/path/to/my_config.yaml:/opt/buscoclade/config/my_config.yaml \
    buscoclade.sif \
    --configfile /opt/buscoclade/config/my_config.yaml
```

### Run on a SLURM cluster

```bash
apptainer run --unsquash \
    ... \
    buscoclade.sif \
    --profile /opt/buscoclade/profile/slurm
```

---

## 4. Pipeline options

All pipeline parameters are passed as `--config key=value` or via a config
file. For a full description of every parameter see [[Configuration]].

The most commonly overridden options at runtime:

```bash
apptainer run --unsquash ... buscoclade.sif \
    --config \
        busco_dataset_path="/opt/ortho_odb12/" \
        alignment="prank" \
        filtration="clipkit" \
        clipkit_params="--mode smart-gap --codon" \
        iqtree_params="-keep-ident -m TESTNEW -bb 1000 -o OUTGROUP_SPECIES" \
        astral_params="--support 2 --root OUTGROUP_SPECIES" \
        tree_visualization_params="--outgroup OUTGROUP_SPECIES" \
    --cores 10
```

The full default config is available inside the container at
`/opt/buscoclade/config/default.yaml` and can be inspected with:

```bash
apptainer exec buscoclade.sif cat /opt/buscoclade/config/default.yaml
```

---

## 5. Conda environment cache

### Minimal image

On first run, Snakemake creates conda environments from `workflow/envs/*.yaml`
and caches them on the host at:

```
$HOME/.cache/buscoclade/conda/
```

Subsequent runs reuse the cache and start immediately.

To override the cache location (useful on clusters with limited `$HOME` quota):

```bash
export BUSCOCLADE_CONDA_PREFIX=/path/to/conda_cache
apptainer run ... buscoclade.sif ...
```

### Full image

Conda environments are pre-built inside the image at `/opt/buscoclade/envs/`.
The `%runscript` detects them automatically — no extra configuration needed.

---

## 6. Directory layout inside the container

```
/opt/buscoclade/          # pipeline code (read-only, bundled at build time)
    workflow/
    config/
        default.yaml      # full list of pipeline options
    input/                # ← mount: -B /host/input:/opt/buscoclade/input
        genomes/
        vcf_reconstruct/
        vcf2phylip/
    results/              # ← mount: -B /host/results:/opt/buscoclade/results
    profile/
        slurm/
    .release_tag          # bundled pipeline version
```

---

## 7. Troubleshooting

**Permission error writing conda environments:**
The container filesystem is read-only. Use `apptainer run buscoclade.sif ...`
(not `apptainer exec ... snakemake ...`) — the `%runscript` sets
`--conda-prefix` to a writable host path automatically.

**No space for conda environments in `$HOME`:**
```bash
export BUSCOCLADE_CONDA_PREFIX=/path/with/enough/space
```

**Cannot build with `--fakeroot`:**
Build the image on a machine with root access and transfer the `.sif` file
to the cluster — `.sif` files are fully portable.

**Check which pipeline version is bundled in the image:**
```bash
apptainer exec buscoclade.sif cat /opt/buscoclade/.release_tag
```
