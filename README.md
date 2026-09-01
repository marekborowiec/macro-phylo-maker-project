# Macro Phylo Maker Project

This repository contains:

- **`MacroPhyloMaker/`** — an R package for phylogenetic grafting
- **`project/`** — input trees and tables used to generate results

The goal is to allow fully reproducible construction of large phylogenies by:
1. Extracting clades from published trees
2. Preparing donor templates
3. Grafting clades and tips onto a backbone chronogram

This file first provides installation instructions followed by a general overview of code functionality. The final portion provides code to reproduce the latest versions of ant trees published at [ant-tree.org](www.ant-tree.org).

This workflow was developed by Marek Borowiec [marek.borowiec@colostate.edu](mailto:marek.borowiec@colostate.edu) and Eddie Pérochon [eddie.perochon@hotmail.com](mailto:eddie.perochon@hotmail.com).

If you use our approach, please cite:
```
Pérochon, E., Bertelsmeier, C., Borowiec, M.L. (2026). Assembling the ant tree of life through scalable phylogenetic synthesis. Journal, volume, pages. doi:XXXX.
```

Parts of this workflow also rely on Chrono-STA:
```
Barba-Montoya, J., Craig, J. M., & Kumar, S. (2025). Integrating phylogenies with chronology to assemble the tree of life. Frontiers in Bioinformatics, 5, 1571568.
```

And Taxonomic Addition for Complete Trees (TACT):
```
Chang, J., Rabosky, D. L., & Alfaro, M. E. (2020). Estimating diversification rates on incompletely-sampled phylogenies: theoretical concerns and practical solutions. Systematic Biology, 69, 602–611. doi:10.1093/sysbio/syz081
```

---

## Setup overview

There are two supported ways to run the workflow:

1. **Recommended: Docker environment**
   - Best for most users.
   - Avoids local R package compilation problems.
   - Avoids local Python dependency problems for Chrono-STA.
   - Includes Docker command-line tools needed for the TACT wrapper.
   - Gives the most reproducible environment.

2. **Alternative: piecemeal local installation**
   - Useful for development and debugging.
   - Requires managing system libraries, R packages, Python packages, Docker, and TACT yourself.
   - More flexible, but more likely to fail on a new machine.

If you only want to replicate the paper workflow, start with the Docker setup. If you are modifying package code or debugging system-specific issues, the local setup may be useful.

---

# Recommended setup: Docker

## Clone the repository

```bash
git clone https://github.com/marekborowiec/macro-phylo-maker-project.git
cd macro-phylo-maker-project
```

## Install Docker on the host system

[Install Docker](https://www.docker.com/get-started/) using the instructions for your operating system. Then test that Docker works:

```bash
docker --version
docker run --rm hello-world
```

If either command fails, fix Docker before continuing.

## Build the MacroPhyloMaker Docker image

From the repository root:

```bash
cd /path/to/macro-phylo-maker-project

docker build -t macrophylomaker:latest .
```

The first build can take several minutes because Docker installs R packages and Python packages. Later builds should reuse cached layers.

The image should include:

- R 4.4.1,
- required R packages,
- Python and packages needed for Chrono-STA,
- Docker command-line tools needed to launch the external TACT Docker image,
- `DOCKER_API_VERSION=1.43` for compatibility with older host Docker daemons,
- `CHRONOSTA_PYTHON=/opt/chronosta-venv/bin/python`.

## Run the Docker environment

From the repository root:

```bash
docker run --rm -it \
  -v "$PWD":"$PWD" \
  -w "$PWD" \
  -v /var/run/docker.sock:/var/run/docker.sock \
  macrophylomaker:latest
```

Alternatively (on Windows computers), you can use this command line. Just adapt "add_your_pathway"

```bash
docker run --rm -it -v "%cd%":/c/add_your_pathway/macro-phylo-maker-project -w /c/add_our_pathway/macro-phylo-maker-project -v /var/run/docker.sock:/var/run/docker.sock macrophylomaker:latest
```

This starts R inside the container.

The mount:

```bash
-v "$PWD":"$PWD" -w "$PWD"
```

is intentional. The project appears at the same absolute path inside the container as on the host. This matters because the TACT wrapper starts another Docker container, and Docker bind mounts are resolved by the host Docker daemon.

The mount:

```bash
-v /var/run/docker.sock:/var/run/docker.sock
```

lets the MacroPhyloMaker container call the host Docker daemon. This is required for `tact_runner = "docker"`.

### Security note

Mounting `/var/run/docker.sock` gives the container access to the host Docker daemon. Use this only with trusted containers that you built yourself. If you do not need TACT, you can omit the Docker socket mount:

```bash
docker run --rm -it \
  -v "$PWD":"$PWD" \
  -w "$PWD" \
  macrophylomaker:latest
```

## Load the package inside Docker

Inside the R session started by Docker:

```r
devtools::load_all("MacroPhyloMaker")
```

## Verify the Python environment for Chrono-STA

This is an important step. Do this before running `run_chronosta_grafting()`.
In the Docker image, the Python executable is defined by the environment variable `CHRONOSTA_PYTHON`.

```r
py <- Sys.getenv("CHRONOSTA_PYTHON")
py
```

Expected:

```text
/opt/chronosta-venv/bin/python
```

Use a temporary Python script to verify imports:

```r
tf <- tempfile(fileext = ".py")

writeLines(
  c(
    "import Bio",
    "import pandas",
    "import numpy",
    "import scipy",
    "import matplotlib",
    "print('Python deps OK')"
  ),
  tf
)

system2(
  py,
  tf,
  stdout = TRUE,
  stderr = TRUE
)
```

Expected output:

```text
Python deps OK
```

## Verify that TACT can run from inside Docker

Inside R:

```r
system2(
  "docker",
  c("run", "--rm", "jonchang/tact", "tact_add_taxa", "--help"),
  stdout = TRUE,
  stderr = TRUE
)
```

Expected output begins with:

```text
Usage: tact_add_taxa [OPTIONS]
```

If this works, the container can run the external TACT Docker image.

## Docker troubleshooting

### TACT cannot mount the work directory

Make sure you ran the container with:

```bash
-v "$PWD":"$PWD" -w "$PWD"
```

Do not mount the project only as `/work` if you plan to run TACT through Docker. The host Docker daemon needs to see the same path that the MacroPhyloMaker container passes to the nested TACT container.

### Output files are owned by root

By default, Docker may create output files owned by root. If this is inconvenient, run the container as your host user:

```bash
docker run --rm -it \
  --user "$(id -u):$(id -g)" \
  --group-add "$(stat -c '%g' /var/run/docker.sock)" \
  -e HOME=/tmp \
  -v "$PWD":"$PWD" \
  -w "$PWD" \
  -v /var/run/docker.sock:/var/run/docker.sock \
  macrophylomaker:latest
```

If Docker access fails with this version, use the simpler root-based command first to confirm that the workflow works.

---

# Alternative setup: piecemeal local installation

Use this route if you do not want to use Docker for the R environment, or if you are developing the package locally.

## Clone the repository

```bash
git clone https://github.com/marekborowiec/macro-phylo-maker-project.git
cd macro-phylo-maker-project
```

## Install system libraries on Ubuntu

On Ubuntu, install development libraries needed by R packages:

```bash
sudo apt update
sudo apt install \
  build-essential \
  gfortran \
  make \
  cmake \
  pkg-config \
  git \
  curl \
  wget \
  ca-certificates \
  libcurl4-openssl-dev \
  libssl-dev \
  libxml2-dev \
  libfontconfig1-dev \
  libcairo2-dev \
  libfreetype6-dev \
  libharfbuzz-dev \
  libfribidi-dev \
  libpng-dev \
  libjpeg-dev \
  libtiff5-dev \
  libwebp-dev \
  libuv1-dev \
  libglpk-dev \
  libgmp3-dev \
  libmpfr-dev \
  libgsl-dev \
  python3 \
  python3-pip \
  python3-venv
```

Other Linux distributions, macOS, and Windows will require equivalent libraries.

## Install R packages

From R in the repository root:

```r
install.packages(c(
  "devtools",
  "ape",
  "phytools",
  "phangorn",
  "stringr",
  "progress",
  "igraph",
  "MonoPhy",
  "here",
  "readr",
  "data.table",
  "dplyr",
  "tidyr",
  "tibble",
  "ggplot2",
  "remotes"
))
```

Then load the package:

```r
devtools::load_all("MacroPhyloMaker")
```

## Set up Python for Chrono-STA locally

Chrono-STA requires Python packages:

- `biopython`
- `pandas`
- `numpy`
- `scipy`
- `matplotlib`

You may use Anaconda, Miniconda, a virtual environment, or a system Python. The important rule is: **choose one Python executable and use the same one for setup, verification, and `run_chronosta_grafting()`**.

### Option A: Anaconda or Miniconda

Example:

```r
py <- "/home/marek/anaconda3/bin/python3"
```

Install packages using conda from the shell:

```bash
/home/marek/anaconda3/bin/conda install -c conda-forge biopython pandas numpy scipy matplotlib
```

or using pip:

```bash
/home/marek/anaconda3/bin/python3 -m pip install biopython pandas numpy scipy matplotlib
```

#### Option B: Python virtual environment

A Python virtual environment keeps the packages required by Chrono-STA separate from other Python installations. Create the environment from a terminal in the repository root:

```bash
python3 -m venv .venv-chronosta
```

On native Windows, the Python launcher may be named `python` or `py` rather than `python3`. For example:

```powershell
python -m venv .venv-chronosta
```

or:

```powershell
py -m venv .venv-chronosta
```

The location of the Python executable inside the virtual environment depends on the operating system:

- Linux, macOS, WSL, and other Unix-like systems:

  ```text
  .venv-chronosta/bin/python
  ```

- Native Windows:

  ```text
  .venv-chronosta/Scripts/python.exe
  ```

Install the required packages using the Python executable inside the virtual environment.

On Linux, macOS, or WSL:

```bash
.venv-chronosta/bin/python -m pip install --upgrade pip setuptools wheel
.venv-chronosta/bin/python -m pip install biopython pandas numpy scipy matplotlib
```

On native Windows Command Prompt or PowerShell:

```powershell
.venv-chronosta\Scripts\python.exe -m pip install --upgrade pip setuptools wheel
.venv-chronosta\Scripts\python.exe -m pip install biopython pandas numpy scipy matplotlib
```

Then select the corresponding Python executable in R. The following code checks the Docker environment variable first, followed by the standard Unix and Windows virtual-environment locations:

```r
python_candidates <- c(
  Sys.getenv("CHRONOSTA_PYTHON"),
  here::here(".venv-chronosta", "bin", "python"),
  here::here(".venv-chronosta", "Scripts", "python.exe")
)

# Remove empty values, such as an unset CHRONOSTA_PYTHON variable.
python_candidates <- python_candidates[nzchar(python_candidates)]

existing_python <- python_candidates[file.exists(python_candidates)]

if (!length(existing_python)) {
  stop(
    paste(
      "Could not find a Python executable for Chrono-STA.",
      "Expected CHRONOSTA_PYTHON,",
      ".venv-chronosta/bin/python, or",
      ".venv-chronosta/Scripts/python.exe."
    ),
    call. = FALSE
  )
}

py <- normalizePath(existing_python[1], mustWork = TRUE)
py
```

The selected path should point to the exact Python executable into which the required packages were installed.

Do not mix operating-system environments. For example:

- Native Windows R should use a native Windows Python executable such as `.venv-chronosta/Scripts/python.exe`.
- R running inside WSL should use a Linux Python executable such as `.venv-chronosta/bin/python`.
- R running inside the MacroPhyloMaker Docker container should use the executable specified by `CHRONOSTA_PYTHON`.

### Verify Python from the R console before Chrono-STA

Before running `run_chronosta_grafting()`, verify that R can start the selected Python executable and that this exact Python environment contains all packages required by Chrono-STA.

This verification is a **test**, not an installation or environment-activation step. It does not create a virtual environment, install packages, or permanently configure R. It only confirms that:

1. the object `py` contains a valid Python executable path;
2. R can launch that executable;
3. the required Python packages are installed in that environment; and
4. those packages can be imported successfully.

First confirm that the selected executable exists:

```r
py
file.exists(py)
```

`file.exists(py)` must return:

```text
TRUE
```

Next, write and run a temporary Python test script:

```r
tf <- tempfile(fileext = ".py")

writeLines(
  c(
    "import sys",
    "print('Python executable:', sys.executable)",
    "print('Python version:', sys.version)",
    "import Bio",
    "import pandas",
    "import numpy",
    "import scipy",
    "import matplotlib",
    "print('Chrono-STA Python dependencies OK')"
  ),
  tf
)

python_test <- system2(
  py,
  tf,
  stdout = TRUE,
  stderr = TRUE
)

cat(python_test, sep = "\n")

python_status <- attr(python_test, "status")

if (!is.null(python_status) && python_status != 0L) {
  stop(
    "Chrono-STA Python verification failed with status ",
    python_status,
    ". Fix the Python environment before running run_chronosta_grafting().",
    call. = FALSE
  )
}
```

Successful output should include:

```text
Chrono-STA Python dependencies OK
```

The reported Python executable should match the value of `py`. If the test fails, do not run `run_chronosta_grafting()` yet.

#### Local installation workflow

For a piecemeal local installation, the complete sequence is:

1. Create or choose a Python environment.
2. Install the required packages into that environment.
3. Set `py` to that environment's Python executable.
4. Optionally run MacroPhyloMaker's setup and check helpers.
5. Run the robust temporary-script verification.
6. Pass the same `py` object to `run_chronosta_grafting()`.

For example:

```r
# Linux, macOS, or WSL virtual environment:
py <- here::here(".venv-chronosta", "bin", "python")

# Native Windows alternative:
# py <- here::here(".venv-chronosta", "Scripts", "python.exe")

py <- normalizePath(py, mustWork = TRUE)

setup_chronosta_env(python = py)
check_chronosta_python(py)

tf <- tempfile(fileext = ".py")

writeLines(
  c(
    "import sys",
    "print('Python executable:', sys.executable)",
    "import Bio",
    "import pandas",
    "import numpy",
    "import scipy",
    "import matplotlib",
    "print('Chrono-STA Python dependencies OK')"
  ),
  tf
)

python_test <- system2(
  py,
  tf,
  stdout = TRUE,
  stderr = TRUE
)

cat(python_test, sep = "\n")

if (!is.null(attr(python_test, "status"))) {
  stop(
    "Python verification failed. Do not run Chrono-STA yet.",
    call. = FALSE
  )
}
```

`setup_chronosta_env()` and `check_chronosta_python()` are package helpers. The temporary-script test is an additional, direct confirmation that R can use the chosen Python executable. Running the test does not activate the environment or make later Chrono-STA calls work by itself. The important action is passing the same verified executable to the wrapper:

```r
res_chronosta <- run_chronosta_grafting(
  ...,
  python = py
)
```

#### Docker installation workflow

The MacroPhyloMaker Docker image already contains a Python virtual environment and the required packages. Therefore, users normally should **not** run `setup_chronosta_env()` inside the Docker container.

Inside the Docker R session, select the packaged Python executable with:

```r
py <- Sys.getenv("CHRONOSTA_PYTHON")

if (!nzchar(py)) {
  stop(
    "CHRONOSTA_PYTHON is not set in this Docker container.",
    call. = FALSE
  )
}

py <- normalizePath(py, mustWork = TRUE)
py
```

Then run the same temporary-script verification shown above. If it succeeds, pass the same object to Chrono-STA:

```r
res_chronosta <- run_chronosta_grafting(
  ...,
  python = py
)
```

The temporary-script verification is useful in both Docker and local installations. The difference is how the Python environment is created:

- In Docker, the environment is built into the image and selected through `CHRONOSTA_PYTHON`.
- In a local installation, the user creates or selects the environment, installs the dependencies, and assigns its executable path to `py`.

## Install or access TACT locally

The recommended route is still Docker:

```bash
docker pull jonchang/tact:latest
docker run --rm jonchang/tact tact_add_taxa --help
```

Then run MacroPhyloMaker TACT with:

```r
tact_runner = "docker"
docker_image = "jonchang/tact"
```

A system TACT installation can also be used:

```r
tact_runner = "system"
tact_bin = "tact_add_taxa"
```

Use:

```r
tact_runner = "none"
```

to prepare TACT input files without running TACT.

---

# Repository structure

```text
MacroPhyloMaker/        # R package and workflow functions
project/
  backbones/            # backbone trees
  published/            # donor phylogenies
  chronosta/            # Chrono-STA donor trees and related inputs
  tables/               # graft tables, taxonomy authority, realm maps
  results/              # outputs: trees, logs, PDFs
```

---

# Important path rule

Use `here::here(...)` throughout the workflow. Do not use `setwd()`, hard-coded absolute paths, or `../` relative paths in scripts that should be portable.

Example:

```r
tree_path <- here::here(
  "project", "published", "attini",
  "hanisch2022",
  "Pbruchi_MC1_SN_mcmctree_combr1-4_rename.nwk"
)
```

---

# Workflow

## Extract clades from source trees

### Example:
```
res <- read.tree(
  here::here(
    "project", "published", "formicidae", "nelsen2018",
    "Dryad_Supplementary_File_7_ML_TREE_treepl_185.tre"
  )
)

extract_clade_with_outgroup(
  res,
  genus = NULL,
  mrca_tips = c("Solenopsis_papuana", "Solenopsis_xyloni"),
  outgroup = "sister_one",
  clean = "none",
  nonmono = "prune_extras",
  resolve_polytomies = TRUE,
  force_positive_lengths = TRUE,
  seed = 42L,
  write_tree = TRUE,
  tree_path = NULL,
  write_renames = TRUE,
  renames_path = NULL,
  write_drops = TRUE,
  drops_path = NULL,

  out_dir = here::here(
    "project", "published", "formicidae", "nelsen2018"
  )
)
```

This:

* extracts an ingroup defined by MRCA
* optionally adds the closest sister outgroup
* writes Newick + log files

The function `extract_clade_with_outgroup()` is used to extract a focal clade from a larger phylogeny and optionally append a single outgroup taxon for downstream grafting.

### Basic usage

You must provide:

- a tree (phylo object), and exactly one of:
- `genus` → extract all species in a genus or
- `mrca_tips` → extract a clade defined by the MRCA of specified tips

### Defining the ingroup

Two modes are available:

**Genus mode (genus = "Genus")**

- selects all tips matching `Genus_species`
- optionally collapses duplicate species names
- if the genus is non-monophyletic: 

  - `"prune_extras"` (default): keeps only genus members within the MRCA,
  - `"error"`: aborts

**MRCA mode (mrca_tips = c("taxon1","taxon2"))**

- extracts the clade subtended by those anchors
- note that `clean = "genus_species"` will not work in this mode; grafting functions will still clean 

### Outgroup selection

- `outgroup = "sister_one"` (default):
  - finds the sister clade to the ingroup
  - selects one tip with minimum patristic distance to the ingroup
  - appends that tip to the output tree

- `outgroup = "none"`:
  - returns the ingroup only

### Tree preprocessing

By default, the function prepares trees for downstream use:

- `resolve_polytomies = TRUE` → randomly resolves polytomies
- `force_positive_lengths = TRUE` → removes zero/negative branch lengths
- `clean = "genus_species"` → standardizes tip labels and collapses duplicates

Set seed to make these steps reproducible.

### Outputs

The function returns a list:

- `tree` → extracted subtree
- `outgroup` → chosen outgroup tip (if any)
- `ingroup_tips` → tips retained
- `paths` → output file paths

If `write_tree = TRUE`, the subtree is written automatically:
```
<stem>_with_outgroup.tre
<stem>_renamed.tsv
<stem>_dropped.tsv
```

All files are written to out_dir unless explicitly overridden.

### Typical use in this project

This function is used to:

- Extract clades from published phylogenies
- Standardize labels and remove duplicates
- Add a single outgroup tip
- Produce donor trees for clade grafting

These outputs are then referenced in the clade grafting plan file.

---

## Grafting tips
The function `run_tip_grafting()` inserts tips into a backbone phylogeny according to a predefined grafting table. This can be the first grafting step, producing a backbone used for downstream clade grafting, or used to attach isolated tips missing from grafted clades.

### Example:
```
run_tip_grafting(
  backbone_path = here::here("project", "backbones", "genus.tre"),
  plan_path     = here::here("project", "tables", "grafted_genera.tsv"),
  out_prefix    = here::here("project", "results", "grafted", "genus"),
  seed_mode     = 42,
  ingroup_anchors = c("Martialis", "Camponotus")
)
```
This applies a TSV-defined grafting plan:

* attaches new tips (in this case, genera missing from Borowiec et al. 2025)
* writes resulting tree + logs

### Inputs

- Backbone tree (`backbone_path`)
  - A phylogeny (Newick)
  - Must already contain the taxa used as anchors for grafting
  - Typically produced or curated prior to grafting

- Grafting table (`plan_path`)
  - Tab-delimited file specifying how each genus should be inserted
  - Each row corresponds to one grafting operation

**Grafting table structure**

Minimum required columns: `GraftedTip`, `Function`, `Sister`

Example:
```
GraftedTip   Function                  Sister
Poneracantha graft sister to clade     Gnamptogenys,Typhlomyrmex
Alfaria      graft sister to tip       Poneracantha
```

**Supported grafting operations**

- graft sister to tip
  - Inserts a new genus as sister to an existing tip
- graft sister to clade
  - Inserts a genus as sister to the MRCA of two or more taxa
- graft within clade random
  - Inserts a genus somewhere within a specified clade
  - Placement is randomized (controlled by `seed_mode`)

### Branch attachment along a distribution
For grafting operations, the exact point where a new lineage attaches along a branch is drawn from a continuous distribution between two relative positions on that branch.

- The attachment location is expressed as a fraction of branch length:
  - `0` = immediately at the parent node
  - `1` = at the descendant tip
- The interval for possible placement is controlled by:
  - `min_frac` → lower bound
  - `max_frac` → upper bound
- The function samples a position randomly within this interval, using a distribution defined by:
  - `shape1`, `shape2` (parameters of a Beta distribution)

**Interpretation**

- Uniform placement
  - `shape1` = 1, `shape2` = 1
  - attachment equally likely anywhere between `min_frac` and `max_frac`
- Bias toward the base of the branch
  - `shape1` < `shape2`
  - attachment closer to the parent node
- Bias toward the tip
  - `shape1` > `shape2`
  - attachment closer to the descendant lineage

**Example**
```
min_frac = 0.1
max_frac = 0.9
shape1 = 1
shape2 = 1
```
→ Attachment is drawn uniformly between 10% and 90% along the branch.

### How graft placement works
For each row:

1. Identify the placement location:
  - single taxon → tip
  - multiple taxa → MRCA of those taxa
2. Insert the new genus (`GraftedTip`) relative to that location
3. Assign branch lengths according to the specified model (or defaults)
4. Repeat for all rows in the table

### Reproducibility

`seed_mode` ensures consistent placement for stochastic operations (e.g., random grafts, branch placement drawn from distribution)
Always set a seed for reproducible pipelines

### Outputs
The function writes:
```
<out_prefix>.tre → updated backbone tree
<out_prefix>_graft_log.tsv → detailed record of graft operations
optional PDF plots if enabled
```
All outputs are written to the directory specified by `out_prefix`.

### Typical workflow
Tip grafting is the first major step in the assembly pipeline:

- Start with a genus-level backbone tree
- Add missing genera using `run_tip_grafting()`
- Use the resulting tree as input for clade grafting

---

## Clade grafting
The function `run_clade_grafting()` integrates full donor phylogenies (clades) into a backbone tree. It takes a TSV plan describing donor trees and their placement, converts donors into standardized templates, and grafts them onto the backbone while preserving branch‑length structure.

Example:
```
run_clade_grafting(
  backbone_path = here::here("project", "backbones", "genus.tre"),
  plan_path = here::here("project", "tables", "clades-to-graft-clean.tsv"),
  authority = here::here("project", "tables", "antwiki-valid-species-2Aug2026.txt"),
  out_prefix = here::here("project", "results", "grafted",
                          "backbone_clade_grafted_new_bby_test"),
  seed_mode = 42,
  chronos_select = "auto",
  ultrametric_final = "none",
  plot_pdf = TRUE,
  pdf_auto = TRUE,
  plot_cex = 0.35
)
```

This:

* reads donor phylogenies
* prepares templates (chronograms if needed)
* grafts onto backbone (tip or MRCA placement)
* ensures ultrametricity of final tree

### Inputs

- Backbone tree (`backbone_path`)
  - An ultrametric genus-level tree (typically output of tip grafting)
  - Must include anchor taxa used for placement

- Clade grafting plan (`plan_path`)
  - TSV file specifying donor phylogenies and where to graft them
  - Each row = one clade graft

- Authority file (`authority`)
  - Optional species list used to standardize and filter donor tips

**Grafting plan structure**

Minimum required columns: `MRCA`, `Phylogeny_file_path`

Example:
```
MRCA                          Phylogeny_file_path
Atta,Acromyrmex               project/published/attini/tree.tre
Camponotus                    project/published/camponotus/tree.tre
```
**Placement modes**
Placement is determined automatically from the MRCA column:

- Clade graft (MRCA mode)
  - Multiple taxa (comma-separated)
  - Donor is grafted at the MRCA of those taxa

- Tip replacement (single label)
  - Single taxon
  - Donor replaces that terminal branch

**How clade grafting works**
For each row:

1. Read donor tree
2. Prepare template:
  - clean tip labels
  - infer ingroup and (optionally) outgroup
  - convert to chronogram if needed
  - compute donor stem fraction r
3. Identify placement in backbone:
  - MRCA of anchors (clade mode), or
  - specific tip (tip mode)
4. Modify backbone:
  - drop existing taxa at the target clade
  - retain one representative tip (internally)
5. Graft donor:
  - insert scaled donor crown along a branch
  - placement depth determined by r or a distribution
6. Repeat for all rows

**Stem vs crown grafting**
Controlled by the optional `Stem_mode` column in the plan:

- outgroup (default)
  - uses donor stem length (more realistic timing)

- crown
  - ignores stem, grafts using crown-only placement

**Time scaling and chronograms**

- If a donor tree is not ultrametric:
  - converted to a chronogram using ape::chronos()
  - best-fitting model selected automatically (chronos_select = "auto")

- If already ultrametric:
  - used directly

**Outputs**
The function writes:
```
<out_prefix>.tre → final grafted tree
<out_prefix>_graft_log.tsv → all graft operations
<out_prefix>.pdf → tree visualization (optional)
<out_prefix>_tips.txt → final tip list
```

---

## Time-informed species grafting with Chrono-STA

After tip grafting and clade grafting, some species may still remain unplaced because they occur in source trees that partially overlap with trees already used for clade grafting. In these cases, choosing a single donor tree for a clade can leave species from alternative, overlapping source trees behind. The function `run_chronosta_grafting()` performs a final time-informed grafting step to recover these species while minimizing disruption to the existing backbone topology.

This step uses Chrono-STA as a controlled gap-filling procedure. The reference tree is the output of the previous MacroPhyloMaker grafting steps. Donor trees are searched for species absent from the reference tree, recalibrated to the reference timescale when needed, split into smaller subtrees, optionally prefused when they overlap, and then merged back into local regions of the reference tree. This allows additional species to be incorporated without rebuilding the entire tree from scratch.

### Note on Chrono-STA compatibility patch

The workflow downloads the upstream `chronosta.py` script from the [Chrono-STA repository](https://github.com/josebarbamontoya/chrono-sta) at runtime. To maintain compatibility with recent NumPy/Pandas versions, MacroPhyloMaker applies a small runtime patch to the temporary copy of `chronosta.py` used for each analysis. The patch replaces an in-place modification of `m.values` with a writable copy before calling `np.fill_diagonal()`.

The original downloaded script is not modified. The patch is applied only to the temporary copy written into the Chrono-STA run directory, and the log records when the patch is applied. This preserves provenance while allowing the workflow to run reproducibly with current Python environments.

Users may alternatively provide their own Chrono-STA script with `chronosta_script = "path/to/chronosta.py"`. If the script already contains the compatibility patch, MacroPhyloMaker will detect this and skip patching.

Chrono-STA is developed independently and distributed by its authors. MacroPhyloMaker does not vendor a modified copy of Chrono-STA. Instead, it downloads or uses a user-provided upstream `chronosta.py` script and applies a documented runtime compatibility patch to the temporary copy used for each analysis. Users should cite Chrono-STA and comply with its license when distributing modified Chrono-STA code.

When using this step of the workflow please cite:
```
Barba-Montoya, J., Craig, J. M., & Kumar, S. (2025). Integrating phylogenies with chronology to assemble the tree of life. Frontiers in Bioinformatics, 5, 1571568.
```

### Example

In order for this step to work, Python environment has to be set up in a way R can see it. See [Verify the Python environment for Chrono-STA](#verify-the-python-environment-for-chrono-sta) and [Set up Python for Chrono-STA locally](#set-up-python-for-chrono-sta-locally) sections above for steps necessary for this function to work under Docker and local installations.

```r
res_chronosta <- run_chronosta_grafting(
  reference_tree = here::here(
    "project", "results", "grafted",
    "genus-reconstituted-2Aug2026-rescaled-for-chronosta.tre"
  ),
  donor_tree_dir = here::here(
    "project", "chronosta", "source_trees"
  ),
  out_prefix = here::here(
    "project", "results", "grafted",
    "chronosta_gapfilled"
  ),
  split_seed = 1998,
  ref_weight = 5,
  recalibrate = TRUE,
  split_gen = TRUE,
  split_sbt = TRUE,
  prefuse = TRUE,
  monoph_restore = TRUE,
  paraph_exc =  c( "Lasius",
  "Camponotus",
  "Colobopsis",
  "Eurhophalothrix",
  "Syllophopsis"),
  python = py
)
```

### Inputs

- `reference_tree`
  - A Newick tree path or a `phylo` object.
  - Usually the output of clade grafting followed by genus reconstitution.
  - This tree should be ultrametric and should contain the species already incorporated by previous grafting steps.

- `donor_tree_dir`
  - Folder containing additional donor trees in `.nwk` or `.tre` format.
  - These trees should already have harmonized tip labels in `Genus_species` format.
  - Donor trees are scanned for species that are absent from the reference tree.

- `out_prefix`
  - Prefix for all output files.
  - Output files, logs, intermediate tables, and final trees are written using this prefix.

- `split_seed`
  - Seed controlling stochastic subtree splitting.
  - Set this value for reproducible runs.

- `ref_weight`
  - Weight assigned to the reference tree during Chrono-STA merging.
  - Internally, this is implemented by duplicating the reference tree when Chrono-STA is run.
  - Larger values make the final merged tree more strongly anchored to the reference topology.

- `paraph_exc`
  - Genera that not are allowed to remain non-monophyletic in the final cleanup step, even though they are found non-monophyletic in reference or source trees.
  - These are cases where non-monophyly is present from source phylogenies, but not accepted.

### What the function does

For each donor tree, `run_chronosta_grafting()`:

1. Identifies species present in the donor tree but absent from the reference tree.
2. Identifies shared nodes between the donor tree and reference tree.
3. Optionally recalibrates donor trees to the reference timescale using `chronos()`.
4. Splits donor trees at the genus level to avoid imposing deep donor-tree topology onto the reference.
5. Further splits donor trees into smaller subtrees containing:
   - missing species,
   - representative shared species,
   - and local outgroups where possible.
6. Optionally prefuses overlapping subtrees with Chrono-STA.
7. Merges each subtree locally with the corresponding region of the reference tree.
8. Rescales the Chrono-STA output to the original reference-tree timescale.
9. Grafts the merged subtree back onto the reference.
10. Optionally restores genus-level monophyly by removing newly introduced conflicts while preserving known non-monophyletic genera.

### Chrono-STA setup

The function can automatically download the Chrono-STA Python script if it is not found. Python must be available on the system, and the following Python packages are required:

```bash
python -m pip install biopython pandas numpy scipy matplotlib
```

Because many systems have more than one Python installation, it is safest to choose one Python executable explicitly and use it consistently for setup, dependency checks, and the final grafting call. For example, on a system using Anaconda:

```r
py <- "/home/marek/anaconda3/bin/python3"

setup_chronosta_env(python = py)
check_chronosta_python(py)
```

Then pass the same Python executable to `run_chronosta_grafting()`:

```r
res <- run_chronosta_grafting(
  reference_tree = here::here(
    "project", "results", "grafted",
    "genus-reconstituted-2Aug2026-rescaled-for-chronosta.tre"
  ),
  donor_tree_dir = here::here(
    "project", "chronosta", "source_trees"
  ),
  out_prefix = here::here(
    "project", "results", "grafted",
    "chronosta_gapfilled"
  ),
  split_seed = 1998,
  ref_weight = 5,
  recalibrate = TRUE,
  split_gen = TRUE,
  split_sbt = TRUE,
  prefuse = TRUE,
  monoph_restore = TRUE,
  paraph_exc = c( "Lasius",
  "Camponotus",
  "Colobopsis",
  "Eurhophalothrix",
  "Syllophopsis"),
  python = py
)
```

If the dependency check reports missing Python modules even after setup, confirm that R and Chrono-STA are using the same Python:

```r
system2(py, c("-c", "import sys; print(sys.executable)"), stdout = TRUE)
system2(py, c("-m", "pip", "--version"), stdout = TRUE)
check_chronosta_python(py)
```

For Anaconda users, dependencies can also be installed with conda:

```bash
/home/marek/anaconda3/bin/conda install -c conda-forge biopython pandas numpy scipy matplotlib
```

If you already have a local copy of `chronosta.py`, you can provide it explicitly:

```r
res <- run_chronosta_grafting(
  reference_tree = here::here(
    "project", "results", "grafted",
    "genus-reconstituted-2Aug2026.tre"
  ),
  donor_tree_dir = here::here(
    "project", "chronosta", "source_trees"
  ),
  out_prefix = here::here(
    "project", "results", "grafted",
    "chronosta_gapfilled"
  ),
  chronosta_script = here::here(
    "software", "chrono-sta", "code", "chronosta.py"
  )
)
```

### Outputs

The function writes several files:

```text
<out_prefix>_final_tree.nwk
<out_prefix>_final_tree_monophyletic.nwk
<out_prefix>_overlap.tsv
<out_prefix>_splits.tsv
<out_prefix>_prefusion_overlap.tsv
<out_prefix>_prefusion_groups.tsv
<out_prefix>_monophyly_removed.tsv
logs/chronosta_grafting_YYYYMMDD_HHMMSS.log
```

The returned object contains:

```r
res$tree              # final tree used for downstream analyses
res$raw_tree          # tree before optional monophyly cleanup
res$overlap           # donor tree overlap summary
res$split_log         # subtree splitting summary
res$monophyly         # monophyly cleanup details
res$paths             # file paths to outputs
res$parameters        # parameters used in the run
```

### Typical workflow

Chrono-STA grafting is used after the main tip and clade grafting steps:

1. Start with a genus-level backbone.
2. Add missing genera with `run_tip_grafting()`.
3. Add species-level clades with `run_clade_grafting()`.
4. Reconstitute genus tips removed during clade grafting, if needed, using `run_tip_grafting()`.
5. Add remaining species from overlapping donor trees using `run_chronosta_grafting()`.

This final step is intended for species that are present in published source trees but could not be incorporated by direct clade grafting because their source trees partially overlap with, but were not selected over, other donor trees.


### Chrono-STA troubleshooting

#### `setup_chronosta_env()` succeeds, but `run_chronosta_grafting()` says Python packages are missing

This usually means that setup installed packages into one Python environment, while the final grafting call is using a different Python executable. Always set a Python path explicitly and use it in all three places:

```r
py <- "/home/marek/anaconda3/bin/python3"

setup_chronosta_env(python = py)
check_chronosta_python(py)

res <- run_chronosta_grafting(
  ...,
  python = py
)
```

#### `sh: 1: Syntax error: word unexpected (expecting ")")`

This indicates that the Python dependency check is being interpreted by the shell rather than passed safely to Python. Use the patched version of `check_chronosta_python()` that writes a temporary Python script and runs that file, rather than sending an inline Python expression through `python -c`.

#### Avoid unnecessary upgrades in a shared Python environment

Chrono-STA requires `biopython`, `pandas`, `numpy`, `scipy`, and `matplotlib`, but it does not require the newest available versions. If you are using an Anaconda base environment that also supports other analyses, avoid unnecessary `--upgrade` installs unless needed. Prefer:

```bash
/home/marek/anaconda3/bin/python3 -m pip install biopython pandas numpy scipy matplotlib
```

or:

```bash
/home/marek/anaconda3/bin/conda install -c conda-forge biopython pandas numpy scipy matplotlib
```

#### Checking imports manually

A direct check outside MacroPhyloMaker is:

```r
py <- "/home/marek/anaconda3/bin/python3"

system2(
  py,
  c(
    "-c",
    "import Bio, pandas, numpy, scipy, matplotlib; print('Chrono-STA imports OK')"
  ),
  stdout = TRUE,
  stderr = TRUE
)
```

If this command works but `check_chronosta_python(py)` fails, the package helper function should be updated.

---

## Taxonomic completion with TACT

After tip grafting, clade grafting, and optional Chrono-STA gap-filling, the tree may still be incomplete relative to the current species-level taxonomy. `run_tact_grafting()` performs a final taxonomic completion step using TACT. MacroPhyloMaker does not distribute TACT; instead, it prepares TACT-ready inputs, calls an external TACT installation or Docker image, restores temporary labels, removes scaffold taxa, and writes validation reports.

When using this step of the workflow please cite:
```
Chang, J., Rabosky, D. L., & Alfaro, M. E. (20120). Estimating diversification rates on incompletely-sampled phylogenies: theoretical concerns and practical solutions. Systematic Biology, 69, 602–611. doi:10.1093/sysbio/syz081
```

### What this step does

The wrapper:

1. Reads the current backbone tree.
2. Reads the AntWiki taxonomy authority file.
3. Uses `TaxonName` as the authoritative name field for AntWiki input.
4. Keeps only binomial species names, i.e. rows of the form `Genus species`.
5. Omits trinomials, subspecies, and AntWiki rows where the genus parsed from `TaxonName` disagrees with the `Genus` column.
6. Detects non-monophyletic genera and temporarily splits them into TACT-safe pseudo-genera.
7. Optionally assigns species to biogeographic realms using type-locality country and a country-to-realm table.
8. Optionally uses biogeographic realms to assign missing species preferentially to same-realm backbone anchors within each genus.
9. Runs TACT to add missing species and simulate branching times.
10. Restores temporary labels in the final tree.
11. Removes scaffold terminals such as code-like species and placeholder labels, e.g. `Eburopone_CM02`, `Recurvidris_TH01`, and `Uwari_sp`.
12. Writes validation reports comparing the final tree with the processed taxonomy.

### AntWiki taxonomy handling

For `taxonomy_format = "antwiki"`, the wrapper uses the `TaxonName` column rather than constructing names from `Genus` and `Species`. This avoids artifacts caused by AntWiki rows where the `Species` column contains an entire binomial. Only two-token `TaxonName` values are retained:

```text
Genus species
```

Three-token names are treated as trinomials or infraspecific names and omitted. Rows where the genus in `TaxonName` does not match the `Genus` column are treated as probable table errors and omitted with a warning.

### Biogeography-aware TACT completion

If `biogeo = TRUE`, the wrapper assigns species to realms using a country-to-realm table. The table should contain either:

```text
Country    UdvardyRealm
```

or:

```text
Country    Realm
```

For this project, the country table is expected at:

```text
project/tables/country_udvardy_realm.tsv
```

The biogeography-aware workflow uses type-locality country as a proxy for biogeographic affinity. Existing backbone representatives are assigned realms when possible. Missing species are then preferentially assigned to temporary same-realm anchors within their genus. If no same-realm anchor exists, the species is assigned among all available anchors for that genus. This is especially useful for globally distributed genera where missing species should not be grafted indiscriminately across distant biogeographic clades.

### Example: biogeography-aware TACT run

```r
res_tact <- run_tact_grafting(
  backbone_tree = here::here(
    "project", "results", "grafted",
    "final_tree_monophyletic_4183sp.tre"
  ),
  taxonomy = here::here(
    "project", "tables",
    "antwiki-valid-species-2Aug2026.txt"
  ),
  out_prefix = here::here(
    "project", "results", "tact",
    "Formicidae_complete_tact_biogeo"
  ),
  taxonomy_format = "antwiki",
  tact_runner = "docker",
  docker_image = "jonchang/tact",
  outgroups = NULL,
  seed = 42,
  nonmono = "split",
  nonmono_allocation = "proportional",
  genus_only = "replace_random_species",
  species_code = "keep",
  drop_code_species_after_tact = TRUE,
  enforce_taxonomy_tip_count = FALSE,
  biogeo = TRUE,
  country_realm_map = here::here(
    "project", "tables",
    "country_udvardy_realm.tsv"
  ),
  biogeo_unknown_realm = "Unknown",
  biogeo_apply_to_all_genera = TRUE,
  write_biogeo_labelled_trees = TRUE,
  keep_temp = TRUE
)
```

### Important arguments

- `backbone_tree`: current tree to complete.
- `taxonomy`: AntWiki or other taxonomy table.
- `out_prefix`: prefix for all output files.
- `tact_runner`: `"docker"`, `"system"`, or `"none"`.
- `seed`: controls wrapper-level stochastic choices, including allocation of missing species to temporary anchors.
- `nonmono`: use `"split"` for normal treatment of non-monophyletic genera.
- `nonmono_allocation`: use `"proportional"`, `"equal"`, or `"random"`.
- `drop_code_species_after_tact`: removes code-like and placeholder scaffold taxa after TACT.
- `enforce_taxonomy_tip_count`: if `FALSE`, writes reports and warns rather than stopping on tree/taxonomy mismatches.
- `biogeo`: enables biogeography-aware allocation.
- `country_realm_map`: country-to-realm table.
- `biogeo_apply_to_all_genera`: if `TRUE`, applies realm-aware temporary anchors to all genera with backbone representatives and missing species.
- `write_biogeo_labelled_trees`: writes inspection trees with realm names appended to tip labels.
- `keep_temp`: keeps the TACT work directory for debugging.

### Main outputs

The main cleaned final tree is:

```text
<out_prefix>_tacted_cleaned.tre
```

Important audit outputs are:

```text
<out_prefix>_validation.tsv
<out_prefix>_taxonomy_mismatch_summary.tsv
<out_prefix>_tree_not_taxonomy.tsv
<out_prefix>_taxonomy_not_tree.tsv
<out_prefix>_dropped_code_species.tsv
<out_prefix>_biogeo_taxonomy_realms.tsv
<out_prefix>_nonmono_temp_name_map.tsv
<out_prefix>_nonmono_genera.tsv
<out_prefix>_skipped_taxa.tsv
<out_prefix>_excluded_taxonomy.tsv
<out_prefix>_excluded_tips.txt
```

If `write_biogeo_labelled_trees = TRUE`, the wrapper also writes:

```text
<out_prefix>_backbone_biogeo_labels.tre
<out_prefix>_tacted_cleaned_biogeo_labels.tre
```

These realm-labelled trees are intended for visual inspection and should not be treated as the primary clean analysis tree.

### Checking the TACT output

```r
ape::Ntip(res_tact$tree)
read.delim(res_tact$paths$validation)
read.delim(res_tact$paths$taxonomy_mismatch_summary)
read.delim(res_tact$paths$biogeo_taxonomy_realms)

grep("TACTTMP|TACTEXCL", res_tact$tree$tip.label, value = TRUE)
```

The last command should return `character(0)` for the cleaned tree.

### Randomized TACT replicates

The wrapper contains stochastic steps, and TACT itself performs stochastic grafting. To generate multiple plausible completed trees, run the TACT step with different seeds and different output prefixes. For stronger wrapper-level randomization, use `nonmono_allocation = "random"`.

```r
seeds <- 1:20

for (seed_i in seeds) {
  run_tact_grafting(
    backbone_tree = here::here(
      "project", "results", "grafted",
      "final_tree_monophyletic_4183sp.tre"
    ),
    taxonomy = here::here(
      "project", "tables",
      "antwiki-valid-species-2Aug2026.txt"
    ),
    out_prefix = here::here(
      "project", "results", "tact",
      sprintf("Formicidae-complete-tact-biogeo-2Aug2026_rep%03d", seed_i)
    ),
    taxonomy_format = "antwiki",
    tact_runner = "docker",
    docker_image = "jonchang/tact",
    seed = seed_i,
    nonmono = "split",
    nonmono_allocation = "random",
    genus_only = "replace_random_species",
    species_code = "keep",
    drop_code_species_after_tact = TRUE,
    enforce_taxonomy_tip_count = FALSE,
    biogeo = TRUE,
    country_realm_map = here::here(
      "project", "tables",
      "country_udvardy_realm.tsv"
    ),
    write_biogeo_labelled_trees = FALSE,
    keep_temp = FALSE
  )
}
```

For exploratory visual checks, use `write_biogeo_labelled_trees = TRUE` and `keep_temp = TRUE`. For large replicate batches, set both to `FALSE`.

### TACT troubleshooting

#### Docker cannot run TACT

Confirm Docker works:

```bash
docker run --rm jonchang/tact tact_add_taxa --help
```

#### Taxonomy and tree tip counts do not match

This does not necessarily mean the tree was not created. Inspect:

```text
<out_prefix>_taxonomy_mismatch_summary.tsv
<out_prefix>_tree_not_taxonomy.tsv
<out_prefix>_taxonomy_not_tree.tsv
```

These reports distinguish scaffold labels, omitted table-error rows, and true name mismatches.

#### Many species have `Unknown` realm

This usually means the type-locality string is not an exact match to a row in the country-to-realm map. Add rows for historical or regional strings such as `Tropical Africa`, `South America`, or `East Indies` if those should map to a realm.

---

## Replicating the ant macrophylogeny workflow from the paper

See [Setup Overview](#setup-overview) at the top of this file and [Verify the Python environment for Chrono-STA](#verify-the-python-environment-for-chrono-sta), [Set up Python for Chrono-STA locally](#set-up-python-for-chrono-sta-locally), [Verify that TACT can run from inside Docker](#verify-that-tact-can-run-from-inside-docker) for steps necessary before proceeding.

The full ant-tree workflow used in the paper can be rerun from the files provided in this repository. The goal is to make each major step explicit, reproducible, and updateable when new phylogenies or taxonomy files become available.

### Step 0. Clone repository and load package

```bash
git clone https://github.com/marekborowiec/macro-phylo-maker-project.git
cd macro-phylo-maker-project
```

Then in R:

```r
devtools::load_all("MacroPhyloMaker")
```

If the workflow will include the Chrono-STA-enabled step, select and verify Python before running Step 4.

When running inside the MacroPhyloMaker Docker container:

```r
py <- Sys.getenv("CHRONOSTA_PYTHON")
py <- normalizePath(py, mustWork = TRUE)
```

For a local virtual environment, use the operating-system-appropriate location:

```r
python_candidates <- c(
  here::here(".venv-chronosta", "bin", "python"),
  here::here(".venv-chronosta", "Scripts", "python.exe")
)

existing_python <- python_candidates[file.exists(python_candidates)]

if (!length(existing_python)) {
  stop(
    "Could not locate the local Chrono-STA Python environment.",
    call. = FALSE
  )
}

py <- normalizePath(existing_python[1], mustWork = TRUE)
```

For an existing Anaconda or Miniconda installation, set `py` directly to its Python executable, for example:

```r
py <- "/home/user/anaconda3/bin/python3"
```

or on native Windows:

```r
py <- "C:/Users/username/anaconda3/python.exe"
```

Run the Python verification described in
#verify-python-from-the-r-console-before-chrono-sta
before running Step 4. Always pass the verified executable explicitly:

```r
res_chronosta <- run_chronosta_grafting(
  ...,
  python = py
)
```

### Step 1. Graft missing genera onto the genus-level backbone

```r
backbone_genus_tree_path <- here::here("project", "backbones", "genus.tre")
backbone_genus_tree <- ape::read.tree(backbone_genus_tree_path)
backbone_noNewGenus <- ape::drop.tip(backbone_genus_tree, "NewGenus") # dropping "NewGenus" tip representing the isolated Neoponera bucki
backbone_noNewGenus_path <- here::here("project", "backbones", "genus-noNewGenus.tre")
ape::write.tree(backbone_noNewGenus, backbone_noNewGenus_path)

run_tip_grafting(
  backbone_path = here::here("project", "backbones", "genus-noNewGenus.tre"),
  plan_path     = here::here("project", "tables", "grafted_genera.tsv"),
  out_prefix    = here::here("project", "results", "grafted", "genus-2Aug2026"),
  seed_mode     = 42,
  ingroup_anchors = c("Martialis", "Camponotus")
)
```

This creates a complete genus-level backbone by inserting genera missing from the backbone chronogram according to the table in `project/tables/grafted_genera.tsv`. See that table for placements, assumptions, and sources.

#### Inspect Step 1 outputs

After this step, confirm that the genus-level grafting tree and graft log exist, then inspect the tree as an R `phylo` object. The file path and the object are different things: the path tells R where the tree is stored, while the `phylo` object is what you plot, count tips from, or pass to downstream functions.

```r
genus_tree_path <- here::here(
  "project", "results", "grafted",
  "genus-2Aug2026.tre"
)

genus_log_path <- here::here(
  "project", "results", "grafted",
  "genus_graft_log.tsv"
)

file.exists(genus_tree_path)
file.exists(genus_log_path)

genus_tree <- ape::read.tree(genus_tree_path)
ape::Ntip(genus_tree)
head(genus_tree$tip.label)

if (file.exists(genus_log_path)) {
  head(read.delim(genus_log_path))
}
```

Useful checks at this stage are:

- Did the expected missing genera appear in `genus_tree$tip.label`?
- Did the graft log record all rows from `project/tables/grafted_genera.tsv`?
- Are the intended anchor taxa still present?

### Step 2. Graft published clade-level phylogenies

```r
grafted_genera_tree_path <- here::here("project", "results", "grafted", "genus-2Aug2026.tre")
grafted_genera_tree <- ape::read.tree(grafted_genera_tree_path)
grafted_noChimaeridris <- ape::drop.tip(grafted_genera_tree, "Chimaeridris") # dropping "Chimaeridris" tip nested in Pheidole (Longino pers. comm.)
grafted_genera_noChimaeridris_path <- here::here("project", "results", "grafted", "genus-2Aug2026-noChimaeridris.tre")
ape::write.tree(grafted_noChimaeridris, grafted_genera_noChimaeridris_path)


run_clade_grafting(
  backbone_path = here::here("project", "results", "grafted", "genus-2Aug2026-noChimaeridris.tre"),
  plan_path = here::here("project", "tables", "clades-to-graft-clean.tsv"),
  authority = here::here(
    "project", "tables",
    "antwiki-valid-species-2Aug2026.txt"
  ),
  out_prefix = here::here(
    "project", "results", "grafted",
    "backbone-clade-grafted-2Aug2026"
  ),
  seed_mode = 42,
  chronos_select = "auto",
  ultrametric_final = "none",
  plot_pdf = TRUE,
  pdf_auto = TRUE,
  plot_cex = 0.35
)
```

This step reads the clade grafting table, filters donor phylogenies against the authority file, converts phylograms to chronograms when needed, prepares graftable templates, and inserts donor clades into the complete generic backbone.

In this step any taxon in donor trees missing valid specific epithet will be dropped. However, in donor trees some genera have only been sequenced for unidenitifiable morphospecies. Because of this, we manually modified certain input trees with placeholder valid names to ensure no genera are orphaned. We should now undo this:  

```r
# Replace orphan/placeholder species names with genus-only terminals in tree files.

files <- list.files(
  here::here("project", "results", "grafted"),
  pattern = "\\.(tre|tree|nwk|newick)$",
  full.names = TRUE,
  ignore.case = TRUE
)

replacements <- c(
  "Eusphinctus_furcatus" = "Eusphinctus",
  "Eburopone_easoana" = "Eburopone",
  "Lividopone_livida" = "Lividopone",
  "Lasiomyrma_gedensis" = "Lasiomyrma",
  "Dicroaspis_cryptocera" = "Dicroaspis",
  "Paratopula_andamanensis" = "Paratopula",
  "Recurvidris_browni" = "Recurvidris",
  "Aenictogiton_attenuatus" = "Aenictogiton"
)

for (f in files) {
  x <- readLines(f, warn = FALSE)
  for (old in names(replacements)) {
    x <- gsub(old, replacements[[old]], x, fixed = TRUE)
  }
  writeLines(x, f)
}
```

Confirm these genera have indeed been restored to genus-only tips:
```r
tr <- ape::read.tree(files[1])

grep(
  "Eusphinctus|Eburopone|Lividopone|Lasiomyrma|Dicroaspis|Paratopula|Recurvidris|Aenictogiton",
  tr$tip.label,
  value = TRUE
)
```

#### Inspect Step 2 outputs

The clade-grafting step writes the grafted tree, a graft log, and, if `plot_pdf = TRUE`, a PDF visualization. Inspect both the file outputs and the returned tree if you assigned the run to an object.

```r
clade_tree_path <- here::here(
  "project", "results", "grafted",
  "backbone-clade-grafted-2Aug2026.tre"
)

clade_log_path <- here::here(
  "project", "results", "grafted",
  "backbone-clade-grafted-2Aug2026_graft_log.tsv"
)

clade_pdf_path <- here::here(
  "project", "results", "grafted",
  "backbone-clade-grafted-2Aug2026.pdf"
)

file.exists(clade_tree_path)
file.exists(clade_log_path)
file.exists(clade_pdf_path)

clade_tree <- ape::read.tree(clade_tree_path)
ape::Ntip(clade_tree)
head(clade_tree$tip.label)

if (file.exists(clade_log_path)) {
  clade_log <- read.delim(clade_log_path)
  nrow(clade_log)
  head(clade_log)
}
```

Useful checks at this stage are:

- Are the major donor clades present in the tree?
- Did the graft log include one row for each intended clade-grafting operation?
- Do any expected donor species appear in the final tip labels?
- A PDF is written by default (`backbone-clade-grafted-2Aug2026.pdf`), does the placement of major clades look reasonable?

### Step 3. Reconstitute genera removed during clade grafting

We removed some genus-level terminals in the previous steps, and others are overwritten when larger clades are grafted wholesale in Step 2. These can be restored using a second tip-grafting step:

```r
run_tip_grafting(
  backbone_path = here::here(
    "project", "results", "grafted",
    "backbone-clade-grafted-2Aug2026.tre"
  ),
  plan_path = here::here(
    "project", "tables",
    "reconstitute-with-genera-post-grafting.tsv"
  ),
  out_prefix = here::here(
    "project", "results", "grafted",
    "genus-reconstituted-2Aug2026"
  ),
  seed_mode = 42,
  plot_cex = 0.35
)
```

This produces the reference tree for the final Chrono-STA-enabled species grafting step. It contains all valid ant taxa present in non-overlapping phylogenetic datasets included (3,387 in 2 August 2026 version).

#### Inspect Step 3 outputs

Check that the file exists, read it into R, and confirm that the genera you intended to restore are present.

```r
reconstituted_tree_path <- here::here(
  "project", "results", "grafted",
  "genus-reconstituted-2Aug2026.tre"
)

reconstituted_log_path <- here::here(
  "project", "results", "grafted",
  "genus-reconstituted-2Aug2026_graft_log.tsv"
)

file.exists(reconstituted_tree_path)
file.exists(reconstituted_log_path)

reconstituted_tree <- ape::read.tree(reconstituted_tree_path)
ape::Ntip(reconstituted_tree)
head(reconstituted_tree$tip.label)

if (file.exists(reconstituted_log_path)) {
  head(read.delim(reconstituted_log_path))
}
```

This tree is the immediate input to Chrono-STA-enabled gap filling. If the tree will be reused in a later R session, read it from `reconstituted_tree_path`. If you assigned the output of `run_tip_grafting()` to an object, use the object directly in the same session.

#### Rescale trees and input species names before next steps

The original backbone chronogram used here had branch length units expressed as 0.01 = 1Mya. Before we proceed, we need to rescale trees for downstream analyses such that branch length 1 = 1Mya:

```r
# This rescales complete genus-level tree
complete_genus <- read.tree(here::here("project", "results", "grafted", 
  "genus-2Aug2026.tre"))
cg_rescaled <- rescale_tree_time_units(
  complete_genus,
  factor = 100,
  out_path = here::here("project", "results", "grafted",
    "genus-2Aug2026-rescaled.tre")
)

# This rescales species tree produced by clade grafting, 
# to be used in Chrono-STA step
clade_grafted <- read.tree(here::here("project", "results", "grafted",
  "genus-reconstituted-2Aug2026.tre"))
clg_rescaled <- rescale_tree_time_units(
  clade_grafted,
  factor = 100,
  out_path = here::here("project", "results", "grafted",
    "genus-reconstituted-2Aug2026-rescaled.tre")
)
``` 

In the resulting trees some terminals are represented only by genus name. It is necessary to add species names to those before using Chrono-STA. We take them from the original Borowiec et al. 2025 chronogram:
```r
# Read clade-grafted tree
gent2ag <- read.tree(here::here(
  "project", "results", "grafted",
  "genus-reconstituted-2Aug2026-rescaled.tre"
))

# Read initial Borowiec et al. 2025 species-level tree
backbone_for_chronosta <- read.tree(here::here("project","chronosta","source_trees",
  "bakb.nwk"))

# Add species names from Borowiec et al. 2025 to monotypic genera that miss species name on the reference
list_reps <- unlist(sapply(gent2ag$tip.label[!str_detect(gent2ag$tip.label, "_")] ,
                           function(x){backbone_for_chronosta$tip.label[str_detect(backbone_for_chronosta$tip.label, x)]}))
for(i in 1:length(list_reps))
{
  gen <- names(list_reps)[i]
  gent2ag$tip.label[gent2ag$tip.label == gen] <- list_reps[i]
}

#Save output
write.tree(phy=gent2ag, here::here(
  "project", "results", "grafted",
  "genus-reconstituted-2Aug2026-rescaled-for-chronosta.tre"))
``` 

### Step 4. Add remaining species using Chrono-STA-enabled time-informed grafting

This step adds taxa with phylogenetic information that could not be included in the previous steps. That happens when multiple trees have been produced for a given clade and only one, most recent, best-resolved, etc. was chosen for the previous step.

Be sure to follow setup instructions before running this function. 

```r
res_chronosta <- run_chronosta_grafting(
  reference_tree = here::here(
    "project", "results", "grafted",
    "genus-reconstituted-2Aug2026-rescaled-for-chronosta.tre"
  ),
  donor_tree_dir = here::here(
    "project", "chronosta", "source_trees"
  ),
  out_prefix = here::here(
    "project", "results", "grafted",
    "chronosta_gapfilled"
  ),
  split_seed = 42,
  ref_weight = 5,
  recalibrate = TRUE,
  split_gen = TRUE,
  split_sbt = TRUE,
  prefuse = TRUE,
  monoph_restore = TRUE,
  paraph_exc = c(
    "Camponotus",
    "Colobopsis",
    "Syllophopsis",
    "Lasius",
    "Eurhopalothrix"
  ),
  python = py
)
```

This step searches additional donor trees for species absent from the reference tree, calibrates them to the reference timescale where possible, divides them into smaller grafting units, uses Chrono-STA to merge missing species with the relevant reference subtrees, rescales the outputs, and grafts them back into the reference tree.

#### Inspect Step 4 outputs

Chrono-STA-enabled grafting returns an R object and writes several files. The object is useful immediately in the same R session, while the files make the run reproducible and allow later inspection.

```r
# In-session object returned by run_chronosta_grafting()
class(res_chronosta$tree)
ape::Ntip(res_chronosta$tree)
head(res_chronosta$tree$tip.label)

# Important output files
res_chronosta$paths
file.exists(unlist(res_chronosta$paths))
```

If starting from a new R session, read the written tree from disk:

```r
chronosta_tree_path <- here::here(
  "project", "results", "grafted",
  "chronosta_gapfilled_final_tree_monophyletic.nwk"
)

chronosta_tree <- ape::read.tree(chronosta_tree_path)
ape::Ntip(chronosta_tree)
head(chronosta_tree$tip.label)
```

Inspect the audit tables:

```r
read.delim(here::here(
  "project", "results", "grafted",
  "chronosta_gapfilled_overlap.tsv"
)) |>
  head()

read.delim(here::here(
  "project", "results", "grafted",
  "chronosta_gapfilled_splits.tsv"
)) |>
  head()

read.delim(here::here(
  "project", "results", "grafted",
  "chronosta_gapfilled_monophyly_removed.tsv"
)) |>
  head()
```

Useful checks at this stage are:

- How many new species did Chrono-STA add?
- Which donor trees contributed missing species?
- Which species or clades were removed during optional monophyly cleanup?
- Does the resulting tree remain on the expected temporal scale?

For a PDF check:

```r
plot_tree_autosize(
  chronosta_tree,
  here::here("project", "results", "grafted", "chronosta_check.pdf"),
  cex = 0.2
)
```


### Step 5. Complete the tree with biogeography-aware TACT

After Chrono-STA-enabled grafting, use TACT to add remaining species from the AntWiki taxonomy. This step retains valid binomials from `TaxonName`, omits malformed AntWiki rows where the genus fields disagree, optionally assigns species to biogeographic realms using type-locality country, and runs TACT to complete the tree.

Be sure to follow Docker setup instructions for TACT before running this function.

```r
res_tact <- run_tact_grafting(
  backbone_tree = here::here(
    "project", "results", "grafted",
    "chronosta_gapfilled_final_tree_monophyletic.nwk"
  ),
  taxonomy = here::here(
    "project", "tables",
    "antwiki-valid-species-2Aug2026.txt"
  ),
  out_prefix = here::here(
    "project", "results", "tact",
    "Formicidae-complete-tact-biogeo-4Aug2026"
  ),
  taxonomy_format = "antwiki",
  tact_runner = "docker",
  docker_image = "jonchang/tact",
  outgroups = NULL,
  seed = 42,
  nonmono = "split",
  nonmono_allocation = "proportional",
  genus_only = "replace_random_species",
  species_code = "keep",
  drop_code_species_after_tact = TRUE,
  enforce_taxonomy_tip_count = FALSE,
  biogeo = TRUE,
  country_realm_map = here::here(
    "project", "tables",
    "country_udvardy_realm.tsv"
  ),
  biogeo_unknown_realm = "Unknown",
  biogeo_apply_to_all_genera = TRUE,
  write_biogeo_labelled_trees = TRUE,
  keep_temp = TRUE
)
```

#### Inspect the final TACT tree

If the TACT run completed in the current R session, the returned object contains both the final tree object and paths to written files:

```r
# R phylo object available in the current session
class(res_tact$tree)
ape::Ntip(res_tact$tree)
head(res_tact$tree$tip.label)

# Path to the same cleaned tree written to disk
res_tact$paths$cleaned_tree
file.exists(res_tact$paths$cleaned_tree)
```

Paths of the output files are written to screen and they can also be retrieved in R. The realm-labelled inspection tree paths are:

```r
res_tact$paths$backbone_biogeo_labels
res_tact$paths$cleaned_tree_biogeo_labels
```

The primary analysis tree is the cleaned TACT tree without realm suffixes should have been written to:

```text
project/results/tact/Formicidae-complete-tact-biogeo-2Aug2026_tacted_cleaned.tre
```

The realm-labelled tree is an inspection tree:

```text
project/results/tact/Formicidae-complete-tact-biogeo-2Aug2026_tacted_cleaned_biogeo_labels.tre
```

Use the realm-labelled tree to visually check whether the biogeography-aware grafting behaved as expected, but use the cleaned tree for downstream comparative analyses unless realm suffixes are explicitly needed.

#### Plot and inspect biogeography-aware placement

To plot the realm-labelled TACT tree use:

```r
t <- ape::read.tree(
  here::here(
    "project", "results", "tact",
    "Formicidae-complete-tact-biogeo-2Aug2026_tacted_cleaned_biogeo_labels.tre"
  )
)

plot_tree_autosize(
  t,
  here::here("project", "results", "tact", "tact_biogeo.pdf"),
  cex = 0.1
)
```

You can also plot the clean final tree without realm suffixes:

```r
t_clean <- ape::read.tree(
  here::here(
    "project", "results", "tact",
    "Formicidae-complete-tact-biogeo-2Aug2026_tacted_cleaned.tre"
  )
)

plot_tree_autosize(
  t_clean,
  here::here("project", "results", "tact", "tact_cleaned.pdf"),
  cex = 0.1
)
```

Once confirmed this worked, we generate 100 replicates of TACT-completed trees:

```r
seeds <- 1:100

for (seed_i in seeds) {
  run_tact_grafting(
    backbone_tree = here::here(
      "project", "results", "grafted",
      "chronosta_gapfilled_final_tree_monophyletic.nwk"
    ),
    taxonomy = here::here(
      "project", "tables",
      "antwiki-valid-species-2Aug2026.txt"
    ),
    out_prefix = here::here(
      "project", "results", "tact",
      sprintf("Formicidae-complete-tact-biogeo-4Aug2026_rep%03d", seed_i)
    ),
    taxonomy_format = "antwiki",
    tact_runner = "docker",
    docker_image = "jonchang/tact",
    seed = seed_i,
    nonmono = "split",
    nonmono_allocation = "random",
    genus_only = "replace_random_species",
    species_code = "keep",
    drop_code_species_after_tact = TRUE,
    enforce_taxonomy_tip_count = FALSE,
    biogeo = TRUE,
    country_realm_map = here::here(
      "project", "tables",
      "country_udvardy_realm.tsv"
    ),
    write_biogeo_labelled_trees = FALSE,
    keep_temp = FALSE
  )
}
```

Then we combined these into a single tree outside of R:
```bash
cat Formicidae-complete-tact-biogeo-4Aug2026_rep*_tacted_cleaned.tre > 100-tact-replicates.tre
```

#### Final output to use in downstream analyses

For most downstream analyses, use this file:

```text
project/results/tact/Formicidae-complete-tact-biogeo-2Aug2026_tacted_cleaned.tre
```

or, within the same R session that ran TACT, use:

```r
res_tact$tree
```

If you want to use the 100 replicate trees, see
```text
project/results/tact/100-tact-replicates.tre
``` 

### Re-running the workflow after adding new data

To update the tree when new source phylogenies become available:

1. Add the new source tree to the appropriate folder under `project/published/` and/or `project/chronosta/source_trees/`.
2. If it should be grafted directly as a clade, add or update its row in:

   ```text
   project/tables/clades-to-graft-clean.tsv
   ```

3. If it should be used only in the Chrono-STA gap-filling step, place it in:

   ```text
   project/chronosta/source_trees/
   ```

4. Re-run the workflow from Step 1 or from the earliest affected step.
5. Compare the new logs and trees.

For reproducible runs, do not use `setwd()` or absolute paths. Use `here::here(...)` throughout, and keep seeds fixed for stochastic steps.