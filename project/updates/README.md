# Ant Macrophylogeny Update List

This directory tracks changes that should be considered for future releases of
the ant macrophylogeny.

The main file is:

`macrophylogeny_updates.tsv`

It records new phylogenetic studies, questionable taxon identities or
placements, corrections to donor and final trees, and updates to the AntWiki
taxonomy.

## General rules

- Use one row for one independently understandable issue.
- Assign every item a permanent identifier such as `MPM-001`.
- Never reuse or renumber an identifier.
- Use repository-relative paths, not paths specific to one computer.
- Separate multiple taxa, sources, or paths with semicolons.
- Leave unknown fields empty rather than guessing.
- Do not delete completed or rejected rows. They document previous decisions.
- Save the file as UTF-8 tab-separated text.
- Check the Git diff after editing the file in Excel or another spreadsheet
  program.

## Columns

### `id`

A permanent identifier for the item.

Examples:

`MPM-001`

`MPM-002`

Use the next available number. Refer to this identifier in commit messages,
GitHub issues, and release notes when useful.

### `status`

The current state of the item.

Allowed values:

- `to_do`: Recorded but not yet being worked on.
- `in_progress`: Someone is actively investigating or implementing it.
- `blocked`: Work cannot continue until evidence or a source becomes available.
- `done`: The issue has been resolved, relevant files have been updated, and
  the result has been checked.
- `not_doing`: The proposal was evaluated and intentionally rejected.

An item should not be marked `done` merely because a file was edited. The
affected workflow should also be rerun or otherwise checked.

### `priority`

How urgently the item should be addressed.

Allowed values:

- `high`: Important correction or source that should normally be handled before
  the next release.
- `medium`: Useful improvement that is not release-critical.
- `low`: Minor cleanup, documentation, or low-impact improvement.

### `category`

The main kind of update.

Allowed values:

- `new_phylogeny`: A newly published phylogenetic study should be evaluated or
  incorporated.
- `taxon_check`: The identity, name, or position of one or more taxa needs to be
  investigated.
- `taxonomy_update`: A new AntWiki taxonomy snapshot or other broad taxonomic
  update should be incorporated.
- `tree_correction`: A known problem in a donor tree, backbone, workflow table,
  or final tree should be corrected.
- `other`: The item does not fit the categories above.

### `taxa`

The species, genus, or larger group affected by the item.

Examples:

`Camponotus_laevigatus`

`Pheidole`

`Allomerus;Wasmannia`

`Formicidae`

Use underscore-separated species names. Separate multiple taxa with
semicolons. If the identity is uncertain, preserve the original source label
rather than assigning a valid name without evidence.

### `summary`

A concise explanation of the issue and why it matters.

Good:

`The source and identity of the terminal currently identified as
Camponotus_laevigatus need verification.`

Too vague:

`Check Camponotus.`

State uncertainty explicitly. Do not describe a suspected error as confirmed
unless the evidence supports it.

### `source`

The publication, taxonomy authority, repository, correspondence, or other
evidence associated with the item.

Examples:

`Ward et al. 2025`

`AntWiki valid species list retrieved 2027-01-15`

`Author et al. 2026; 10.xxxx/xxxxx`

Include a DOI or stable URL when useful. For a taxonomy update, record the
authority and exact retrieval date.

### `files_affected`

Repository-relative files or directories likely to be inspected or changed.

Examples:

`project/published/camponotus/ward2025/`

`project/chronosta/source_trees/camp.nwk`

`project/tables/clades-to-graft-clean.tsv`

Separate multiple paths with semicolons. Do not use machine-specific paths
such as `/home/name/...` or `C:\Users\name\...`.

### `action`

A concrete description of what needs to be done.

For a new phylogeny, this might include:

1. obtain the original tree and supplement;
2. determine whether it is a phylogram or chronogram;
3. check rooting and branch-length units;
4. reconcile source labels with valid names;
5. decide whether to use tip grafting, clade grafting, or Chrono-STA;
6. rerun the affected downstream stages.

For a taxon check, this might include:

1. identify the original source terminal;
2. inspect voucher, accession, and locality information;
3. compare the name with current taxonomy;
4. rename, retain as an operational label, reposition, or remove the terminal;
5. rerun the affected stage and inspect the result.

For a taxonomy update, this might include:

1. obtain a dated AntWiki taxonomy snapshot;
2. compare it with the previous snapshot;
3. review added and removed species;
4. review synonymies and changed genus combinations;
5. update the active taxonomy file;
6. rerun TACT and review mismatch reports.

### `outcome`

What was ultimately decided and changed.

Leave this blank while the item is open.

A completed outcome should state:

- what evidence was found;
- what decision was made;
- which label, source, or placement changed;
- whether the affected analysis was rerun;
- any important remaining limitation.

Examples:

`The terminal was traced to voucher CASENT0123456 and renamed to
Camponotus_example. The Camponotus donor tree and Chrono-STA stage were rerun.`

`The terminal could not be assigned confidently to a valid species and was
removed from the donor tree.`

`The new study was not incorporated because no machine-readable tree and no
reproducible topology were available.`

### `date_added`

The date the item was added, using:

`YYYY-MM-DD`

Do not change this date when the status changes.

### `date_completed`

The date the item became `done` or `not_doing`, using:

`YYYY-MM-DD`

Leave it empty for open items.

## Working with the registry

When adding an item:

1. Search for an existing row about the same source or taxon.
2. Assign the next `MPM-###` identifier.
3. Complete the issue, action, source, and affected-path fields as fully as
   possible.
4. Set the status to `to_do`, `in_progress`, or `blocked`.
5. Commit the new row.

When completing an item:

1. Make the required changes.
2. Rerun the affected part of the workflow.
3. Inspect the relevant logs, mappings, and resulting tree.
4. Fill in `outcome`.
5. Change `status` to `done` or `not_doing`.
6. Enter `date_completed`.
7. Commit the registry and associated files together when practical.

## Relationship to other project files

- `TODO.md` tracks software and repository-development tasks.
- `project/updates/macrophylogeny_updates.tsv` tracks scientific and data
  updates to the ant tree.
- `project/tables/` tells the workflow what operations to perform.
- `project/published/` stores and documents published phylogenetic sources.
- `project/chronosta/source_trees/` contains active Chrono-STA donor trees.
- `project/results/` contains generated trees, logs, and validation output.
- Release notes summarize completed changes included in a particular published
  version of the ant macrophylogeny.