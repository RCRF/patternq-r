# patternq (R)

Query and analysis tools for the Pattern Data Commons. One package replacing
datalogr, wick and (incrementally) luminance: Datalog queries over the Pattern
Data Commons query service, canned queries for every kind of data in a dataset
database, and plotly plots.

The patternq family: [patternq](https://github.com/RCRF/patternq) (Python), [patternq-r](https://github.com/RCRF/patternq-r) (R), [patternq-clj](https://github.com/RCRF/patternq-clj) (Clojure) and [PatternQ.jl](https://github.com/RCRF/PatternQ.jl) (Julia) share one function catalog, the same result columns and the same plots.

## Install

```r
# install.packages("remotes")
remotes::install_github("RCRF/patternq-r")
```

On macOS with a conda/miniforge Python on `PATH`, build the `curl`/`xml2`
dependencies with miniforge removed from `PATH`, or their shared libraries
link against conda's `libz` and fail to load.

## Configure

```sh
export PATTERNQ_ENDPOINT=https://data-commons.rcrf-dev.org
export PATTERNQ_API_KEY=...   # user settings page of the Pattern Data Commons dashboard
```

or in a session: `set_query_server()`, `set_token()`.

## Quick start

Every dataset is its own database. Resolve a dataset name to its current
database and set it as the session default, or pass `db =` to any function.

```r
library(patternq)
list_datasets()                          # dataset, db, counts, assays, tags
set_db(resolve_db("prince-2022"))

dataset_summary()                        # assays and measurement sets
measurement_types("PICI CyTOF Immune Profiling")
subjects(); samples(); timepoints()
subject_outcomes()                       # BOR / PFS / OS per subject

cy <- measurements("percent-of-parent", "PICI CyTOF Immune Profiling")
cy <- add_sample_context(cy, include.outcomes = TRUE)
plot_by_timepoint(cy[cy$cell_population == "CD8 T cells", ], group = "bor")

plot_survival(subject_outcomes(), "os", "os_event", group = "bor")
plot_mutation_landscape(variants(), n.genes = 20)

# CNV queries always take a subset: genes, samples or subjects
cnv_segments(db = resolve_db("H37004"), genes = c("TP53", "BAP1"))
```

### A sample against a reference cohort

```r
h <- resolve_db("H37001"); uvm <- resolve_db("tcga-uvm")
cmp <- compare_to_cohort("H37001-003", db = h, cohort.db = uvm, measurement = "tpm")
top_by_zscore(cmp, 25)
plot_zscores(cmp)
examine_geneset(c("BAP1", "GNAQ", "PRAME", "PMEL"), samples = "H37001-003", db = h, cohort.dbs = uvm)

ch <- compare_samples("H37001-003", "H37001-001", db = h)   # two-sample change
plot_ma(ch); plot_fold_change(ch)
```

Only compare measurements that are comparable across the two databases (same
units and normalization; TPM is the usual common ground), and check
`cohort_observed`: a usually-expressed gene with few stored cohort values
points to an import or annotation problem in the cohort.

### Writing queries

Queries are the JSON form of Datomic Datalog, written as R data:

```r
q <- dq(find = c("?sample-id", "?vaf"),
        where = list(c("?m", ":measurement/vaf", "?vaf"),
                     c("?m", ":measurement/sample", "?s"),
                     c("?s", ":sample/id", "?sample-id"),
                     list(list(">", "?vaf", 0.3))))
do_query(q)
```

or with the datalogr DSL:

```r
genes <- c("TP53", "KRAS")
do_query(query(find(?sample-id, ?hgnc, ?vaf),
               where(d(?m, measurement/vaf, ?vaf),
                     d(?m, measurement/variant, ?v),
                     d(?v, variant/gene, ?g),
                     d(?g, gene/hgnc-symbol, ?hgnc),
                     d(?m, measurement/sample, ?s),
                     d(?s, sample/id, ?sample-id)),
               args(?hgnc <- genes)))
```

Every canned query has a `*_query()` companion returning the query as data.
Results carry `provenance()`: database, basis t and query time.

`cache = TRUE` (the default) uses the service's S3 result cache; `cache =
FALSE` returns results inline and skips it; `refresh.cache = TRUE` recomputes.

## Tests

```sh
R CMD INSTALL .
Rscript -e 'testthat::test_local()'   # live tests run when PATTERNQ_API_KEY is set
```

## Examples

`examples/` has a tutorial and a PRINCE trial demo as Quarto notebooks
(`quarto render examples/tutorial.qmd`); the same notebooks exist for every
language in the family.

## Coming from wick / datalogr

- No `dataset.name` arguments: pass the database (`db =`) instead.
- `get_all_x()` / `all_x()` are now `x()` / `x_query()`.
- datalogr's `with()` is `with_vars()` so `base::with` stays usable.
- Enum values come back as names (`"female"`, `"high"`), no `resolve_db_idents()` needed.

## Contributing

patternq is one library in four languages: [patternq](https://github.com/RCRF/patternq)
(Python), [patternq-r](https://github.com/RCRF/patternq-r) (R),
[patternq-clj](https://github.com/RCRF/patternq-clj) (Clojure) and
[PatternQ.jl](https://github.com/RCRF/PatternQ.jl) (Julia). All four are generated
from a common source and published here, so this repository does not accept pull
requests.

Please report bugs and feature requests as
[issues](https://github.com/RCRF/patternq-r/issues). Code is welcome in an issue: a
minimal example (the call or query, the dataset, what you expected and what you
got), or a proposed change as a snippet, is the most useful way to suggest one.
