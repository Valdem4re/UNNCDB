# UNNCDB

UNNCDB builds a PostgreSQL database from CIF files and prepares datasets from the
structures stored in it. It has a C++ importer and a Python module that you can
use directly in a notebook.

The current dataset workflow focuses on centrosymmetry prediction: unit-cell
parameters and calculated X-ray intensities are the features, and the presence
of an inversion centre is the target label.

## How it works

```text
CIF files → PostgreSQL → structure selection → diffraction calculation → pandas / CSV
```

The importer keeps the original CIF blocks together with parsed cell parameters,
atoms, symmetry operations and refinement data. Problems found during import are
recorded in the database.

The Python module reads these structures, applies selection rules, checks the
atomic models and calculates features on a shared reflection grid. You can work
with the result in memory or save it as one CSV for later experiments.

## Build and import

You will need:

- A C++20 compiler and CMake 3.25 or later.
- PostgreSQL, its command-line tools and libpq development files.
- OpenSSL development files (`libcrypto`).

CMake downloads Gemmi during configuration. Run these commands from the project
root:

```sh
cmake -S parser -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

With PostgreSQL running, create an empty database and load a CIF file:

```sh
createdb unncdb
psql -X -v ON_ERROR_STOP=1 -d unncdb -f unncdb_schema.sql

build/cif_parser structures.cif "dbname=unncdb"
```

The schema file creates all tables, indexes and views in one transaction. It is
intended for an empty database; an already populated database does not need to
be recreated.

Replace `structures.cif` with your file. The importer reads uncompressed CIF 1.x
files, including files with multiple data blocks. Connection settings follow the
usual libpq format: `host`, `port`, `user` and `.pgpass` are supported.

Data is read block by block and written in batches. The file checksum is checked
before and during import, so the first progress message may take a while to
appear. Keep the input file unchanged during loading.

If the process is interrupted, run the same command again. Already committed
blocks are reused. Progress messages show how many blocks have been read, stored
and reused, along with the time spent parsing and writing.

Batch limits are configurable:

```sh
build/cif_parser structures.cif "dbname=unncdb" \
  --batch-size 1000 \
  --batch-bytes 16777216 \
  --batch-rows 25000
```

These are the defaults. A batch is written when any limit is reached. The byte
limit counts source data; parsed rows and write buffers also take memory.

To check a file without connecting to the database:

```sh
build/cif_parser structures.cif --validate-only --max-blocks 2000
```

## Use the Python module

Python 3.10 or later is required. Install the package into the environment used
by your notebook:

```sh
python -m pip install -e .
```

In Jupyter, use `%pip install -e "/path/to/unncdb"`. The module can then be
imported from any notebook directory.

```python
from unncdb_dataset import Database, Filters, ReflectionGrid, Calculation

db = Database("dbname=unncdb")

dataset = db.build(
    n=1000,
    balance=True,
    seed=42,
    filters=Filters(r_max=0.1),
    grid=ReflectionGrid.cube(4),
    calculation=Calculation(wavelength=0.71073, anomalous=True),
)

X, y = dataset.X, dataset.y
print(dataset.report)
dataset.frame.head()
```

This requests 500 centrosymmetric and 500 noncentrosymmetric structures after
all checks. The database connection is read-only, and candidates are fetched in
batches.

| Attribute | Contents |
| --- | --- |
| `dataset.frame` | Full table with features, labels and source metadata |
| `dataset.X` | Cell parameters and reflection features |
| `dataset.y` | Labels: `0` for noncentrosymmetric, `1` for centrosymmetric |
| `dataset.groups` | CSD refcode families for grouped splits |
| `dataset.mask` | Which reflections are available in the calculation range |
| `dataset.rejections` | Skipped candidates and the reasons they were skipped |
| `dataset.report` | Sample size, class counts and selection statistics |

Use `dataset.X` for training. The full table also contains the known space group
and other metadata that should stay outside the model input.

### Selection rules

The default filters require `R1_gt <= 0.1`, a known target label and no import
errors. They reject known disorder, ambiguous labels such as `Ca1?`, unknown
coordinates and partial or explicitly unknown occupancies.

Coordinates are wrapped into `[0, 1)`. Symmetry is expanded using the operations
or Hall symbol from the CIF. Repeated atoms of the same element at the same
periodic position are removed. If different elements share a position, the whole
structure is rejected.

Missing disorder information does not establish that a structure is ordered.
The checks use the descriptions, atom-site fields and coordinates available in
the source file.

You can change the filters without editing the module:

```python
filters = Filters(
    r_max=0.08,
    required_elements=("Ca",),
    excluded_elements=("Hg",),
    max_volume=20000,
)

dataset = db.build(
    n=1000,
    filters=filters,
    seed=17,
    max_per_space_group=200,
    max_candidates=20000,
)
```

To inspect candidate counts before the atom checks and deduplication:

```python
audit = db.audit(filters)
print(audit.groupby("is_centrosymmetric")["candidates"].sum())
```

If too few structures pass, `build()` raises `InsufficientDataError`. Pass
`allow_partial=True` to inspect a smaller result. The report shows whether the
requested size and balance were reached.

By default, only one structure per CSD refcode family is kept: `ABCDEF`,
`ABCDEF01` and `ABCDEF02` belong to the same family. Repeated CIF blocks and
matching atomic models in the same cell are also removed. Refcode grouping does
not identify every chemically related structure, and model matching does not
cover arbitrary changes of basis or origin.

Set `one_per_family=False` to keep multiple family members, then use
`dataset.groups` when splitting the data. Splits are left to the training code,
so one dataset can be used with different splits and cross-validation schemes.
Fit scaling and imputation on the training subset only.

### Calculated features

`ReflectionGrid.cube(4)` creates 728 reflections with indices from −4 to 4,
excluding `000`. Every structure uses the same column order. Friedel pairs are
kept separately, systematic absences stay as zeros, and reflections are not
merged using the target space group.

The calculation uses IT92 atomic scattering factors, with anomalous corrections
enabled by default at the chosen wavelength. The default features are intensities
divided by their mean over the available reflections of each structure. This is
`I / mean(I)`, not crystallographic E² normalization.

Use `Calculation(normalization="none")` for raw intensities, or
`Calculation(quantity="amplitude")` for amplitudes. Thermal motion is disabled by
default. `Calculation(thermal="isotropic")` uses isotropic U or B values and
requires them for every atom. Full anisotropic scattering is not implemented.
These are calculated diffraction data from experimental structure models.

A reflection outside the chosen wavelength and angular range is stored as `NaN`.
Its column stays in place, and `dataset.mask` distinguishes it from zero intensity.

### Save and load

```python
from unncdb_dataset import Dataset

dataset.save_csv("datasets/baseline.csv")

loaded = Dataset.load_csv("datasets/baseline.csv")
X, y = loaded.X, loaded.y
```

This writes `baseline.csv` and `baseline.manifest.json`. The CSV has one row per
structure. The JSON stores the shared `hkl` grid, feature order, filters,
calculation settings, source information and package versions. The CSV checksum
is checked on loading.

Saving is optional. Existing files are kept unless `overwrite=True` is passed.
The same seed and settings reproduce a selection while the source records and
module version stay the same.

## Project layout

| Path | Purpose |
| --- | --- |
| `parser/` | C++ CIF importer |
| `unncdb_schema.sql` | Complete PostgreSQL schema |
| `python/unncdb_dataset/` | Dataset module |
| `examples/` | Notebook examples |
| `tests/` | Importer and dataset tests |

The database stores original CIF bytes as well as parsed structures, atoms,
symmetry operations, displacement parameters and import diagnostics. Fields
without a dedicated SQL column can still be read from the saved CIF block.

## Tests

```sh
ctest --test-dir build --output-on-failure
python -m unittest discover -s tests -p test_dataset.py -v
```

The tests cover atom cleanup, diffraction calculations, sampling and CSV export.
Comparison with the calculation notebook also uses `pymatgen` and `tabulate`.

Set `UNNCDB_DATASET_TEST_DSN='dbname=unncdb'` to include the read-only dataset test
against a populated database. Importer integration tests are described in
[tests/integration.py](tests/integration.py); they need a separate PostgreSQL test
server and create their own throwaway databases.

## License

[MIT](LICENSE).
