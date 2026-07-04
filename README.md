# SpaceGroupCohomology — a GAP package

Mod-2 cohomology rings and Lieb–Schultz–Mattis anomaly classes for the
230 crystallographic space groups.

Reference: Chunxiao Liu and Weicheng Ye, *Crystallography, group cohomology,
and Lieb–Schultz–Mattis constraints*, SciPost Phys. **18**, 161 (2025)
([arXiv:2410.03607](https://arxiv.org/abs/2410.03607)).

## Installation

Copy or symlink this directory into a path on your GAP `pkg/` search list,
typically `~/.gap/pkg/SpaceGroupCohomology/`:

```bash
mkdir -p ~/.gap/pkg
ln -s "$(pwd)/SpaceGroupCohomology" ~/.gap/pkg/SpaceGroupCohomology
```

The package depends on [HAP](https://gap-packages.github.io/hap/) (≥ 1.30),
which is loaded automatically.

### Optional: external linear-algebra accelerator

Large computations can offload mod-2 linear algebra to a small bundled C++
tool. Building it is optional — without it the package behaves exactly as
before:

```bash
make -C cpp        # builds bin/sgclinalg (no external dependencies)
make -C cpp test   # self-tests
```

Matrices whose size (rows × columns) exceeds `SGC_LINALG_THRESHOLD`
(default `10^6`) are then handled by `bin/sgclinalg` via a sparse-matrix
file interface; everything smaller stays on GAP's native routines, so
published outputs are unchanged. Set the GAP global
`SGC_LINALG_THRESHOLD := 0;` to force everything through the external tool,
or point the environment variable `SGC_LINALG_PATH` at a binary outside the
package directory.

## Usage

```gap
gap> LoadPackage("SpaceGroupCohomology");
true
gap> SpaceGroupCohomologyRingGapInterface(191);   # any IT number in 1..230
```

That single call prints the mod-2 cohomology ring presentation and the
LSM anomaly classes for the given space group.

## Layout

```
SpaceGroupCohomology/
├── PackageInfo.g       package metadata + HAP dependency
├── init.g              forward declarations of data globals
├── read.g              reads gap/data.gi then gap/functions.gi
├── gap/
│   ├── data.gi         PGGens230, IWP, GENNAMES, funcs230 (230 entries each)
│   └── functions.gi    library of cohomology-ring routines
└── tst/
    └── smoke.tst       package self-test (Test or TestPackage)
```

## Testing

`tst/smoke.tst` does two things in order:

1. A constant-time sanity check that all four data tables have 230 entries.
2. Three real cohomology computations — space groups No. 1, 16, 89 — that
   exercise the full HAP-backed pipeline end-to-end and pin down the
   expected output. Any regression in the cup product, relations, or
   LSM-class code will then surface as a `Test` failure instead of passing
   silently.

```gap
gap> TestPackage("SpaceGroupCohomology");
```

The first item is instantaneous; the three computations take a few seconds
to ~30 s total on a modern laptop, with group 89 (tetragonal `P4/mmm`) the
slowest of the three.

## Reference

If this repository is useful for your research, please consider citing the [Scipost](https://www.scipost.org/SciPostPhys.18.5.161) article:

```bibtex
@Article{10.21468/SciPostPhys.18.5.161,
	title={{Crystallography, group cohomology, and Lieb–Schultz–Mattis constraints}},
	author={Chunxiao Liu and Weicheng Ye},
	journal={SciPost Phys.},
	volume={18},
	pages={161},
	year={2025},
	publisher={SciPost},
	doi={10.21468/SciPostPhys.18.5.161},
	url={https://scipost.org/10.21468/SciPostPhys.18.5.161},
}
```