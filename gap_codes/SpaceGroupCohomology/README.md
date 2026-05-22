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

```gap
gap> TestPackage("SpaceGroupCohomology");
```

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