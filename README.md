# PDL::IO::NYHead

Read the **New York Head** forward model (`sa_nyhead.mat`, MATLAB v7.3 / HDF5)
into [PDL](https://pdl.perl.org). A focused reader that exposes electrodes,
leadfield, cortical/scalp surfaces, the Harvard-Oxford atlas, coordinate
transforms, and nearest-vertex lookup with a small, PDL-native API.

This distribution was split out of the `PDL-EEG` project so it can be reused
independently. It has no dependency on `PDL-EEG`.

## Why not PDL::IO::Matlab?

`PDL::IO::Matlab` (matio-backed) can open MAT v7.3 files, but its Perl API only
surfaces real numeric arrays: it does not traverse the nested `sa` struct and
cannot decode MATLAB cell-of-strings (the electrode and atlas names, stored as
HDF5 object references). This distribution reads numeric data with
`PDL::IO::HDF5` and recovers the string fields via `h5dump`
(`PDL::IO::NYHead::H5Dump`), which is a genuinely general MATLAB-v7.3
cell-of-strings helper.

## Requirements

- PDL, PDL::IO::HDF5
- `h5dump` (from the HDF5 tools) on `PATH` — **only** needed for the label
  fields (`electrode_labels`, `ho_labels`). Everything numeric works without it.
  MacPorts: `port install hdf5`.

## Install

    perl Makefile.PL
    make
    make test
    make install

## Synopsis

    use PDL::IO::NYHead;
    my $ny = PDL::IO::NYHead->new('sa_nyhead.mat');

    my $labels = $ny->electrode_labels;       # 231 names
    my $pos    = $ny->electrode_pos;          # (231,3) MNI mm
    my $L      = $ny->leadfield;              # (231,74382) normal-oriented CAR

    my $ctx    = $ny->surface('cortex75K');   # { vc, tri, tri0 }
    my $lo     = $ny->cortex('cortex2K');     # reconstructed from cortex75K

    my $areas  = $ny->ho_labels;              # 97 Harvard-Oxford area names
    print $ny->area_of_vertex(1234), "\n";    # e.g. 'Left Frontal Pole'

    my $iv = $ny->nearest_vertex($pos, 'head'); # each electrode -> scalp vertex

    # export a surface to OBJ / USDA
    $ny->write_mesh('cortex75K', 'obj',  'cortex.obj');
    $ny->write_mesh('head',      'usda', 'head.usda');

## Axis convention

h5dump reports shapes `{A,B,...}`; PDL reads them reversed as `(...,B,A)`. The
New York Head stores geometry as `{3,N}`, so `vc` / `normals` / `locs` / `tri`
come back as `(N,3)` / `(Nf,3)` directly (no transpose). `tri` is 1-based
(kept for OBJ export; `tri0` is the 0-based copy). `V_fem_normal` is
`(electrode, source)`.

## See also

`examples/check_nyhead.pl` — runs the full reader against a real
`sa_nyhead.mat` and prints a summary with PASS/FAIL structural checks.
