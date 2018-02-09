# DM
A tool to generate distance matrix and sorted species relations from the given minimal or relative absent words.

## Installation
Download the repository and then, in the root directory, run:
```bash
c++ -o dm main.cpp Set.cpp
```

## Usage
`dm` takes 4 arguments (or, at least three arguments) namely *Absent Word Type*, *Dissimilarity Index*, *Source Directory* and (optional) *Target Directory*. If no *Target Directory* is supplied, *Source Directory* is assumed as the *Target Directory*. You can also find this simply by running `dm` without any argument.
```
dm <Absent Word Type> <Dissimilarity Index> <Source Directory> [<Target Directory>]
```

If all the conditions are fullfilled, output files are generated in the *Target Directory*.

### Absent Word Types
| Argument Name | Abbreviation          |
| ------------- | --------------------- |
| MAW           | Minimal Absent Words  |
| RAW           | Relative Absent Words |

### Dissimilarity Indexes
For Minimal Absent Words:
| Argument Name     | Abbreviation
| ----------------- | ------------
| MAW_LWI_SDIFF     | Length weighted index of symmetric difference 
| MAW_LWI_INTERSECT | Length weighted index of intersection
| MAW_GCC_SDIFF     | GC content of symmetric difference
| MAW_GCC_INTERSECT | GC content of intersection
| MAW_JD            | Jaccard distance
| MAW_TVD           | Total variation distance

For Relative Absent Words:
| Argument Name | Abbreviation
| ------------- | ------------
| RAW_LWI       | Length weighted index
| RAW_GCC       | GC content

### Source Directory format
The source directory should contain a `SpeciesFull.txt` where
each line should contain the names of the species
(`CRLF` is not supported right now, so
don't use Notepad to edit `SpeciesFull.txt`.
Use Notepad++ on Windows.),
and DNA/RNA sequences of the species in the format
described below. The species header is to be ignored.

For Minimal Absent Words:
The sequence of a species name should contain in the `<species name>.maw.txt` file,
e.g. human.maw.txt.

For Relative Absent Words:
There has to be folders named after the species name (let's suppose `human`).
Each folder contains files which is to be named
as the base species name and the species named it was compared to
joined by an underscore and the extension is to be `.raw.txt`, e.g. `human_goat.raw.txt`

> See the `Examples` folder for details. (RC = Reverse Complement)

### Outputs
- `SpeciesRelation.txt` : Contains sorted species relations (Readable)
- `SpeciesRelation.json` : Contains sorted species relations (Programmable)
- `DistanceMatrix.txt` : A plain text of Distance Matrix (Programmable & Readable)
- `Output.txt` : A formatted Distance Matrix to be used with `Match7`

> The `Match7` folder contains the modified `Match7` where the contents of the
> `Output.txt` can be supplied as an array.
sl