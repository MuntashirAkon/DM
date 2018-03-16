# DM
A tool to generate distance matrix and sorted species relations from the given minimal or relative absent words.

## Installation
Download/clone the repository and set root directory as the current directory.
### Using cmake
```bash
cmake ./ && make
```
### Manual Installation
```bash
g++ -o dm main.cpp Set.cpp Process.cpp
```

## Usage
`dm` takes 4 arguments (or, at least three arguments) namely *Absent Word Type*, *Dissimilarity Index*,
*Source Directory* and (optional) *Target Directory*. If no *Target Directory* is supplied, *Source Directory* is
assumed to be the *Target Directory*. You can also find this help by running `dm` without any argument.

```bash
dm <Absent Word Type> <Dissimilarity Index> <Source Directory> [<Target Directory>]
```

If all the conditions are fulfilled, output files will be generated in the *Target Directory*.

### Absent Word Types
| Argument Name | Abbreviation
| ------------- | ------------
| MAW           | Minimal Absent Words
| RAW           | Relative Absent Words

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
The source directory should contain a `SpeciesFull.txt` where each line should contain the names of the species
(`CRLF` is not supported right now, so don't use Notepad to edit `SpeciesFull.txt`. Use Notepad++ on Windows.),
and DNA/RNA sequences of the species in the format described below. **No description header should be included.**

**For Minimal Absent Words:**
The MAW sequence of a species should be contained in the `<species name>.maw.txt` file, e.g. human.maw.txt.

**For Relative Absent Words:**
The folders should be named after the referenced species name (let's suppose `human`). Each folder contains files
which are to be named after the referenced species name and the species names it was compared to, joined by an underscore
and the extension should be `.raw.txt`, e.g. `human_goat.raw.txt`.

> See the `Examples` folder for details. (RC = Reverse Complement)

### Outputs
- `SpeciesRelation.txt` : Contains sorted species relations (Readable)
- `SpeciesRelation.json` : Contains sorted species relations (Programmable)
- `DistanceMatrix.txt` : A plain text of Distance Matrix (Programmable & Readable)
- `Output.txt` : A formatted Distance Matrix to be used with `Match7`

> The `Match7` folder contains the modified `Match7` file. The only argument that it takes is the directory containing
> `DistanceMatrix.txt` and `SpeciesFull.txt` files.
