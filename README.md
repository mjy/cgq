
# C(halcid) G(genome) Q(uality)

Utilities to compute some metrics from a custom set of metadata.

_Very preliminary_

## Code

Written in Ruby

See `src/`

* `go.rb` is a playground that drives generation/experimentation for now 

### Quick start

Clone this repo to your laptop.  Assume you have Ruby installed (pretty safe).

```
cd cgq
bundle install

ruby src/go.rb
```

## Results

All derived data are on `out/`.

### CSV

#### scores.csv
`scores.csv` - summary table of individual metrics, and column that combines some of them

* `overap_type` - A-G 

* `s_locus-difference` (1) same, (0) different
* `s_proportional_length` (1) length match bp / tlength > curoff, (0) < 
* `s_taxon_difference` (1) different families, (0) same families
* `s_plate_difference` (1) same plate, (0) different plate
* `s_concentration_difference` (1) absolute (query qbit - target qbit) > 3, (0) different plates OR < 3 OR no qbit values
* `s_proportional_difference` (1) %similar == 100.00, (0) !=

#### overlap_by_loci_by_genus
* `overlap_by_loci_by_genus.csv` - unique overlap types per unique loci pair per unique genus pair

#### locus_overlap_by_i_num
_Idea was 3 or more indicates an issue._
* `locus_overlap_by_i_num.csv` - unique target loci per query locus/i# pair 

### HTML
To visualize a report clone this repo to your desktop, then click on an HTML file.

