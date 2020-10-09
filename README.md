
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
* `scores.csv` - summary table of individual metrics, and column that combines some of them
* `overlap_by_loci_by_genus.csv` - unique overlap types per unique loci pair per unique genus pair

### HTML
To visualize a report clone this repo to your desktop, then click on an HTML file.

