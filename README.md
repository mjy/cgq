
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
* `scores.csv` - summary table of individual metrics, and column that combines some of them
* `ucd_genus_id.csv` - result of querying the UCD by string name
