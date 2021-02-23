# C(halcid) G(genome) Q(uality)

Utilities to compute some metrics from a custom set of metadata.

## Download

* Click the green 'Code' option on this very page (above to right)
* Click the "Download as zipfile" option

## Materials and Methods

The goal is to use an algorithm to identify sequences to eliminate because they are predicted to be contaminated.  
The core raw data are [possible_contaminations_supercontig_95p.csv](data/working/possible_contaminations_supercontig_95p.csv), a spreadsheet indexing query/target pairs. These data are cross-linked to a [file containing qbit and plate metadata](data/working/plates.csv). 

_Resolving nomenclature_. Prior to the parameterized analysis two additional files are generated by querying a [TaxonWorks API](https://sfg.taxonworks.org/api/v1/) that exposes the Universal Chalcidoidea Database data. The Genus names observed in `possible_contaminations_supercontig_95p.csv` are sent to the API using the parameters `&nomenclature_group=Genus&exact=true&name=#{genus}`, where `#{genus}` is the name of the genus (the project token is also sent). A genus name is matched to the database only if the database record contains and author and year (preventing homonymous, and incomplete data in the database creating an incorrect match). The matches are recorded as a [file](data/working/dervied/ucd_genus_ids.csv). A second pass is used to determine the family of the genus in question, the parameter queries are `&ancestors=true&nomenclature_group=FamilyGroup::Family&taxon_name_id[]=#{taxon_name_id}`, where `#{taxon_name_id}` is the id of the genus. These data are recorded as a [file](data/working/derived/family_metadata.json).

_Contamination metric_. A metric (`composite_score` in [scores.csv](out/csv/scores.csv) is calculated by summing individual scores (columns prefixed with `s_`), each score compares two (or more) properties of the target and query sequences. The scores are:

1. `Locus difference ("s_locus_overlap")` - The nature of the overlap of query and target sequences is calculated (see Addendum, `Locus overlap types`).  If the type is D, E, or F, then the score is 1, otherwise it is 0.  _Likely uninformative as all loci are currently the same!!_
2. `Same locus` ("s_same_locus") - If the locus is identical in query and target, then the score is 1, otherwise the score is 0.
3. `Proportional length difference ("s_proportional_length")` - The length of the of the match (`length_of_match`) is divided by the length of the target (`target_length`). Both lengths are in base pairs.  If the result is > 0.95 then the score is 1, otherwise it is 0.  
4. `Taxon difference ("s_taxon_difference")` - The family and subfamily names of the query and target are compared.  If both family and subfamily are same values then the score is 0.  If either subfamily or family is different the score is 1. Records like "Pteromalidae-Pteromilinae" and "Pteromalidae-blank" are considered _different_.  If either genus name can not be matched to a family then the score is 0.
5. `Plate similarity ("s_plate_similarity")` - The plates of the query and target are compared. If the plate number is the same then the score is 1.  If the plate number is different then the score is 0.  If the plate can not be determined the score is 0.
6. `Concentration ratio difference ("s_concentration_ratio")` - If the `Plate difference` is 0, then the score is zero.  If the ratio of the smallest qbit concetration to the largest is less than 0.3 then the score is 1.  If it is greater than 0.3 then the score is 0.  If the score can not be calculated then the score is 0.
7. `Proportional difference ("s_proportional_difference")` - If the '%similarity' field in the original data is == '100.00' then the score is 1.  If it is some other value (including not being provided) then the score is 0.
8. `Plate column identity ("s_column_identity")` - If the query and target are on the same plate, and they share the same column (compare `query_plate_x` with `target_plate_x`), then the score is 1, otherwise it is 0. 

_Flagging sequences as contaminated_. To determined whether sequence should be eliminated due to contamination the contamination metric is compared to the `Concentration ratio difference` in the following way:
1. A range of Contamination metric scores is defined (`composite_cutoff` in code), by default this is [3,4,5].
2. A concentration ratio cutoff is defined (`concentration_cutoff` in code), by default this is 0.3.
3. If a row has as contamination metric in the range of the composite cutoff the test continues, otherwise that rows has no contaminants.
4. If the concentration ratio (`s_concentration_ratio`) of the row is less than the `concentration_cutoff` then the test continues, otherwise the row has no contiminants.
5. The contaminated sequence is the sequences with the smaller qbit, concentration.  The other sequence is not contaminated.

_TODO: We need to resolve the non-indepdence of this comparison as the `Concentration ratio difference` cutoff is used in 2 places._

## Results

All derived data and reports are in `out/`.

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

### CSV
* `scores.csv` - summary table of individual metrics, and column that combines some of them
* `overlap_by_loci_by_genus.csv` - unique overlap types per unique loci pair per unique genus pair

_Idea was 3 or more indicates an issue._
* `locus_overlap_by_i_num.csv` - unique target loci per query locus/i# pair 

### HTML
To visualize a report clone this repo to your desktop, then click on an HTML file.

## Addendum 

### Locus overlap types

##### Columns
* `overap_type`: (`-` query, `=` target, `|` center of target)
```
A -----
     ==|==

B ----- 
   ==|==

C ------- 
   ==|==

D   --- 
   ====|====

E   ----- 
   ===|===

F       --- 
   ====|====

G   ------
   ==|== 

H     -----
   ==|== 
``` 
