require 'byebug'
require 'csv'
require 'amazing_print'

require_relative 'lib/cgq'


=begin
NOTES

# Another variable -> loci same, loci different (dataset report)
# Identify overlapping LOCI -> the same target and query locus AND a different target and locus
#   Graph of locus matches -> query taxon loci matches target taxon (broad overlap) (less contamination) 
# 
# Query_seq (e.g L420) -> filter on that -> multiple targets 
#     NOT on same plate
#         on same plate
# query loci target multiple taxa 
#
# !! LOCI are sequential, CONTIGS -> could be overlapping, might not be, some might be single locus

# Sum differences of different taxa per LOCI <- good one 
#
#   per loci -> sum different taxa 
#            -> sum of end matches
#            -> sum of start matches
#            -> close matches (on plate)
#            -> far matches (on plate)
#
# predict contaminator ---  target genus is contaminator (flips either way) 
=end 


# Open the data
file = CSV.open( File.expand_path('../data/working/possible_contamination_data.csv', __dir__), col_sep: ',', headers: true)
plates = CSV.open( File.expand_path('../data/working/plates.csv', __dir__), col_sep: ',', headers: true)

# Create a new data object
data = Cgq::Records.new

# Populate the data object
file.each do |r|
  data.all_rows.push Cgq::Row.new(r.to_h)
end
file.close

# Index the records against their plates
current_plate = '0.5' # why? dunno
current_pos = 1
plates.each_with_index do |r|
  h = r.to_h
  if h['AE plate #'] != current_plate
    current_plate =  h['AE plate #']
    current_pos = 1
  end

  h.merge!('position' => current_pos)
  data.plate_rows[h['AE I#']] = h 
  current_pos += 1
end
plates.close

# At this point data are all loaded

=begin

# Do things like:
# puts data.loci.sort.collect{|l| l.join(',')}.join("\n")

# Or:

scores = []
data.all_rows.each do |r|
  s = data.concentration_difference(r)
  print s 
  print ':'
  print data.score_concentration_difference(r)
  print "\n"
  scores.push s.to_f
end 
puts 'min: ' + scores.compact.min.to_s
puts 'max: ' + scores.compact.max.to_s
=end


# Echo a count dataset to the console, pipe (`|`) it to a file
# puts Cgq::Report.count_records_heatmap(data, '1')

# Re-write the genus_ids file
# Cgq::Report.write_genus_ids(data)

# Cgq::Report.write_scores(data)

# Write count heatmaps
# Cgq::Report.count_heatmaps(data)

# puts data.foo

Cgq::Report.write_overlap_loci_by_genera(data)

# puts ap data.locus_overlap_by_i_num 
Cgq::Report.locus_overlap_by_i_num(data)

#puts ap data.overlap_type_per_locus_pair 
Cgq::Report.overlap_type_per_locus_pair(data)






