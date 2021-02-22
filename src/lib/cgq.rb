require 'csv'
require 'URI'
require 'net/http'
require 'JSON'
require 'fileutils'

=begin
contamination
query_genus
target_genus
nil
I# query
I# target
query_seq
target_seq
nil
%similar
length of match (bp)
tlength
qlength
mismatch
gapopen
qstart
qend
nil
tstart
tend
nil
evalue
note     # "withinloci" -> ignore
=end


=begin
# Target stats

taxon differece (0/1) 				0=same
plate difference <10 (vary) (1/1)         	0=far
locus same/different (0/1)			0=different
percent difference (100, 98-99, 95-97)	0=95-99    	1=100
proportional difference (<20%) (0/1)		0=<20

target/query query/target - realistically - look at half -> anything flipped should be deleted!!

# Another variable -> loci same, loci different (dataset report)
# Identify overlapping LOCI -> the ame target and query locus AND a different target and logus
#   Graph of locus matches -> query taxon loci matches target taxon (broad overlap) (less contanmination) 
# 
# Query_seq (e.g L420) -> filter on that -> multiple targets 
#     NOT on same plate
#         on same plate
# query loci target multiple taxa 
# 1 0 ??
#
# !! LOCI are seqeuntial !! CONTIGS -> could be overlapping, might not be, some might be single locus
# 

# Junxia had a lot of loci that were deleted (within taxa) - remove both
#   contaminator, contaminatee <- predicted          

# less than 100% 
# More missmatches 
#
# SUM differences of different taxa per LOCI <- good one 
#
#   per loci -> sum different taxa 
#            -> sum of end matches
#            -> sum of start matches
#            -> close matches (on plate)
#            -> far matches (on plate)
#
#  within a column -> 1-12 on a plate positions # x, y # 
#     
#     same plate
#     same column
#     adjacent 
#
# !! smaller left # -> more likely contaminants
#
# thing with the lowest concentration has a higher
#   
#    predict contaminator ---  target genus is contaminator (flips either way) 
=end 


module Cgq
  # Local
  # UCD_API = 'http://127.0.0.1:3000/api/v1/taxon_names?'
  # PROJECT_TOKEN = 'Q2eJtYrqIfa9hXTrAvkVnQ'

  # Remote
  UCD_API = 'https://sfg.taxonworks.org/api/v1/taxon_names?'
  PROJECT_TOKEN = 'adhBi59dc13U7RxbgNE5HQ' 

  # 96 - top left to bottom right
  PLATE_WIDTH = 8  # (columns)
  PLATE_HEIGHT = 12  

  # TODO: not used
  DISTANCE_CUTTOFF = 2
end

require_relative 'cgq/records'
require_relative 'cgq/row'
require_relative 'cgq/report'
