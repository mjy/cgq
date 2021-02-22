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

module Cgq
  class Row
 
   # Default classification of overlaps that 
    # might indicate contamination
    # See `score_overlap`  
   BAD_OVERLAPS = %w{D E F}

    attr_accessor :d # the hash

    def initialize(h)
      @d = h
    end

    def genus_pair
      [d['query_genus'], d['target_genus']]
    end

    def locus_pair
      [query_locus, target_locus]
    end

    def query_locus
      d['query_seq'].split('_').first
    end

    def target_locus
      d['target_seq'].split('_').first
    end

    def sequence_middle(kind = 'query_seq')
      d['qlength'].to_f / 2
    end

    # @return String
    #   A-H
    #
    # Characterize target/query relationship
    # with respect to overlaping regions.
    #
    # A -----
    #      ==|==
    #
    # B ----- 
    #    ==|==
    #
    # C ------- 
    #    ==|==
    # 
    # D   --- 
    #    ====|====
    #
    # E   ----- 
    #    ===|===
    # 
    # F       --- 
    #    ====|====
    # 
    # G   ------
    #    ==|== 
    #
    # H     -----
    #    ==|== 
    #
    def overlap_type
      m1 = sequence_middle('query_seq')
      m2 = sequence_middle('target_seq')

      ml = d['length of match (bp)'].to_i

      ql = d['qlength'].to_i
      tl = d['tlength'].to_i

      qs = d['qstart'].to_i
      qe = d['qend'].to_i

      ts = d['tstart'].to_i
      te = d['tend'].to_i

      if ml < ql  # no full overlapp (A B C G H)
        if ts == 1 # A B C
          if te == tl 
            return 'C'
          else # A B
            if qe > m2
              return 'A'
            else
              return 'B'
            end
          end
        else # G H
          if qs < m2
            return 'G'
          elsif qs > m2
            return 'H'
          else
            return 'OOPS'
          end
        end
      else # full overlapp ( D E F )
        if qe > m2 # E F
          if qs > m2 
            return 'F'
          else
            return 'D'
          end
        else # D
          return 'D'
        end
      end
    end

    def identical_seqs?
      d['%similar'] == '100.00'
    end

    #
    # Scores
    #
    #   
    # TODO: consider
    # match_length (similarity) -> 'length of match (bp)'
    # match_score - 'evalue'?

    # TODO: bad classification is completely arbitrary, see constant
    def score_locus_overlap(bad = BAD_OVERLAPS)
      s = overlap_type
      case s 
      when bad.include?(s)
        0
      else
        1
      end 
    end

    # same == 1 (bad)
    # different = 0 
    def score_locus_difference
      query_locus == target_locus ? 1 : 0
    end

    def score_proportional_length(cutoff = 0.95)
      proportional_length > cutoff ? 1 : 0
    end

    def proportional_length
      d['length of match (bp)'].to_f / d['tlength'].to_f
    end

    def composite_score_exact_match_different_locus(bad = BAD_OVERLAPS)
      if score_locus_difference == 1 # different loci 
        return  score_locus_overlap(bad) == 1 ? 1 : 0 # with full overlap
      else
        1
      end
    end

  end
end
