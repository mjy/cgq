require 'URI'
require 'net/http'
require 'JSON'



module Cgq
  class Records

    DERIVED_PATH = File.expand_path('../../data/working/derived/', __dir__)

    # All rows in possible_contaminants
    attr_accessor :all_rows

    # Rows excluding inverse queries (none?)
    attr_accessor :rows

    # @return [Hash]
    # { 
    # "I9733" => {
    #                      "AE I#" => "I9733",
    #                 "AE plate #" => nil,
    #   "DNA concentration: qubit" => nil,
    #                      "genus" => "Timioderus",
    #                    "species" => "peridentatus",
    #                   "position" => 138
    # }, ...}
    # 
    attr_accessor :plate_rows

    attr_accessor :families
    attr_accessor :genera
    attr_accessor :genus_ids
    attr_accessor :loci

    def initialize
      @all_rows = []
      @rows = []
      @plate_rows = {}
      @families = {}
      @genera = []
      @genus_ids = []
      @loci = []
    end

    #
    # Accessors
    #

    def rows
      return @rows if !@rows.empty?
      presence = {}
      all_rows.each do |r|
        k = [r.d['I# query'], r.d['I# target']]
        i = k.reverse
        if presence[k] || presence[i]
          next
        else
          presence[k] = nil
          @rows.push r
        end
      end
      @rows
    end

    # @return [Hash]
    #   genus_id: { name: 'name', id: 'id' }
    def families
      return @families if !@families.empty?

      p = DERIVED_PATH + '/family_metadata.json'
      if File.exists?(p)
        @families = JSON.parse(File.read(p))
      else
        @families = {}
        genus_ids.each do |k, v|
          if n = get_family(v)
            @families[k] = { name: n['name'], id: n['id'] }
          else
            @families[k] = { name: nil, id: nil}
          end
        end
      end
      @families
    end

    # @return [Hash]
    #   { 'name' => 'id', ... }
    def genus_ids
      return @genus_ids if !@genus_ids.empty?
      @genus_ids = {}
      p = DERIVED_PATH + '/ucd_genus_ids.csv'
      if File.exists?(p)
        CSV.read(p, col_sep: ",").each do |r|
          @genus_ids[r.first] = r.last
        end
      else
        # Should write first
        genera.each do |g|
          @genus_ids[g] = get_genus_id(g).to_s
        end
      end

      @genus_ids
    end

    def genera
      return @genera if !@genera.empty?
      g = {}
      unique_genus_pairs.each do |p|
        g[p.first] = nil
        g[p.last] = nil
      end
      @genera = g.keys.sort
    end

    def loci
      return @loci if !@loci.empty?
      g = {}
      all_rows.each do |r|

        [r.query_locus, r.target_locus].each do |l|
          if g[l]
            g[l] += 1
          else
            g[l] = 1
          end
        end
      end
      @loci = g
    end

    #
    # Plate metadata
    #

    # @param kind ['I# query', 'I# target']
    def plate_name(row, kind = 'I# query')
      i = plate_rows[ row.d[kind] ]
      i['AE plate #']
    end

    def plate_names
      plate_rows.collect{|k,v| v['AE plate #']}.uniq.compact
    end

    # @param kind
    #   one of 'I# query' or 'I# target'
    def plate_xy(row, kind = 'I# query')
      position = plate_rows[
        row.d[kind]
      ]
      if p = position['position']
        p.divmod(PLATE_HEIGHT)
      else
        return [nil, nil]
      end
    end

    def plate_cell_distance(row)
      return nil if score_plate_difference(row) == 0  # They are not on the same plate
      q = plate_xy(row, 'I# query')
      t = plate_xy(row, 'I# target')

      Math.sqrt((q[0] - t[0])**2 + (q[0] +t[0])**2).round(3)
    end

    def plate_data(i_number)
      plate_rows(i_number)
    end

    # @return [Float, -1]
    #    query - target
    def concentration_difference(row)
      c1 = query_qubit(row)
      c2 = target_qubit(row) 

      return nil if (c1 =~ /fail/) || (c2 =~ /fail/)
      c1.to_f - c2.to_f
    end

    def query_qubit(row)
      p1 = row.d['I# query']
      plate_rows[p1]['DNA concentration: qubit']
    end

    def target_qubit(row)
      p2 = row.d['I# target']
      plate_rows[p2]['DNA concentration: qubit']
    end

    #
    # Scores
    #

    # TODO
    def score_query_matching_multiple_target(row)
      # query_loci = more than 2 target loci
    end

    def score_taxon_difference(row)
      families[ row.genus_pair.first ]['id'] == families[ row.genus_pair.last ]['id'] ? 0 : 1
    end

    # i.e. same plate or different
    def score_plate_difference(row)
      p1 = plate_rows[ row.d['I# query'] ]['AE plate #']
      p2 = plate_rows[ row.d['I# target'] ]['AE plate #']
      p1 == p2 ? 1 : 0
    end

    # Scores range from -14.3 to 14.3
    # TODO: what is a reasonable cuttoff
    def score_concentration_difference(row, cutoff = 3) # consider using ratio rather than absolute
      # they are on different plates, shouldn't be a contamination
      return 0 if score_plate_difference(row) == 1
      if s = concentration_difference(row)
        return s.abs > cutoff ? 1 : 0
      else
        0 # Hmm, likely not here
      end
    end

    # def absolute difference vs. min/max
    #   2x the volume as a cutoff 
    # end

    def score_proportional_difference(row)
      row.d['%similar']  == '100.00' ? 1 :0
      # 100% -> significant (then possible below) -> match to same species, or match to contaminants -> one of the biggest flags !!
    end

    def composite_score_proportional_overlap(cutoff_length = 100) # basepairs
      if score_proportional_difference(row) == 1
        if row.d['length of match (bp)'] > cutoff_length
          1
        else
          0
        end
      end
    end

    def composite_score_difference(row)
      row.score_locus_difference + 
        row.score_proportional_length +  
        score_taxon_difference(row) + 
        score_plate_difference(row) + 
        score_concentration_difference(row) + 
        score_proportional_difference(row)
    end

    def composite_score_identical_and_different_family(row)
      score_taxon_difference(row) == 1 && row.identical_seqs?
    end

    #
    # Helpers, summaries, etc.
    #

    # @return [Array]
    #   A list of all unique pairs of genera
    def unique_genus_pairs
      taxa_pairs = {}
      all_rows.each do |r|
        v = r.genus_pair.sort
        taxa_pairs[v] = nil
      end

      taxa_pairs.keys.sort
    end

    # @return [Array]
    #   A list of all unique locus_pairs 
    def unique_locus_pairs
      locus_pairs = {}
      all_rows.each do |r|
        v = r.locus_pair.sort
        locus_pairs[v] = nil
      end
      locus_pairs.keys.sort
    end


    def offenders
      offenders = {}
      unique_genus_pairs.each do |p|
        p.each do |o|
          if offenders[o]
            offenders[o] += 1
          else
            offenders[o] = 1
          end
        end
      end

      offenders.sort{|a,b| b[1] <=> a[1]}.collect{|r| puts r[0].to_s + ': ' + r[1].to_s}
    end

    # If pteromalidae get subfamily
    def get_family(taxon_name_id)
      b = UCD_API + 'project_token=' + PROJECT_TOKEN + "&ancestors=true&nomenclature_group=FamilyGroup::Family&taxon_name_id[]=#{taxon_name_id}"
      u = URI(b)
      n = JSON.parse(::Net::HTTP.get(u))
      n[0]
    end

    def get_genus_id(genus)
      b = UCD_API + 'project_token=' + PROJECT_TOKEN + "&nomenclature_group=Genus&exact=true&name=#{genus}"
      u = URI(b)
      r = JSON.parse(::Net::HTTP.get(u))

      r.each do |n|
        # puts n['name_string']
        if (!n['cached_author_year'].nil?) &&  (n['cached_author_year'] =~ /\d\d\d\d/)
          return n['id']
        end
      end
      return nil
    end

    # @plate [String]
    #   the name of the plate
    def heatmap_score_differences(plate = '0.5')
      t = 'I# query'
      viz = {}
      rows.each do |r|
        next unless plate_name(r, t) == plate
        x, y = plate_xy(r, t).collect{|v| v.to_s.rjust(2, '0')}
        viz.push( { group: x, variable: y,  value: score_difference(r) } )
      end
      viz.uniq!
      viz
    end

    # @plate [String]
    #   the name of the plate
    def heatmap_count_records(plate = '0.5')
      t = 'I# query'

      v = {}
      # all_rows ?!
      rows.each do |r|
        next unless plate_name(r, t) == plate
        if a = plate_xy(r,t)
          if v[a]
            v[a][:total] += 1
          else
            v[a] = { total: 1, query_genus: r.d['query_genus'], i_num: r.d['I# query'], score: composite_score_difference(r) }
          end
        end
      end

      viz = []

      v.keys.each do |k|
        x, y = k.collect{|j| j.to_s.rjust(2, '0')}

        viz.push(
          { group: x,
            variable: y,
            value: v[k][:total],
            query_genus: v[k][:query_genus],
            i_num: v[k][:i_num],
            score: v[k][:score] })
      end
      viz
    end

    def runs_per_genus
      m = {}
      all_rows.each do |r|
        genera.each do |g|
          if r.d['query_genus'] == g
            if m[g]
              m[g].push r.d['I# query']
            else
              m[g] = [ r.d['I# query'] ]
            end
          end

          if r.d['target_genus'] == g
            if m[g]
              m[g].push r.d['I# target']
            else
              m[g] = [ r.d['I# target'] ]
              endputs data.loci.sort.collect{|l| l.join(',')}.join("\n")
            end
          else
            next
          end
        end
      end
      m.keys.each do |k|
        m[k].uniq!
        m[k].sort!
      end
      m
    end

    def loci_per_run
      m = {}
      all_rows.each do |r|
        genera.each do |g|
          if r.d['query_genus'] == g
            if m[g]
              m[g].push r.d['I# query']
            else
              m[g] = [ r.d['I# query'] ]
            end
          end

          if r.d['target_genus'] == g
            if m[g]
              m[g].push r.d['I# target']
            else
              m[g] = [ r.d['I# target'] ]
            end
          else
            next
          end

        end
      end
      m.keys.each do |k|
        m[k].uniq!
        m[k].sort!
      end
      m
    end

    def overlap_by_loci_by_genera
      d = {}
      rows.each do |r|
        g = r.genus_pair.sort
        l = r.locus_pair.sort
        o = r.overlap_type

        if d[g]
          if d[g][l]
            d[g][l][o] = nil
          else
            d[g][l] = {o => nil}
          end
        else
          d[g] = { l => {o => nil} }
        end
      end
      d
    end 

    def locus_overlap_by_i_num
      d = {}
      rows.each do |r|
        lq, lt = r.locus_pair
        o = r.overlap_type
        i = r.d['I# query']

        if d[lq]
          if d[lq][i]
            d[lq][i][lt] = nil
          else
            d[lq][i] = {lt => nil}
          end
        else
          d[lq] = { i => {lt => nil} }
        end
      end
      d
    end 

    def overlap_type_per_locus_pair
      d = {}
      rows.each do |r|
        a = r.locus_pair.sort
        b = r.overlap_type

        if d[a]
          d[a][b] = nil
        else
          d[a] = {b => nil}
        end
      end
      d
    end 

  end
end
