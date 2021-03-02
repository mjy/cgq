require 'URI'
require 'net/http'
require 'JSON'


module Cgq
  class Records

    # Assumes you run `ruby src/go.rb` (TODO: better resolve absolute path)
    DERIVED_PATH = File.expand_path('../../../data/working/derived/', __dir__)

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
            @families[k] = { 
              family: {name: n['name'], id: n['id']}
            }
          else
            @families[k] = { 
              family: { name: nil, id: nil}
            }
          end

          if n1 = get_subfamily(v)  
            @families[k][:subfamily] = { 
              name: n1['name'],
              id: n1['id']
            }
          else
            @families[k][:subfamily] = { 
              name: nil, id: nil
            }
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
      return i['AE plate #'] if i
      nil
    end

    def plate_names
      plate_rows.collect{|k,v| v['AE plate #']}.uniq.compact
    end

    # @return Array
    #   [nil, nil] if not determined
    #   [column, row] if determined
    #
    # @param kind
    #   one of 'I# query' or 'I# trget'
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
      return nil if score_plate_similarity(row) == 0  # They are not on the same plate
      q = plate_xy(row, 'I# query')
      t = plate_xy(row, 'I# target')

      Math.sqrt((q[0] - t[0])**2 + (q[0] +t[0])**2).round(3)
    end

    def plate_data(i_number)
      plate_rows(i_number)
    end

    # !! Can not be move to Row because they reference plate data

    # @return [Float, nil]
    #    query - target
    def concentration_difference(row)
      c1 = query_qubit(row)
      c2 = target_qubit(row)

      return nil if (c1 =~ /fail/i) || (c2 =~ /fail/i) || c1.nil? || c2.nil?

      (c1.to_f - c2.to_f).round(3)
    end

    # @return [Float, nil]
    #    query - target
    def concentration_ratio(row)
      c1 = query_qubit(row)
      c2 = target_qubit(row)

      return nil if (c1 =~ /fail/i) || (c2 =~ /fail/i)  || c1.nil? || c2.nil?

      a, b = [c1.to_f, c2.to_f].sort # numerator is always on top
      
      # print " #{a} - #{b} "
      (a / b).round(3)
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
    # Exclude
    #
    #
    def exclude_score_new(row, focus: :target, qubit_ratio_cutoff: 0.33, qubit_consideration_cutoff: 3.0)

      return 0 if contaminated_score(row) != 3 # yes this should be cached

      qq = query_qubit(row)
      tq = target_qubit(row)

      return 1 if (qq =~ /fail/i) || (tq =~ /fail/i)  || qq.nil? || tq.nil?

      qq = qq.to_f
      tq = tq.to_f

      return 1 if qq <= qubit_consideration_cutoff && (tq <= qubit_consideration_cutoff)

      qr = concentration_ratio(row)

      if qr <= qubit_ratio_cutoff
        if focus == :query
            qq < tq  ? 1 : 0
        elsif focus == :target
          tq < qq  ? 1 : 0
        end
      else
        return 1 # TODO: uncertain
      end 
    end


    # Old methods

    def exclude_score(row, focus: :target, type: :ratio, concentration_cutoff: nil, composite_cutoff: [3,4,5] )
      send("exclude_score_#{type}", row, focus: focus, type: type, concentration_cutoff: concentration_cutoff, composite_cutoff: composite_cutoff)
    end

    def exclude_score_difference(row, focus: :target, type: :ratio, concentration_cutoff: 3.0, composite_cutoff: [3,4,5] )
      concentration_value = concentration_difference(row)

      return nil if concentration_value.nil?

      cs = composite_score(row, concentration_cutoff, :difference)

      v = nil

      if composite_cutoff.include?(cs)
        if focus == :query
          v = (concentration_value < 0)
        elsif focus == :target
          v = (concentration_value > 0)
        end

        if v && (concentration_value.abs > concentration_cutoff)
          return 1
        end
      end
      return nil
    end

    # @return [1, nil]
    #   whether to exclude this row based on the ratio of query/target qbit
    def exclude_score_ratio(row, focus: :target, type: :ratio, concentration_cutoff: 0.3, composite_cutoff: [3,4,5])
      cr = concentration_ratio(row)
      return nil if cr.nil?

      cs = composite_score(row, concentration_cutoff, :ratio)

      v = nil
      
      if composite_cutoff.include?(cs)
         d = concentration_difference(row) 

        if focus == :query
          v = (d < 0)
        elsif focus == :target
          v = (d > 0)
        end

        if v && (cr < concentration_cutoff)
          return 1
        end
      end
      return nil
    end

    #
    # Scores
    #

    # TODO
    def score_query_matching_multiple_target(row)
      # query_loci = more than 2 target loci
    end

    # @return [0, 1]
    #   If both subfamily and family match return 0, else 1.
    def score_taxon_difference(row)

      a = row.genus_pair.first
      b = row.genus_pair.last

      c = families[a]
      d = families[b]

      puts "#{a} is not found" if c.nil?
      puts "#{b} is not found" if d.nil?

      return 0 if c.nil? or d.nil?

      if (c['family']['id'] == d['family']['id']) && ( c['subfamily']['id'] == d['subfamily']['id'])
        0
      else
        1
      end
    end

    # @return [1,0]
    #   returns 1 if the plate is *the same*
    def score_plate_similarity(row)
      p1 = plate_rows[ row.d['I# query'] ]['AE plate #']
      p2 = plate_rows[ row.d['I# target'] ]['AE plate #']
      if p1 && p2
        p1 == p2 ? 1 : 0
      else
        0 
      end
    end

    def score_column_identity(row)
      if score_plate_similarity(row) == 1
        a = plate_xy(row, 'I# query')
        b = plate_xy(row, 'I# target')
        if (a[0] && b[0]) && (a[0] == b[0]) # Comparison is column to column
          1 
        else
          0
        end
      else
        0
      end
    end

    def score_concentration_difference(row, cutoff = 3)
      # they are on different plates, shouldn't be a contamination
      return 0 if score_plate_similarity(row) == 1

      if s = concentration_difference(row)
        cutoff = 3 if cutoff.nil?

        return s.abs > cutoff ? 1 : 0
      else
        0 # Hmm, likely not here
      end
    end

    def score_concentration_ratio(row, cutoff = 0.3)
      # they are on different plates, shouldn't be a contamination
      return 0 if score_plate_similarity(row) == 1

      if s = concentration_ratio(row)
        cutoff = 0.3 if cutoff.nil? # handle nil coming in

        return s.abs < cutoff ? 1 : 0
      else
        0
      end
    end

    # def absolute difference vs. min/max
    #   2x the volume as a cutoff
    # end

    # 100% -> significant (then possible below) -> match to same species, or match to contaminants -> one of the biggest flags !!
    def score_proportional_difference(row)
      if row.d['%similarity']  == '100.00' 
        1
      else
        0
      end
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

    # @param concentration_method [Symbol]
    #   one of :difference or :ratio
    def composite_score(row, concentration_cutoff = nil, concentration_method = :ratio)
      q = 0
      if concentration_method == :difference
        q = score_concentration_difference(row, concentration_cutoff)
      elsif concentration_method == :ratio
        q = score_concentration_ratio(row, concentration_cutoff)
      else
        raise 'bad concentration_method'
      end

      row.score_locus_overlap +
        row.score_locus_difference +
        row.score_proportional_length +
        score_taxon_difference(row) +
        score_plate_similarity(row) +
        score_column_identity(row) +
        q +
        score_proportional_difference(row)
    end

    def contaminated_score(row)
      row.score_proportional_length +
        score_taxon_difference(row) +
        score_plate_similarity(row)
        # score_proportional_difference(row) (all records)
    end

    def composite_score_identical_and_different_family(row)a
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

    def get_family(taxon_name_id)
      b = UCD_API + 'project_token=' + PROJECT_TOKEN + "&ancestors=true&nomenclature_group=FamilyGroup::Family&taxon_name_id[]=#{taxon_name_id}"
      u = URI(b)
      n = JSON.parse(::Net::HTTP.get(u))
      n[0]
    end

    def get_subfamily(taxon_name_id)
      b = UCD_API + 'project_token=' + PROJECT_TOKEN + "&ancestors=true&nomenclature_group=FamilyGroup::Subfamily&taxon_name_id[]=#{taxon_name_id}"
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
        viz.push( { group: x, variable: y,  value: score_difference(r) } ) # defaults to ratio
      end
      viz.uniq!
      viz
    end

    # @plate [String]
    #   the name of the plate
    def heatmap_count_records(plate = '0.5')
      t = 'I# query'

      v = {}
      # rows?
      all_rows.each do |r|
        next unless plate_name(r, t) == plate
        if a = plate_xy(r,t)
          if v[a]
            v[a][:total] += 1
          else
            v[a] = { total: 1, query_genus: r.d['query_genus'], i_num: r.d['I# query'], score: composite_score(r) }
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
