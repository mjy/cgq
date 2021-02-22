require 'fileutils'
require 'erb'

module Cgq
  module Report

    CSV_EXPORT_PATH = File.expand_path('../../../out/csv/', __dir__).freeze
    HTML_EXPORT_PATH = File.expand_path('../../../out/html/', __dir__).freeze

    TEMPLATE_PATH = File.expand_path('viz_templates/', __dir__)

    class << self

      def write_scores(data, options = {})

        concentration_method = options[:concentration_method]
        concentration_method ||= :ratio

        concentration_cutoff = options[:concentration_cutoff]

        composite_cutoff ||= options[:composite_cutoff]
        composite_cutoff ||= [3,4,5]

        FileUtils.mkdir_p(CSV_EXPORT_PATH)
        CSV.open(CSV_EXPORT_PATH + "/scores.csv", "w") do |csv|
          csv << %w{
          exclude
          exclude_query
          query_genus
          query_family
          exclude_target
          target_genus
          target_family
          i_query
          i_target
          composite_score
          query_qubit
          target_qubit
          concentration_difference
          concentration_ratio
          plate_query
          plate_target
          query_plate_x
          query_plate_y
          target_plate_x
          target_plate_y
          qt_plate_distance
          overlap_type
          s_locus_overlap
          s_proportional_length
          s_taxon_difference
          s_plate_difference
          x_concentration_difference 
          s_concentration_ratio
          s_proportional_difference
          cs_exact_match_different_locus
          query_locus
          target_locus
          percent_similarity
          length_of_match
          query_length
          target_length
          unmatched_basepairs
          gap_open
          qstart
          qend
          tstart
          tend
          evalue
          ucd_genus_query_id
          ucd_genus_target_id
        }

          data.rows.each do |r|
            gq = r.d['query_genus']
            gt = r.d['target_genus']
            gq_id = data.genus_ids[gq] 
            gt_id = data.genus_ids[gt]
            qx,qy = data.plate_xy(r) 
            tx,ty = data.plate_xy(r,'I# target') 

            fam_q = data.families[gq] ? data.families[gq]['name'] : 'UNKNOWN'
            fam_t = data.families[gt] ? data.families[gt]['name'] : 'UNKNOWN'

            exclude_query = data.exclude_score(r, focus: :query, type: concentration_method,  concentration_cutoff: concentration_cutoff, composite_cutoff: composite_cutoff)
            exclude_target = data.exclude_score(r, focus: :target, type: concentration_method, concentration_cutoff: concentration_cutoff, composite_cutoff: composite_cutoff)

            exclude = (exclude_query || 0) + (exclude_target || 0) # presently never 2

            csv << [
              exclude,
              exclude_query,
              gq,
              fam_q, 
              exclude_target,
              gt,
              fam_t, 
              r.d['I# query'],
              r.d['I# target'],
              data.composite_score(r, concentration_cutoff, concentration_method),
              data.query_qubit(r),
              data.target_qubit(r),
              data.concentration_difference(r),
              data.concentration_ratio(r),
              data.plate_name(r),
              data.plate_name(r, 'I# target'),
              qx, 
              qy,
              tx, 
              ty,
              data.plate_cell_distance(r),
              r.overlap_type,
              r.score_locus_difference,
              r.score_proportional_length, 
              data.score_taxon_difference(r),
              data.score_plate_difference(r),
              data.score_concentration_difference(r),
              data.score_concentration_ratio(r),
              data.score_proportional_difference(r),
              r.composite_score_exact_match_different_locus,
              r.query_locus,
              r.target_locus,
              r.d['%similar'],
              r.d['length of match (bp)'],
              r.d['qlength'],
              r.d['tlength'],
              r.d['mismatch'],
              r.d['gapopen'],
              r.d['qstart'],
              r.d['qend'],
              r.d['tstart'],
              r.d['tend'],
              r.d['evalue'],
              gq_id, 
              gt_id 
            ] 
          end
        end
      end

      def write_genus_ids(data)
        FileUtils.mkdir_p(CSV_EXPORT_PATH)
        CSV.open(CSV_EXPORT_PATH + "/ucd_genus_ids.csv", "w") do |csv|
          data.genus_ids.each do |r|
            csv << r
          end
        end
      end

      def write_family_metadata(data)
        FileUtils.mkdir_p(CSV_EXPORT_PATH)
        File.open(CSV_EXPORT_PATH + "/family_metadata.json", "w") do |f|
          f.write(data.families.to_json)
        end
      end

      def write_overlap_loci_by_genera(data)
        FileUtils.mkdir_p(CSV_EXPORT_PATH)

        d = data.overlap_by_loci_by_genera

        CSV.open(CSV_EXPORT_PATH + "/overlap_by_loci_by_genus.csv", "w") do |csv|

          csv << %w{
            genus1
            genus2
            locus1
            locus2
            overlap_count
            overlap 
          }
          
          d.keys.each do |k|
            d[k].keys.each do |l|
              csv << k + l + [ d[k][l].keys.count ] + [d[k][l].keys.sort.join(' ')] 
            end
          end
        end
      end

      def locus_overlap_by_i_num(data)
        FileUtils.mkdir_p(CSV_EXPORT_PATH)

        d = data.locus_overlap_by_i_num

        CSV.open(CSV_EXPORT_PATH + "/locus_overlap_by_i_num.csv", "w") do |csv|

          csv << %w{
            query_locus 
            query_i_num 
            target_loci 
            count_target_loci 
          }

          d.each do |k, v|
            v.keys.each do |i|
              csv << [k] + [i] + [ d[k][i].keys.sort.join(',') ] + [ d[k][i].count ] 
            end
          end
        end
      end

      def overlap_type_per_locus_pair(data)
        FileUtils.mkdir_p(CSV_EXPORT_PATH)

        d = data.overlap_type_per_locus_pair

        CSV.open(CSV_EXPORT_PATH + "/overlap_type_per_locus_pair.csv", "w") do |csv|

          csv << %w{
            query_locus
            target_locus 
            overlap_types 
            overlap_type_count 
          }

          d.keys.sort.each do |k|
            csv << k + [ d[k].keys.sort.join(',') ] + [ d[k].keys.count ] 
          end
        end
      end

      # @param concentration_cutoff_range [Int]
      #   min 0, max 5
      def count_exclusion(data, composite_score_cutoff_range = [3,4,5], concentration_cutoff_range = [3] )
        CSV.open(CSV_EXPORT_PATH + "/count_contamination_per.csv", "w", col_sep: ',') do |csv|

          csv << %w{
            concentration_cutoff 
          } + composite_score_cutoff_range.collect{|c| "score_#{c}"}

          concentration_cutoff_range.each do |j|
            values = []
            composite_score_cutoff_range.each do |i|

              t = 0

              data.rows.each do |r|
                a = data.exclude_score(r, focus: :query, type: :difference, composite_cutoff: [i], concentration_cutoff: j)
                b = data.exclude_score(r, focus: :target, type: :difference, composite_cutoff: [i], concentration_cutoff: j) # TODO: make option
                t = t + a if !a.nil?
                t = t + b if !b.nil?
              end 
              values.push t

            end
            csv << [j] + values
          end
        end
      end

      # @param concentration_cutoff_range [Int]
      #   min 0, max 5
      def count_exclusion_ratio(data, composite_score_cutoff_range = [3,4,5], concentration_cutoff_range = [0.3] )
        CSV.open(CSV_EXPORT_PATH + "/count_contamination_ratio_per.csv", "w", col_sep: ',') do |csv|
          csv << %w{
            concentration_cutoff 
          } + composite_score_cutoff_range.collect{|c| "score_#{c}"}

          concentration_cutoff_range.each do |j|
            values = []
            composite_score_cutoff_range.each do |i|
              t = 0

              data.rows.each do |r|
                a = data.exclude_score_ratio(r, type: :ratio, focus: :query, composite_cutoff: [i], concentration_cutoff: j)
                b = data.exclude_score_ratio(r, type: :ratio, focus: :target, composite_cutoff: [i], concentration_cutoff: j)
                t = t + a if !a.nil?
                t = t + b if !b.nil?
              end 
              values.push t

            end
            csv << [j] + values
          end
        end
      end 

      def count_heatmaps(data)
        p = HTML_EXPORT_PATH  + '/heatmaps/'  
        FileUtils.mkdir_p(p)

        t = ERB.new(
          File.read(TEMPLATE_PATH + '/heatmap.html.tt').to_s
        )

        data.plate_names.each do |plate_name|
          filename = p + "plate_#{plate_name}.html"
          d = data.heatmap_count_records(plate_name)

          f = File.open( filename, "w" ) 
          f.puts t.result(binding)          
          f.close
        end
      end

      def score_difference_heatmap(data, plate = '0.5')
        viz = data.heatmap_count_records(plate)
        puts 'group,variable,value'
        puts viz.sort{|a,b,c| a <=> b}.collect{|r| r.join(',')}.join("\n")
      end

      def count_records_heatmap(data, plate = '0.5')
        viz = data.heatmap_count_records(plate)
        puts 'group,variable,value'
        puts viz.sort{|a,b,c| a <=> b}.collect{|r| r.join(',')}.join("\n")
      end

      def score_concentration_heatmap(data, cutoff = 2)
        # map concentrations to plate
      end

      def foo(data)
        puts data.loci.sort.collect{|l| l.join(',')}.join("\n")
      end

      # https://stackoverflow.com/questions/10236049/including-one-erb-file-into-another 
      def render(path, binding)
        content = File.read( TEMPLATE_PATH + '/' + path )
        t = ERB.new(content, nil, nil, '_sub01')
        t.result(binding)
      end

    end
  end
end
