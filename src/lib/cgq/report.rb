require 'fileutils'
require 'erb'

module Cgq
  module Report

    CSV_EXPORT_PATH = File.expand_path('../../../out/csv/', __dir__).freeze
    HTML_EXPORT_PATH = File.expand_path('../../../out/html/', __dir__).freeze

    TEMPLATE_PATH = File.expand_path('viz_templates/', __FILE__)

    class << self

      def write_scores(data, options = {})
        FileUtils.mkdir_p(CSV_EXPORT_PATH)
        CSV.open(CSV_EXPORT_PATH + "/scores.csv", "w") do |csv|
          csv << %w{
          query_genus
          target_genus
          ucd_genus_query_id
          ucd_genus_target_id
          query_family
          target_family
          i_query
          i_target
          plate
          plate_x
          plate_y
          qt_plate_distance
          query_qubit
          target_qubit
          concentration_difference
          overlap_type
          s_locus_difference
          s_proportional_length
          s_taxon_difference
          s_plate_difference
          s_concentration_difference 
          s_proportional_difference
          cs_sum_difference
          cs_exact_match_different_locus
          cs_predicted_contaminant
          }
          data.rows.each do |r|
            gq = r.d['query_genus']
            gt = r.d['target_genus']
            gq_id = data.genus_ids[gq] 
            gt_id = data.genus_ids[gt]
            x,y = data.plate_xy(r) 

            csv << [
              gq,
              gt,
              gq_id, 
              gt_id, 
              data.families[gq]['name'],
              data.families[gt]['name'],
              r.d['I# query'],
              r.d['I# target'],
              data.plate_name(r),
              x, 
              y,
              data.plate_cell_distance(r),
              data.query_qubit(r),
              data.target_qubit(r),
              data.concentration_difference(r),
              r.overlap_type,
              r.score_locus_difference,
              r.score_proportional_length, 
              data.score_taxon_difference(r),
              data.score_plate_difference(r),
              data.score_concentration_difference(r),
              data.score_proportional_difference(r),
              data.composite_score_difference(r),
              r.composite_score_exact_match_different_locus,
              'todo_predicted_contaminant', # TODO: query or target
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
              csv << k + l + [ d[k][l].keys.count ] + [d[k][l].keys.join(' ')] 
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
              csv << [k] + [i] + [ d[k][i].keys.join(',') ] + [ d[k][i].count ] 
            end
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
