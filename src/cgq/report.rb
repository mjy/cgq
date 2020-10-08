require 'fileutils'

module Cgq
  module Report

    CSV_EXPORT_PATH = File.expand_path('../../out/csv/', __dir__)

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

      def score_difference_heatmap(data, plate = '0.5')
        viz = data.score_difference_heatmap(plate)
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
      
    end
  end
end
