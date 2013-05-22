#!/usr/bin/env ruby

if ARGV.length != 2
  $stderr.puts "Usage: extract_from_multihetsep.rb <positions> <multihetsep_file> > <outfile>"
  exit
end

pos_str = ARGV[0]
input_file_name = ARGV[1]
allele_pos = pos_str.split(',').map {|s| s.to_i}

$stderr.puts "extracting positions #{allele_pos} from file #{input_file_name}"

File.open(input_file_name, "r") do |file|
  file.each_line do |line|
    fields = line.split
    allele_strings = fields[3].split(",")
    new_allele_strings = allele_strings.map do |allele_str|
      allele_pos.map {|i| allele_str[i]}.join
    end.join(",")
    puts "#{fields[0]} #{fields[1]} #{fields[2]} #{new_allele_strings}"
  end
end
