#!/usr/bin/env ruby

if ARGV.size < 2
    puts "USAGE: ./average_dist.rb <input files> <output file>"
    exit
end

ifile = Array.new

for i in (0..ARGV.length-2)
    ifile.insert(-1, File.open(ARGV[i], "r"))
end

ofile = File.open(ARGV[-1], "w+")

before = Hash.new(0)
after = Hash.new(0)

nvals = 0

ifile.each do |f|
    nvals += 1
    f.readline
    f.each do |line|
        break if line == "\n"
        before[line[/^[0-9]+:/][0..-2]] += line[/ [0-9]*$/][1..-1].to_i
    end
    f.readline
    f.each do |line|
        after[line[/^[0-9]+:/][0..-2]] += line[/ [0-9]*$/][1..-1].to_i
    end
end

ofile.puts "Before simplification"

before.each do |k,i|
    ofile.puts "#{k}: #{i/nvals}"
end

ofile.puts "\nAfter simplification"

after.each do |k,i|
    ofile.puts "#{k}: #{i/nvals}"
end
