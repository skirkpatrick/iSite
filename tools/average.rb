#!/usr/bin/env ruby

if ARGV.size < 2
    puts "USAGE: ./average.rb <input files> <output file>"
    exit
end

ifile = Array.new

for i in (0..ARGV.length-2)
    ifile.insert(-1, File.open(ARGV[i], "r"))
end

ofile = File.open(ARGV[-1], "w+")
sums = Array.new(19, 0)

nvals =0

header = "JOBID subfuncProb asymmetry selfloopLoss fusionProb fissionProb actualAsymmetry selfloops order size tris trips CC numComponents lgComponentOrder lgComponentSize lgComponentTris lgComponentTips lgComponentCC\n"

ifile.each do |f|
    f.readline
    f.each do |line|
        nvals += 1
        vals = line.split(' ')
        (0...sums.length).each do |i|
            sums[i] += vals[i].to_f
        end
    end
end

ofile.print header
(0...sums.length).each do |i|
    if [1, 2, 3, 4, 5, 6, 12, 18].include? i
        ofile.print "%.3f " % [sums[i]/nvals]
    else
        ofile.print "%d " % [sums[i]/nvals]
    end
end
