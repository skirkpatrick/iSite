#!/usr/bin/env ruby

if ARGV.size < 2
    puts "USAGE: ./collect.rb <input file> <output file>"
    exit
end

header = "JOBID subfuncProb asymmetry selfloopLoss fusionProb fissionProb actualAsymmetry selfloops order size tris trips CC numComponents lgComponentOrder lgComponentSize lgComponentTris lgComponentTips lgComponentCC\n"

ofile = File.open(ARGV[-1], "w+")

ofile.print header

for i in (0..ARGV.length-2)
    ifile = File.open(ARGV[i], "r")
    ifile.readline
    ofile.puts ifile.readline
end

