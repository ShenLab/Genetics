#!/usr/bin/env ruby

require 'getoptlong'
require 'zlib'

### 

## filter vcf to call de novo / germline mutations
# rules:
# 1. genotype call: 0/0 in parents, 0/1 or 1/1 in offspring
# 2. PL cutoffs
# 3. min AD in offspring
# 4. min DP in parents
# 5. max alternative allele fraction in parents
# 6. max allele frequency in populations: reference panel AND within-cohort
# 7. 

def main 

  optHash = getopt()
  input = optHash["--input"]
  metrics = optHash["--metrics"]

  
  
  genes = loadMetrics(metrics)
  
  $stderr.puts genes.keys.length

  index = 1
  
  File.new(input, 'r').each do |line|
    line.chomp!
    cols = line.split(/\s+/)
    if line.match("^Chr")   # header
        i = 1
        cols[1..-1].each do |col|
            if col.match?("Gene")
                index = i
                $stderr.puts index
            end
            i = i + 1
        end
        puts line + "\tpLI"
  
        
    else 
   #     $stderr.puts index
        gene = cols[index]
#        $stderr.puts gene
        pLI = "-1"
        if genes.key?(gene)
            pLI = genes[gene]
            $stderr.puts gene + "\t" + pLI
        end
        puts line + "\t" + pLI
    end
 
  end
end


def loadMetrics(m)
  genes = {}
  index = 1
  File.new(m, 'r').each do |line|
    cols = line.split(/\s+/)
    if line.match("^gene")   # header
        i = 1
        cols[1..-1].each do |col|
            if col.match?("pLI")
                index = i
                $stderr.puts index
            end
            i = i + 1
        end
       
    else 
        gene, pLI = cols[0], cols[index]
        genes[gene] = pLI
    end
 
  end
  return genes
  


end


def getopt
  
  opts = GetoptLong.new(
    ["--input", "-i", GetoptLong::REQUIRED_ARGUMENT],
    ["--metrics", "-m", GetoptLong::REQUIRED_ARGUMENT],
    ["--column", "-c", GetoptLong::OPTIONAL_ARGUMENT],
    ["--output", "-o", GetoptLong::OPTIONAL_ARGUMENT],
    ["--help", "-h", GetoptLong::NO_ARGUMENT]
                        )
  optHash = {}
  opts.each do |opt, arg|
    optHash[opt] = arg
  end
  
  
  if optHash.key?("--help") or !optHash.key?("--metrics") 
    $stderr.puts "Usage: ruby __.rb -i input -m gene-metrics [-o prefix]"
    $stderr.puts " options: "
    $stderr.puts "       --input [-i]          input tab/space delimited file"
    $stderr.puts "       --metrics [-m]      gene metrics file"
    $stderr.puts "       --column [-c]       the index of column which specify gene symbol"
   
    $stderr.puts "       --output [-o] prefix   [stdout] prefix of output file; print to stdout by default"
    $stderr.puts "       --help  [-h]           show help information"
    exit
  end
  return optHash
  
end


main()
