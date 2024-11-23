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
  settings = {}
  settings["--minSpan"] = 0.9
  optHash = getopt()
  vcf = optHash["--vcf"]
  
  settings["--output"] = vcf  
  settings["--header"] = vcf
  settings.keys.sort.each do |s|
    if optHash.key?(s) 
      if s == "--output" or s == "--header"
        settings[s] = optHash[s]
      else
        settings[s] = optHash[s].to_f
      end
    end
  end
  
  
  samples, markers = loadVCF(vcf)
  
  $stderr.puts "number of samples: " + samples.size.to_s
  $stderr.puts "number of SNP markers: " + markers.size.to_s
 ## find SNPs that span 90% of samples
  minSet = []
  notCovered = samples.keys
  $stderr.puts "not covered: " + notCovered.length.to_s
  snps = findMaxSpan(samples, markers, notCovered, minSet, settings["--minSpan"])
  $stderr.puts snps
  puts snps.length
end

def loadVCF(vcf)
  
  #  if settings.key?("--output") 
#    prefix = settings["--output"]
  #  end
  
  samples = {}
  markers = {}
  fields = {:gt => 0}
  
  n, k = 0, 0
  
#  $stderr.puts vcf
  
  if vcf.match(/.gz$/)
    vio = Zlib::GzipReader.new(File.open(vcf))
  else
    vio = File.new(vcf, 'r')
  end
  
  flag = 0
#   while !vio.eof? 
#  File.new(vcf, 'r').each do |line|
  #     line = vio.gets(sep="\n")
  while line = vio.gets(sep="\n")
#    $stderr.puts line

#
    if line.match("^#CHROM") 
#      $stderr.puts line
      cols = line.chomp.split(/\s+/)
      #     $stderr.puts cols[-5..-1].join(";")
      #     $stderr.puts cols.size
      k = 8
      cols[9..-1].each do |sID|
        k = k + 1
        samples[sID] = {:col => k}
      end
    #      flag = 1
    elsif line.match("^#")
      next
    else
      cols=line.chomp.split(/\s+/)
#      $stderr.puts cols.size
      chr, pos, ref, alt, qual, pass, info  = cols[0], cols[1], cols[3], cols[4],  cols[5].to_f, cols[6], cols[7].split(';')
      #10      61852   .       T       C       999     .       DP=76;AF1=0.5;AC1=2;DP4=9,17,11,15;MQ=40;FQ=999;PV4=0.78,0.17,0.42,0.043        GT:PL:DP:SP:GQ  0/1:215,0,221:28:4:99   0/1:219,0,208:24:0:99
      #10      68575   .       C       T       999     .       DP=72;AF1=1;AC1=4;DP4=4,6,24,25;MQ=29;FQ=-60.5;PV4=0.73,1,1,1   GT:PL:DP:SP:GQ  1/1:235,29,0:31:4:63    1/1:239,41,0:28:2:75
      
      genotypes = cols[9..-1]
      coor = chr + "-" + pos
     
      if ref.size != alt.size
        ## do not consider indel
        next
      end
     
     
      s = {:AF => 0, :nhet => 0,  :ref => ref, :alt => alt, :id => coor, :hetind => {}}
      
 #     $stderr.puts genotypes.size
      s[:nhet] = countHet(genotypes)

#       format = cols[8].split(':')
      
	  # puts s[:id] + "\t" + s[:nhet].to_s

     
      
      info.each do |item|
        if item =~ /^AF\=(\S+)/
           s[:AF] = $1.to_f 
        end
      end
      
      ## HWE approximate filter
      if s[:nhet] >  samples.keys.size * 0.5
      	next
	  end 
	  
      markers[coor] = s
      
      samples.each_key do |sID|
        gt = cols[samples[sID][:col]]
        samples[sID][coor] = gt
        if gt == "0|1" or gt == "1|0" or gt == "0/1" or gt == "1/0"
        	s[:hetind][sID] = 1
        end	
        #       $stderr.puts "#{sID}\t#{samples[sID][:call].join("\t")}"
      end
    end
  end
  
  vio.close
#  $stderr.puts "top SNP: " + markers[0][:nhet].to_s 
#  $stderr.puts "top 100 SNP: " + markers[100][:nhet].to_s 
  
#  $stderr.puts "top SNP: " + nm[0][:nhet].to_s 
#  $stderr.puts "top 100 SNP: " + nm[100][:nhet].to_s 
  return samples, markers
end

def countHet(genotypes) 

#
  nhet = 0

  genotypes.each do |gt|
    if gt == "0|1" or gt == "1|0" or gt == "0/1" or gt == "1/0"
      nhet += 1
    end
  end
# $stderr.puts nhet
  return nhet
end

def findMaxSpan(samples, markers, notCovered, minSet, minSpan)


#   $stderr.puts "remaining N of SNPs: " + markers.keys.size.to_s
  totalSize = samples.keys.size
  topSNP = markers.keys.sort_by {|pos| markers[pos][:nhet]}.reverse[0]
  $stderr.puts "top SNP: " + topSNP + ", covering " + markers[topSNP][:nhet].to_s + " samples"
  #	notCovered = []
  #	if SNPs[topSNP][:nhet] > 0
  
  newTBD = []
  covered = []
  ## notCovered is an array of keys for samples not covered by previous minSet SNPs
  notCovered.each do |sID| 
  #	$stderr.puts samples[sID][topSNP]
    if samples[sID][topSNP] != "0|1" and samples[sID][topSNP] != "1|0"
      newTBD << sID
    else	
      covered << sID
		## remove the sample from the remaining markers      
    end
    
  end
  $stderr.puts "samples to be covered:" + notCovered.size.to_s + "->" + newTBD.size.to_s
    
  if newTBD.size < notCovered.size 
    minSet << topSNP
  end
  
  remainRatio = newTBD.size.to_f / totalSize.to_f
#  $stderr.puts "Ramining ratio: " + remainRatio.to_s 

	


  if remainRatio < (1 - minSpan) or markers[topSNP][:nhet] < 2
    return minSet
  else 
  
    ### resort remaining markers based on nhet among remaining samples

    covered.each do |sID| 
    ##        	s[:hetind][sID] = 1
    	markers.each_key do |coor|
			if markers[coor][:hetind].key?(sID)
				markers[coor][:hetind].delete(sID)
			end
			markers[coor][:nhet] = markers[coor][:hetind].keys.size
		end
	end		
  
    findMaxSpan(samples, markers, newTBD, minSet, minSpan)
  end
  
end

def getopt
  
  opts = GetoptLong.new(
    ["--vcf", "-v", GetoptLong::REQUIRED_ARGUMENT],
    ["--minSpan", "-m", GetoptLong::OPTIONAL_ARGUMENT], 
    ["--output", "-o", GetoptLong::OPTIONAL_ARGUMENT],
    ["--header", "-H", GetoptLong::OPTIONAL_ARGUMENT],
    ["--help", "-h", GetoptLong::NO_ARGUMENT]
                        )
  optHash = {}
  opts.each do |opt, arg|
    optHash[opt] = arg
  end
  
  $stderr.puts optHash
  
  if optHash.key?("--help") or !optHash.key?("--vcf") 
    $stderr.puts "Usage: ruby __.rb -v VCF -p ped [-o prefix]"
    $stderr.puts " options: "
    $stderr.puts "       --vcf [-v] VCF         input VCF file"
    $stderr.puts "       --maxSpan [-m] FLOAT    [0.9] minimum coverage of fraction of samples"
    $stderr.puts "       --output [-o] prefix   [stdout] prefix of output file; print to stdout by default"
    $stderr.puts "       --header [-H] header   [] optional VCF CHROM header"
    $stderr.puts "       --help  [-h]           show help information"
    exit
  end
  return optHash
  
end


main()
