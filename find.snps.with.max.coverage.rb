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
  settings["--minAD"] = 6  ## min # of reads carrying alternative allele for SNV
  settings["--minADIndel"] =  8  ## min # of reads carrying alternative allele for indel
  settings["--minDP"] = 12  # min depth in parents
  settings["--minPL"] = 70
  settings["--minPLP"] = 30
  settings["--minPLIndel"] = 80
  settings["--maxAAF"] = 0.015
  settings["--maxFreq"] = 0.001
  settings["--maxAC"]  = 3

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
  

  samples, SNPs=loadVCF(vcf)

    #  $stderr.puts samples
 ## find SNPs that span 90% of samples
  minSet = []
  notCovered = samples.keys
  snps = findMaxSpan(samples, SNPs, notCovered, minSet, 0.9)

end

def loadVCF(vcf, settings)

#  if settings.key?("--output") 
#    prefix = settings["--output"]
#  end

	samples = {}
	SNPs = {}
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

    if line.match("^#CHROM") 
 #     $stderr.puts line
      cols = line.chomp.split(/\s+/)
 #     $stderr.puts cols[-5..-1].join(";")
 #     $stderr.puts cols.size
      k = 8
      cols[9..-1].each do |sID|
        k = k + 1
        samples[sID] = {:col => k }
      end
#      flag = 1
	elsif line.!match("^#")
        cols=line.chomp.split(/\s+/)
      $stderr.puts cols.size
      chr, pos, ref, alt, qual, pass, info  = cols[0], cols[1], cols[3], cols[4],  cols[5].to_f, cols[6], cols[7].split(';')
      #10      61852   .       T       C       999     .       DP=76;AF1=0.5;AC1=2;DP4=9,17,11,15;MQ=40;FQ=999;PV4=0.78,0.17,0.42,0.043        GT:PL:DP:SP:GQ  0/1:215,0,221:28:4:99   0/1:219,0,208:24:0:99
      #10      68575   .       C       T       999     .       DP=72;AF1=1;AC1=4;DP4=4,6,24,25;MQ=29;FQ=-60.5;PV4=0.73,1,1,1   GT:PL:DP:SP:GQ  1/1:235,29,0:31:4:63    1/1:239,41,0:28:2:75

	genotypes = cols[9..-1]
		coor = chr + "-" + pos
	 s = {:AF => 0, :nhet => 0,  :ref => ref, :alt => alt}
	 s[:nhet] = countHet(genotypes)
      SNPs[coor] = s
      format = cols[8].split(':')

      if ref.size != alt.size
        indelflag = 1
      else
        indelflag = 0
      end
      

      # if firstline == 1 
      i = 0
      format.each do |k|
        if k == "GT"
          fields[:gt] = i
        elsif k == "PL" 
          fields[:pl] = i
        elsif k == "AD"
          fields[:ad] = i
        end
        i+=1
      end
      #end

#      #$stderr.puts fields

     # firstline = 2
      
      freq = 0 
      hetcount = 0

      info.each do |item|
        if item =~ /^AF\=(\S+)/
          freq = [$1.to_f, freq].max
        end
      end
      
      if freq < settings["--minFreq"]
#        e.puts line
        next
      end
      


#      gt_tumor,gt_norm =  cols[9].split(':'), cols[10].split(':')
      
      samples.each_key do |sID|
        gt = cols[samples[sID][:col]]
        samples[sID][coor] = gt
       
 #       $stderr.puts "#{sID}\t#{samples[sID][:call].join("\t")}"
      end
      


    end
  end

  vio.close
	return samples, SNPs
end

def countHet(genotypes) 
	nhet = 0
	genotypes.each |gt|
	    if gt == "0|1" or gt == "1|0" or gt == "0/1" or gt == "1/0"
			nhet += 1
		end	
	end
	return nhet
end

def findMaxSpan(samples, SNPs, notCovered, minSet, minSpan)
	
	
	totalSize = samples.keys.size
	topSNP = SNPs.keys.sort {|s1, s2| SNPs[s2][:nhet] <=> SNPs[s1][:nhet]}[0]
#	notCovered = []
#	if SNPs[topSNP][:nhet] > 0
		
	newArray = []
		## notCovered is an array of keys for samples not covered by previous minSet SNPs
	notCovered.each do |sID| 
		if samples[sID][topSNP] != "0|1" and samples[sID][topSNP] != "1|0"
			newArray << sID
		end
	end
	SNPs.delete(topSNP)
	if newArray.size < notCovered.size 
		minSet << topSNP
	end
		
	if newArray.size / totalSize < minSpan  or SNPs.keys.size < 1
		return minSet
	else 
		findMaxSpan(samples, SNPs, newArray, minSet, minSpan)
	end

end

def getopt
  
  opts = GetoptLong.new(
                        ["--vcf", "-v", GetoptLong::REQUIRED_ARGUMENT],
                        ["--minAD", "-m", GetoptLong::OPTIONAL_ARGUMENT],
                        ["--minADIndel", "-M", GetoptLong::OPTIONAL_ARGUMENT],
                        ["--minDP", "-d", GetoptLong::OPTIONAL_ARGUMENT],
                        ["--minPLIndel", "-T", GetoptLong::OPTIONAL_ARGUMENT],
                        ["--minFreq", "-f", GetoptLong::OPTIONAL_ARGUMENT],
                        ["--minPL", "-t", GetoptLong::OPTIONAL_ARGUMENT], 
                        ["--minPLP", "-l", GetoptLong::OPTIONAL_ARGUMENT],
                        ["--maxAAFSigma", "-a", GetoptLong::OPTIONAL_ARGUMENT], 
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
    $stderr.puts "       --minAD [-m] INT       [5] min number of alt-allele reads for SNVs"
    $stderr.puts "       --minADIndel [-M] INT  [6] min number of alt-allele reads for INDELs"
    $stderr.puts "       --minDP [-d] INT       [10] min depth in parents/control"
    $stderr.puts "       --minPL [-t] INT       [60] min PL diff for SNVs in proband/case"
    $stderr.puts "       --minPLP [-l] INT      [30] min PL diff for SNVs in parents/control"
    $stderr.puts "       --minPLIndel [-T] INT  [60] min PL diff for INDELs"
    $stderr.puts "       --minFreq [-f] FLOAT   [0.01] min allele frequency"
    $stderr.puts "       --maxAAF [-a] FLOAT    [0.03] max ratio of alternative allele reads in control/parents " 
    $stderr.puts "       --output [-o] prefix   [stdout] prefix of output file; print to stdout by default"
    $stderr.puts "       --header [-H] header   [] optional VCF CHROM header"
    $stderr.puts "       --help  [-h]           show help information"
    exit
  end
  return optHash
  
end


main()
