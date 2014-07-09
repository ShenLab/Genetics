#!/usr/bin/env ruby

require 'getoptlong'

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
  settings["--minADIndel"] =  6  ## min # of reads carrying alternative allele for indel
  settings["--minDP"] = 12  # min depth in parents
  settings["--phenotype"] = ""
  settings["--minPL"] = 70
  settings["--minPLP"] = 30
  settings["--minPLIndel"] = 70
  settings["--maxAAF"] = 0.015
  settings["--maxFreq"] = 0.001
  settings["--maxAC"]  = 3

  optHash = getopt()
  vcf = optHash["--vcf"]
  
  settings["--output"] = vcf  
  settings.keys.sort.each do |s|
    if optHash.key?(s) 
      if s == "--phenotype"  or s == "--output"
        settings[s] = optHash[s]
      else
        settings[s] = optHash[s].to_f
      end
    end
  end
  
  samples=countSamples(settings["--phenotype"], vcf)
#  $stderr.puts samples
  
  filterVCF(vcf,settings,samples)  # gt: gene -> pos -> sample -> genotype, 

end

def filterVCF(vcf, settings, samples)

#  if settings.key?("--output") 
#    prefix = settings["--output"]
#  end

  o = File.new(settings["--output"] + ".mut.vcf", 'w')
#  e = File.new(prefix + ".non-mut.vcf", 'w')
  firstline = 1
  fields = {:gt => 0}


  probands = []
#   $stderr.puts samples
  samples.each_key do |s|
    if samples[s][:pheno] == "2" and samples[s][:parents].size == 2
      probands << s
    end
  end

  $stderr.puts "Number of complete trios:  #{probands.size}"
  $stderr.puts "All probands: #{probands.join("\t")}"
   
  File.new(vcf, 'r').each do |line|
    if line.match("^#")
#      if line.match("^#CHROM") 
#        o.puts "\#MutSample\t#{line}"
#      else
      o.puts line
#      end
#      e.puts line
    else
      cols=line.chomp.split(/\t/)
      ref, alt, qual, pass, info  = cols[3], cols[4],  cols[5].to_f, cols[6], cols[7].split(';')
      #10      61852   .       T       C       999     .       DP=76;AF1=0.5;AC1=2;DP4=9,17,11,15;MQ=40;FQ=999;PV4=0.78,0.17,0.42,0.043        GT:PL:DP:SP:GQ  0/1:215,0,221:28:4:99   0/1:219,0,208:24:0:99
      #10      68575   .       C       T       999     .       DP=72;AF1=1;AC1=4;DP4=4,6,24,25;MQ=29;FQ=-60.5;PV4=0.73,1,1,1   GT:PL:DP:SP:GQ  1/1:235,29,0:31:4:63    1/1:239,41,0:28:2:75

      
      format = cols[8].split(':')

      if info[0]== "INDEL" or ref.size != alt.size
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
      alleleCounts = 0

      info.each do |item|
        if item =~ /^1KG\.score\=(\S+)/
          freq = [$1.to_f, freq].max
        elsif item =~ /^ESP\d+\.score\=(\S+)/
          freq = [$1.to_f, freq].max
        elsif item =~ /^AC\=(\S+)/
          alleleCounts = $1.to_i
        end
      end
      
      if freq > settings["--maxFreq"]
#        e.puts line
        next
      end
      

      if alleleCounts > settings["--maxAC"]
        next 
      end

#      gt_tumor,gt_norm =  cols[9].split(':'), cols[10].split(':')
      
      samples.each_key do |sID|
        k = samples[sID][:col]
        samples[sID][:call] = cols[k].split(':')
 #       $stderr.puts "#{sID}\t#{samples[sID][:call].join("\t")}"
      end
      

#      exclusionarray = []
      
      mut = []
      ada = []
      probands.each do |sid|
        #        $stderr.puts sid
      #  $stderr.puts "#{sid}"
 #       $stderr.puts "#{samples[sid][:call].join("\t")}"
        gt = samples[sid][:call][fields[:gt]] 

        if samples[sid][:call][fields[:gt]] == "0/0" or samples[sid][:call][fields[:gt]] == './.'
          #          $stderr.puts "case not 0/1"
          next
        end
#         $stderr.puts gt
        
        
        if plf = samples[sid][:call][fields[:pl]] 
          pl = plf.split(',')
        else
          next
        end
        if adf = samples[sid][:call][fields[:ad]]
          ad = adf.split(',')

        else
          next
        end

        if (indelflag == 0 and (pl[0].to_i < settings["--minPL"] or ad[1].to_i < settings["--minAD"] )) or (indelflag == 1 and (pl[0].to_i < settings["--minPLIndel"] or ad[1].to_i < settings["--minADIndel"] ))
        #  $stderr.puts "case not mut"
          next
        end
        
        exclusionflag = 0                             
        samples[sid][:parents].each do |p|
          if samples[p][:call][fields[:gt]] != "0/0" 
         #   $stderr.puts "control not 0/0???"
            exclusionflag = 1
          else
            pl = samples[p][:call][fields[:pl]].split(',')
            ad = samples[p][:call][fields[:ad]].split(',')
            dp = ad[0].to_f + ad[1].to_f
            aaf = ad[1].to_f / dp

            if  (indelflag == 0 and pl[1].to_i < settings["--minPLP"]) or (indelflag == 1 and pl[1].to_i < settings["--minPLIndel"])  or aaf > settings["--maxAAF"]  or dp < settings["--minDP"]
              #             $stderr.puts "control not 0/0"
              exclusionflag = 1
            end
          end
        end
  #      $stderr.puts exclusionflag
        if exclusionflag == 0
          mut << sid
          ada << adf
        end
      end
        
      
      if mut.size > 0 
      #  o.puts "#{mut.join(',')}\t#{line}"
        str = cols[0..6].join("\t") + "\t" + mut.join(":") + ";" + ada.join(":")  + ";" + cols[7..-1].join("\t")
        o.puts str
    #  else
    #    e.puts line
      end
    end
  end
  o.close

end

def countSamples(phenotype, vcf)
##example .fam file:
# #family sample paternalID maternalID gender phenotype(1/2)
# 1018 BW2 0 0 0 1
# 1018 BW1 0 0 0 2
# 1018 WB13 0 0 0 2
  samples = {}
  n, k = 0, 0
  File.new(vcf, 'r').each do |line|
    n += 1
    if line.match("^#CHROM") 
      cols = line.chomp.split(/\s+/)
      k = 8
      cols[9..-1].each do |sID|
        k = k + 1
        samples[sID] = {:col => k, :fID => nil, :pheno => "0" }
      end
    elsif n > 2000
      break
    end
  end

  File.new(phenotype, 'r').each do |line|
    next if line.match(/^#/)
    
    cols = line.chomp.split(/\s+/)
#    $stderr.puts cols.join("\t")
    fID, sID, pID, mID, gender, pheno = cols[0], cols[1], cols[2], cols[3], cols[4], cols[5]
    if samples.key?(sID)
      samples[sID][:fID] = fID
      samples[sID][:pID] = pID
      samples[sID][:mID] = mID
      samples[sID][:gender] = gender
      samples[sID][:pheno] = pheno
    end
  end

#  $stderr.puts samples

  samples.each_key do |sID|
    if samples[sID][:pheno] == "2"  ## case / proband
      samples[sID][:parents] = []
      samples.each_key do |s|
        if samples[s][:fID] == samples[sID][:fID] and samples[s][:pheno] == "1" 
          samples[sID][:parents] << s
        end
      end
    end
  end
  return samples

end

def getopt
  
  opts = GetoptLong.new(
                        ["--vcf", "-v", GetoptLong::REQUIRED_ARGUMENT],
                        ["--minAD", "-m", GetoptLong::OPTIONAL_ARGUMENT],
                        ["--minADIndel", "-M", GetoptLong::OPTIONAL_ARGUMENT],
                        ["--minDP", "-d", GetoptLong::OPTIONAL_ARGUMENT],
                        ["--minPLIndel", "-T", GetoptLong::OPTIONAL_ARGUMENT],
                        ["--phenotype", "-p", GetoptLong::REQUIRED_ARGUMENT],
                        ["--maxFreq", "-f", GetoptLong::OPTIONAL_ARGUMENT],
                        ["--minPL", "-t", GetoptLong::OPTIONAL_ARGUMENT], 
                        ["--minPLP", "-l", GetoptLong::OPTIONAL_ARGUMENT],
                        ["--maxAAFSigma", "-a", GetoptLong::OPTIONAL_ARGUMENT], 
                        ["--output", "-o", GetoptLong::OPTIONAL_ARGUMENT],
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
    $stderr.puts "       --phenotype [-p]       file specify pedigree or phenotypes [PLINK .fam format]"
    $stderr.puts "       --minAD [-m] INT       [5] min number of alt-allele reads for SNVs"
    $stderr.puts "       --minADIndel [-M] INT  [6] min number of alt-allele reads for INDELs"
    $stderr.puts "       --minDP [-d] INT       [10] min depth in parents/control"
    $stderr.puts "       --minPL [-t] INT       [60] min PL diff for SNVs in proband/case"
    $stderr.puts "       --minPLP [-l] INT      [30] min PL diff for SNVs in parents/control"
    $stderr.puts "       --minPLIndel [-T] INT  [60] min PL diff for INDELs"
    $stderr.puts "       --maxFreq [-f] FLOAT   [0.001] max allele frequency in 1kG or GO-ESP"
    $stderr.puts "       --maxAAF [-a] FLOAT    [0.03] max ratio of alternative allele reads in control/parents " 
    $stderr.puts "       --output [-o] prefix   [stdout] prefix of output file; print to stdout by default"
    
    $stderr.puts "       --help  [-h]           show help information"
    exit
  end
  return optHash
  
end


main()
