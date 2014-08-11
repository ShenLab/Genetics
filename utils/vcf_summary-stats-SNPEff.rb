## calcualte transition/transversion ratio from VCF file
## per sample
## also count the number of synonymous, missense, silent var per sample

## 
###  139 SNPEFF_EFFECT=CODON_CHANGE_PLUS_CODON_DELETION
#     26 SNPEFF_EFFECT=CODON_CHANGE_PLUS_CODON_INSERTION
#    271 SNPEFF_EFFECT=CODON_DELETION
#    101 SNPEFF_EFFECT=CODON_INSERTION
#   4451 SNPEFF_EFFECT=DOWNSTREAM
#   1770 SNPEFF_EFFECT=EXON
#    702 SNPEFF_EFFECT=FRAME_SHIFT
#     53 SNPEFF_EFFECT=INTERGENIC
#   2383 SNPEFF_EFFECT=INTRAGENIC
#   3366 SNPEFF_EFFECT=INTRON
#  54382 SNPEFF_EFFECT=NON_SYNONYMOUS_CODING
#      6 SNPEFF_EFFECT=NON_SYNONYMOUS_START
#    101 SNPEFF_EFFECT=SPLICE_SITE_ACCEPTOR
#    121 SNPEFF_EFFECT=SPLICE_SITE_DONOR
#    248 SNPEFF_EFFECT=START_GAINED
#     97 SNPEFF_EFFECT=START_LOST
#    665 SNPEFF_EFFECT=STOP_GAINED
#     54 SNPEFF_EFFECT=STOP_LOST
#  49581 SNPEFF_EFFECT=SYNONYMOUS_CODING
#     50 SNPEFF_EFFECT=SYNONYMOUS_STOP
#    253 SNPEFF_EFFECT=UPSTREAM
#    130 SNPEFF_EFFECT=UTR_3_PRIME
#    196 SNPEFF_EFFECT=UTR_5_PRIME




def main
  vcf = ARGV[0]
  codingOnly = ARGV[1]

  codingSwitch = 1
  if codingOnly != nil and codingOnly.to_i <= 0 ## 
    codingSwitch = -1
  end

  novelti = {}
  noveltv = {}
  knownti = {}
  knowntv = {}
  synon = {}
  missense = {}
  nonsense = {}
  splice = {}
  funcunknown ={}
  homo = {}
  het  = {}
  sid = [] # sample array
#  indel =0 

  lastpos = -1 
#  while line=ARGF.gets do 
  File.new(vcf, 'r').each do |line|
    next if line.match(/^\##/)
    cols=line.chomp.split(/\s+/)
    if line.match(/^#CHROM/)  # header
      cols[9..-1].each do |cc|
        cc = cc.sub(" ","_")
        sid << cc
        novelti[cc] = 0
        noveltv[cc] = 0
        knownti[cc] = 0 
        knowntv[cc] = 0
        synon[cc] = 0
        missense[cc] = 0
        nonsense[cc] = 0
        splice[cc] = 0
	funcunknown[cc] = 0
        homo[cc] = 0
        het[cc] = 0
      end
    else
      pos,name,ref,alt, passflag, info, gt = cols[1].to_i,cols[2],cols[3],cols[4].split(",")[0], cols[6], cols[7], cols[9..-1]
#	puts "#{pos}	#{name}	#{ref}	#{alt}	#{passflag}	#{info}	#{gt}"
      next if passflag != "PASS"
      if ref.size != alt.size ## indel   || pos == lastpos ## indels or same var
#	indel +=1
	next
      end
      ## SYN AND NONSYN

      funcunknownFlag, known, synonFlag, missenseFlag, nonsenseFlag, knowntiFlag, knowntvFlag ,noveltiFlag, noveltvFlag, spliceFlag = 0,0,0, 0, 0, 0,0, 0, 0,0
      coding = -1
      info.split(';').each do |l|  
	if l == "DB"		##GATK reports DB in the Info field if variant in the DBSNP version given to it in the --dbsnp option during unifiedgenotyper, we used dbsnp.132_excluding after_129 for better novel ti/tv etc. Also recommended by GATK. 
		known = 1
	end
        k,v = l.split('=')
     
        if  k =~ /SNPEFF_EFFECT/ 
          # functionalClass=nonsynonymousSNV
          #functionalClass=stopgainSNV
          #functionalClass=stoplossSNV
	#functionalClass=synonymousSNV
          #functionalClass=unknown
          coding = -1
          if  v == "SYNONYMOUS_CODING" or v == "SYNONYMOUS_STOP"  ## syn (Annovar reports "synonymousSNV" or "stopgainSNV" etc..)
            coding = 1
            synonFlag = 1
          elsif  v == "NON_SYNONYMOUS_CODING" or v == "NON_SYNONYMOUS_START"   # missense
            coding = 1
            missenseFlag = 1
          elsif  v == "STOP_GAINED"  or v == "STOP_LOST" 
            coding = 1
            nonsenseFlag = 1
          elsif  v == "SPLICE_SITE_DONOR" or v == "SPLICE_SITE_ACCEPTOR" 
            coding = 1
            spliceFlag = 1
          elsif v == "EXON"
            coding = 1
            funcunknownFlag = 1
          end
        end
      end
      
      next unless coding * codingSwitch > 0  # only take coding var
      
      next if pos == lastpos       

      ti = judge(ref,alt)
      if known == 1   #      if name =~ /^rs/  # known
        if ti == 1
          knowntiFlag = 1
        else
          knowntvFlag = 1
        end
      else  # novel
        if ti == 1
          noveltiFlag = 1
        else
          noveltvFlag = 1
        end
      end

      
      
#      $stderr.puts gt
      i = 0 
      gt.each do |genotypeinfo|
        subject = sid[i]
        genotype = genotypeinfo.split(":")[0]
        if genotype == '0/1'  ## het
          het[subject] += 1
          
        elsif genotype == '1/1' ## homo
          homo[subject] += 1
        end
        
        if genotype == '0/1' or genotype == '1/1'
          novelti[subject] += noveltiFlag
          noveltv[subject] += noveltvFlag
          knownti[subject] += knowntiFlag
          knowntv[subject] += knowntvFlag
          synon[subject] += synonFlag
          missense[subject] += missenseFlag
          nonsense[subject] += nonsenseFlag
          splice[subject] += spliceFlag
	  funcunknown[subject] += funcunknownFlag  
        end
        i += 1
      end
      lastpos = pos
    end
#    puts "#knownti/tv: #{knownti}/#{knowntv} = #{knownti/knowntv.to_f}"
#    puts "#novelti/tv: #{novelti}/#{noveltv} = #{novelti/noveltv.to_f}"
  end

  puts "\Samples\t#{sid.join("\t")}"
  print "all_SNV"
  sid.each {|s| print "\t#{knownti[s] + knowntv[s] + novelti[s] + noveltv[s]}"}
  print "\n"
  
  print "all_known"
  sid.each {|s| print "\t#{knownti[s] + knowntv[s]}"}
  print "\n"


  print "all_novel"
  sid.each {|s| print "\t#{novelti[s] + noveltv[s]}"}
  print "\n"

  print "known:ti/tv"
  sid.each {|s| print "\t#{knownti[s]}/#{knowntv[s]}"}
  print "\n"

  print "known:ti/tv-ratio"
  sid.each {|s| print "\t#{(knownti[s]/knowntv[s].to_f).round(3)}"}
  print "\n"


  print "novel:ti/tv"
  sid.each {|s| print "\t#{novelti[s]}/#{noveltv[s]}"}
  print "\n"

  print "novel:ti/tv-ratio"
  sid.each {|s| print "\t#{(novelti[s]/noveltv[s].to_f).round(3)}"}
  print "\n"

  print "silent"
  sid.each {|s| print "\t#{synon[s]}" }
  print "\n"
  
  print "missense"
  sid.each {|s| print "\t#{missense[s]}" }
  print "\n"

  print "nonsense"
  sid.each {|s| print "\t#{nonsense[s]}" }
  print "\n"

  print "splice"
  sid.each {|s| print "\t#{splice[s]}" }
  print "\n"

  print "function unknown"
  sid.each {|s| print "\t#{funcunknown[s]}" }
  print "\n"

  print "homozygous"
  sid.each {|s| print "\t#{homo[s]}" }
  print "\n"

  print "heterozygous"
  sid.each {|s| print "\t#{het[s]}" }
  print "\n"

#  print "all_Indels"
#  sid.each {|s| print "\t#{indel}" }
#  print "\n"

end


def judge(a1, a2)
  ti = 0
  if a1 == 'A'
    if a2 == 'G'  # ti
      ti = 1
    end
  elsif a1 == 'C'
    if a2 == 'T'
      ti = 1
    end
  elsif a1 == 'T'
    if a2 == 'C'
      ti = 1
    end
  elsif a1 == 'G'
    if a2 == 'A'
      ti = 1
    end
  end
  return ti
  
end

main()

