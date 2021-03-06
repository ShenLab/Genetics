## calcualte transition/transversion ratio from VCF file
## per sample
## also count the number of synonymous, missense, silent var per sample



def main
  vcf = ARGV[0]
  minGQ = ARGV[1]

  if minGQ == nil 
    minGQ = 40
  else
    minGQ = minGQ.to_i
  end

  sid = [] # sample array
#  indel =0 

  lastpos = -1 

#  while line=ARGF.gets do 
  File.new(vcf, 'r').each do |line|
    next if line.match(/^\##/)
   #  $stderr.puts line
    cols=line.chomp.split(/\s+/)
    if line.match(/^#CHROM/)  # header
      sampleCols = []
      cols[9..-1].each do |cc|
        cc = cc.sub(" ","_")
        sid << cc
        sampleCols << "REF_#{cc}\tALT_#{cc}"
      end
      outputline = "\#type\tCHR\tPOS\t#{sampleCols.join("\t")}"
      puts outputline
    else
      
      chr,pos,name,ref,alt, passflag, info, fields,  gt = cols[0], cols[1].to_i,cols[2],cols[3],cols[4].split(",")[0], cols[6], cols[7], cols[8], cols[9..-1]
#	puts "#{pos}	#{name}	#{ref}	#{alt}	#{passflag}	#{info}	#{gt}"
      next if passflag != "PASS"
      outputline = "SNV\t"
      if ref.size != alt.size ## indel   || pos == lastpos ## indels or same var
#	indel +=1
        outputline = "INDEL\t"
#	next
      end


      outputline << "#{chr}\t#{pos}\t"
      

      fa = fields.split(":")
      i = 0
      indexGT, indexAD, indexDP, indexGQ = 0, 1, 2, 3
      fa.each do |f|
        if f == "GT"
          indexGT = i
        elsif f == "AD"
          indexAD = i
        elsif f == "DP"
          indexDP = i
        elsif f == "GQ"
          indexGQ = i
        end
        i += 1
      end
      
      gt.each do |genotypeinfo|
        genotypeFields = genotypeinfo.split(":")
        genotype, ad, gq = genotypeFields[indexGT], genotypeFields[indexAD].split(','), genotypeFields[indexGQ].to_i
        if genotype == '0/1' and gq >= minGQ  ## het
          outputline << "#{ad[0]}\t#{ad[1]}\t"
        else
          outputline << "NA\tNA\t"
        end
      end
      puts outputline
      lastpos = pos
    end
    #    puts "#knownti/tv: #{knownti}/#{knowntv} = #{knownti/knowntv.to_f}"
#    puts "#novelti/tv: #{novelti}/#{noveltv} = #{novelti/noveltv.to_f}"
  end

end


main()

