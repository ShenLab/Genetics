## input: a VCF file

## output: adding a field to info fields, listing the carriers (ID, GT, AD)

def main
  vcf = ARGV[0]
  
  idindex = []

  File.new(vcf, 'r').each do |line|
    cols = line.chomp.split(/\t/)
    if line.match(/^\#CHROM/) # header
      total = cols.size
      idindex = cols[9..-1]
    elsif (!cols[0].match(/^\#/)) 
      carriers = []
      gt = []
      ad = []
      adindex = 0
      i = 0
      cols[8].split(':').each do |j|
        if j == "AD"
          adindex = i
        end
        i += 1
      end
      i = 0
      cols[9..-1].each do |s|
        items = s.split(":")
        if items[0] == "0/1" or items[0] == "1/1"  # carrier
          carriers << idindex[i]
          gt << items[0]
          ad << items[adindex]
        end
        i += 1
      end
      cols[7] = carriers.join("|") + ";" + gt.join("|") + ";"  + ad.join("|") + ";" + cols[7]
    end

    puts cols.join("\t")
  end
end



main()



