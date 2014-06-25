def main
  annovar = ARGV[0]

  File.new(annovar, 'r').each do |line|
    cols = line.split(/\t/)
    chr,pos, gene, ref, alt,  aa, info, func = cols[0], cols[1], cols[6].gsub(/\"/, ""), cols[3], cols[4], cols[8], cols[-1], cols[7]
    next if pos == "Start"
    if gene =~ /^(\S+)\((\S+)\)/
      geneName = $1
      aachange = "splicing;#{$2}"
      func = "splicing"
    else
      geneName = gene
      aachange = aa
    end
    
    infof = info.gsub(/\"/,"").split(';')
    sid, ad = infof[0], infof[1].split(',')[0,2]

    if func == ""
      func = "unknown"
    end
    func.gsub!(/\s+/, "_")
    puts "#{sid}\thg19\t#{geneName}\t#{pos}\t#{ref}\t#{alt}\t0/1\t#{aachange}\tde_novo\tNA\t#{ad.join("\t")}\t#{func}"
  end
end

main()
