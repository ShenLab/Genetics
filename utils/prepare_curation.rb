def main
  annovar = ARGV[0]

  File.new(annovar, 'r').each do |line|
    cols = line.split(/\t/)
    chr,pos, gene, ref, alt,  aa, freqesp, freq1kg, gerp, info, func = cols[0], cols[1], cols[6].gsub(/\"/, ""), cols[3], cols[4], cols[8], cols[10], cols[11], cols[37], cols[-1], cols[7]
    next if pos == "Start"
    if gene =~ /^(\S+)\((\S+)\)/
      geneName = $1
      aachange = "splicing;#{$2}"
      func = "splicing"
    else
      geneName = gene
      aachange = aa
    end
    
    if freqesp == ""
      freqesp = 0
    end
    
    if freq1kg == ""
      freq1kg = 0
    end
    
    freq = [freqesp.to_f, freq1kg.to_f].max
    
    if gerp == ""
      gerp = "NA"
    end

    infof = info.gsub(/\"/,"").split(';')
    sid, ad = infof[0], infof[1].split(',')[0,2]

    if func == ""
      func = "unknown"
    end
    func.gsub!(/\s+/, "_")
    puts "#{sid}\thg19\t#{geneName}\t#{chr}\t#{pos}\t#{ref}\t#{alt}\t0/1\t#{aachange}\tde_novo\tNA\t#{freq}\t#{gerp}\t#{func}\t#{ad.join("\t")}"
  end
end

main()
