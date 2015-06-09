
def main
  while line = ARGF.gets do 
    cols = line.chomp.split(/\t/)
    chr, pos, s, p = cols[0], cols[1], cols[2], cols[3]
    
    id = "NA"
    if s =~ /^\"(\S+)\(/
      id = $1
    end
    
    pa = p.split(/\,/)
    codonc = {}
    pa.each do |pp|
      x = pp.split(/\:/)[-1].sub("\"", "")
      codonc[x] = 1 
    end

    puts "#{chr}\t#{pos}\t#{id}\t#{codonc.keys.join(",")}"
  end
end

main()
