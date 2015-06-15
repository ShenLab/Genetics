## find shared genes that are mutated among different groups of cases

# calculate jaccard index!

def main
  phenofile = ARGV[0]
  mutfile = ARGV[1]

  groups = readPheno(phenofile)

 #  $stderr.puts groups.keys.size
  s = countSharing(mutfile, groups)
  ji = jaccardIndex(s, groups) 
#  s.each do |pair, n|
#    puts "#{pair}\t#{n}"
#  end

#  puts ji
  ji.each do |pair, v|
    $stderr.puts "#{pair}\t#{v}\t#{s[pair]}"
  end

  ## input for heatmap

  m = groups.values.max
  
  
  for i in 0..m 
  #  print "#{i}"
    for j in 0..m
      if j > i 
        pair = "#{i}-#{j}"
        print "\t#{ji[pair]}"
      elsif j == i
        print "\t1"
      else
        pair = "#{j}-#{i}"
        print "\t#{ji[pair]}"
      end
    end
    print "\n"
  end
  

end

def jaccardIndex(sharing, groups)
  ji = {}
  g = []
  for i in 0..groups.values.max
    g << i
  end

  g.each do |g1|
    g.each do |g2|
      next if g2 <= g1
      pair = "#{g1}-#{g2}"
   #   puts sharing[pair]
      ji[pair] = sharing[pair].to_f / (sharing["#{g1}-#{g1}"] + sharing["#{g2}-#{g2}"] - sharing[pair])
    end
  end
      
  return ji
end

def countSharing(mutfile, groups)
  mut = {}
  sharing = {}
  m = groups.values.max
  for i in 0..m
    for j in i..m
      sharing["#{i}-#{j}"] = 0
    end
  end

#  puts sharing
  File.new(mutfile, 'r').each do |line|
    cols = line.chomp.split(/\t/)
    next if line=~ /^Case/
    cid, gene, type  = cols[0], cols[8], cols[9]
#    puts type
    mut[gene] = {} unless mut.key?(gene)
    if type != "silent" and type != "nonframeshift" and type != "Bmis"   and type  != "." and groups.key?(cid)
      g = groups[cid]
      if mut[gene].key?(g)
        mut[gene][g] += 1
      else
        mut[gene][g] = 1
      end
#      puts mut[gene]
    
    end
  end
  

  mut.each_key do |gene|
    gid = mut[gene].keys.sort
    l = gid.length
 #   puts "#{gene}\t#{l}"
    for i in 0..(l-1)
      g1 = gid[i]
      for j in i..(l-1)
        g2 = gid[j]
#        puts "#{i}\t#{j}\t#{g1}-#{g2}"
        sharing["#{g1}-#{g2}"] += 1
      end
    end
  end
  return sharing
end

def readPheno(p)
  g = {}
  File.new(p, 'r').each do |line|
    cols = line.chomp.split(/\s+/)
    if line =~ /^Case/ # header
      $stderr.puts line
    else
      cid = cols[0]
      p = []
      g[cid] = 0
      j = 0
      cols[1..-1].reverse.each do |i|
        if i == "Y"
          g[cid] += 2**j
        end
        j += 1
      end
    end
  end
  return g
end


main()
