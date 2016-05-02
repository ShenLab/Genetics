
def main
  genes = {}
  while line = ARGF.gets do 
    cols = line.chomp.split(/\s+/)
    gene, type = cols[0], cols[1]
    next if gene == "Gene.refgene"
      
    if !genes.key?(gene)
      genes[gene] = {}
      genes[gene][:LGD] = 0
      genes[gene][:DMIS] = 0
      genes[gene][:MIS] = 0
      genes[gene][:SYN] = 0
    end
    
    if type == "LOF" or type == "LGD"
      genes[gene][:LGD] += 1
    elsif type == "DMIS" or type == "D-mis" 
      genes[gene][:DMIS] += 1
      genes[gene][:MIS] += 1
    elsif type == "MIS"
      genes[gene][:MIS] += 1
    elsif type == "silent"
      genes[gene][:SYN] += 1
    end
    
  end
  puts "gene\tN_LGD\tN_Dmis\tN_Mis"
  genes.each_key do |gene|
    puts "#{gene}\t#{genes[gene][:LGD]}\t#{genes[gene][:DMIS]}\t#{genes[gene][:MIS]}"
  end
    
  
end
main()
