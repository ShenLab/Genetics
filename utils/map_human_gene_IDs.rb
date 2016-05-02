## ensembl ID map:
# A1BG    ENSG00000121410


def main
  input=ARGV[0]
  mapping = ARGV[1]

  genes = readMap(mapping)
  
  addGeneNameID(input, genes)


end

def addGeneNameID(input, genes)
  File.new(input, "r").each do |line|
    cols = line.chomp.split(/\s+/)
    gid = cols[0]
    if genes.key?(gid)
      gsymbol = genes[gid]
      puts  "#{gsymbol}\t#{cols[1]}"
    end
  end
end 


def readMap(mapping)
  genes = {}

  File.new(mapping, "r").each do |line|
    cols = line.chomp.split(/\t/)
    next unless cols.size == 2
    geneSymbol, geneID  = cols[0], cols[1]
    genes[geneID] = geneSymbol
  end

  return genes
end


main()
