## ensembl ID map:
# ENSMUSG00000064372      ENSMUST00000082423      mitochondrially encoded tRNA proline [Source:MGI Symbol;Acc:MGI:102478] mt-Tp

## Kallisto output:
# target_id       length  eff_length      est_counts      tpm

def main
  input=ARGV[0]
  mapping = ARGV[1]

  transcripts = readMap(mapping)
  
  addGeneNameID(input, transcripts)


end

def addGeneNameID(input, transcripts)
  File.new(input, "r").each do |line|
    if line =~ /^target/  # header
      puts "#{line.chomp}\tGeneID\t#GeneName"
    else
      cols = line.chomp.split(/\s+/)
      tid = cols[0]
      if transcripts.key?(tid)
        gid = transcripts[tid][:geneID]
        gname = transcripts[tid][:geneName]
      else
        gid = "NA"
        gname = "NA"
      end
      puts  "#{line.chomp}\t#{gid}\t#{gname}"
    end
  end
end 


def readMap(mapping)
  transcripts = {}

  File.new(mapping, "r").each do |line|
    cols = line.chomp.split(/\t/)
    geneID, transcriptID, geneDesc , geneName = cols[0], cols[1], cols[2], cols[3]
    transcripts[transcriptID] = {} unless transcripts.key?(transcriptID)
    
    transcripts[transcriptID][:geneID] = geneID
    transcripts[transcriptID][:geneName] = geneName
    transcripts[transcriptID][:geneDesc] = geneDesc

  end

  return transcripts
end


main()
