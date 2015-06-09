## cosmic format:
# Gene name	Accession Number	Gene CDS length	HGNC ID	Sample name	ID_sample	ID_tumour	Primary site Site subtype Primary histology Histology subtype Genome-wide screen Mutation ID Mutation CDS Mutation AA Mutation Description Mutation zygosity Mutation GRCh37 genome position Mutation GRCh37 strand SNP FATHMM prediction Mutation somatic status Pubmed_PMID ID_STUDY Sample sourceTumour origin Age Comments

# column 16: mutation type
# output: gene  CDS_length total nonsense missense silent frameshift inframe deletion other

def main
  
  cosmic = ARGV[0]

  genes = {}
  total = {}
  nonsense = {}
  missense = {}
  silent = {}
  frameshift = {}
  inframe = {}
  deletion = {}
  other = {}

  File.new(cosmic, "r").each do |line|
    cols = line.split(/\t/)
    
    name, mut = cols[0], cols[15]
    if !genes.key?(name)
      genes[name] = cols[2]
      total[name] = 0
      nonsense[name] = 0
      missense[name] = 0
      silent[name] = 0
      frameshift[name] = 0
      inframe[name] = 0
      deletion[name] = 0
      other[name] = 0
    end

    if mut.match("Missense")
      missense[name] += 1
    elsif mut.match("Nonsense")
      nonsense[name] += 1
    elsif mut.match("silent")
      silent[name] += 1
    elsif mut.match("Frameshift") or mut.match("frameshift")
      frameshift[name] += 1
    elsif mut.match("inframe")
      inframe[name] += 1
    elsif mut.match("Whole gene deletion")
      deletion[name] += 1
    else
      other[name] += 1
    end
    total[name] += 1

  end
  
  puts "#gene\tCDS_length\ttotal\tnonsense\tmissense\tsilent\tframeshift\tinframe\tgene_deletion\tother"

  total.sort_by {|a , b| b}.reverse.each do |gene, t|
    puts "#{gene}\t#{genes[gene]}\t#{t}\t#{nonsense[gene]}\t#{missense[gene]}\t#{silent[gene]}\t#{frameshift[gene]}\t#{inframe[gene]}\t#{deletion[gene]}\t#{other[gene]}"
  end

end


main()



