## pick variants by chr and position

def main 
  vcf = ARGV[0]
  list = ARGV[1]

  variants = getList(list)
  
  selectV(vcf, variants)


end

def selectV(vcf, v)
  File.new(vcf, "r").each do |line|
    if line.match("^#")
      puts line
    else
      if line =~ /^(\S+)\s+(\d+)/
        chr, pos = $1, $2
        if v.key?(chr) and v[chr].key?(pos)
          puts line
        end
      end
    end
  end
end

def getList(l)
  v = {}
  File.new(l, "r").each do |line|
    if line=~ /^(\S+)\s+(\d+)/
      chr, pos = $1, $2
      
      if !v.key?(chr)
        v[chr] = {}
      end

      v[chr][pos] = 1
      
    end
  end
  return v

end

main()
