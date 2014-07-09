## grep by a column, keep order of the query list, return NA for null 
require 'getoptlong'


def main
  col = 0
  optHash = getopt()
  if optHash.key?("--col")
    col = optHash["--col"].to_i - 1
  end

  list = getQuery(optHash["--query"])
  hash = {}
  list.each do |i|
    hash[i] = "NA"
  end

  doGrep(hash, optHash["--source"], col)

  list.each do |i|
    puts "#{i}\t#{hash[i]}"
  end
  

end

def doGrep(h, s, col)
#  puts col
  File.new(s,'r').each do |line|
    cols = line.lstrip.chomp.split(/\s+/)
#    puts cols.join("\t")
    k = cols[col]
 #   puts k
    if h.key?(cols[col])
      h[cols[col]] = line.lstrip.chomp
   #   puts line
    end
  end
  return h
end


def getQuery(f)
  array = []
  File.new(f,'r').each do |line|
    if line =~ /^(\S+)/
      array << $1
    end
  end
  return array

end


def getopt
  
  opts = GetoptLong.new(
                        ["--query", "-f", GetoptLong::REQUIRED_ARGUMENT],
                        ["--source", "-s", GetoptLong::REQUIRED_ARGUMENT],
                        ["--col", "-c", GetoptLong::OPTIONAL_ARGUMENT],
                        ["--help", "-h", GetoptLong::NO_ARGUMENT]
                        )
  optHash = {}
  opts.each do |opt, arg|
    optHash[opt] = arg
  end
  
 # $stderr.puts optHash
  
  if optHash.key?("--help") or !optHash.key?("--query") or !optHash.key?("--source")
    $stderr.puts "Usage: ruby __.rb -f query_file -s source_file [-c column_number]"
    exit
  end
  return optHash
  
end




main()
