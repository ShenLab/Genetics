f1 = ARGV[0]
f2 = ARGV[1]

list1 = {}
list2 = {}

File.new(f1, "r").each do |line|
  cols = line.chomp.split(/\s+/)
  list1[cols[0]] = line.chomp
end

File.new(f2, "r").each do |line|
  cols = line.chomp.split(/\s+/)
  if list1.key?(cols[0])
    o = list1[cols[0]]
  else
    o = "#{cols[0]}\tNA"
  end
  puts o

end
