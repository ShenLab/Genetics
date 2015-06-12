# problematic input format:

# #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO:Case;Accession;Gender;DOB;Autism;ID;Seizures;Genotype;Quality;ControlQ1;ControlQ2;TotalCov;XomeFreq
# 17      59556019        .       C       T       99      PASS    FAM140084;1402431;F;8/3/11;N;N;N;het;99;99;99;736;0
# 6       152420025       .       C       T       99      PASS    FAM150139;1456838;M;5/31/13;N;N;N;het;99;99;99;212;0
# 2       232458736       .       C       T       99      PASS    FAM131438;1331779;M;9/23/03;N;Y;N;het;99;99;99;266;0
# 1       214209146       .       C       T       99      PASS    FAM141226;1425284;F;10/15/10;N;N;N;het;99;99;99;369;0
# 20      472926  .       T       C       99      PASS    FAM141226;1425284;F;10/15/10;N;N;N;het;99;99;99;499;0
# 1       220088869       .       -       GCTA

while line = ARGF.gets do 
      cols = line.chomp.split(/\t/)
      next if line.match("^#CHROM") or cols.size < 4
      offset = cols[3].size  - 1
      posend = cols[1].to_i +	offset	

      puts "#{cols[0]}\t#{cols[1]}\t#{posend}\t#{cols[3]}\t#{cols[4]}"
end   
