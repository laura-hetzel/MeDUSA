#!/usr/env bash
chmod 777 /usr/hmdb_metabolites.xml
chmod 777 /usr/lipids.sdf;

#HMDB
egrep "<.?hmdb( xml.*)?>|<.?metabolite>|^\s\s<accession>.*|^\s\s<name>.*|<monisotopic_molecular_weight>.*" \
   /usr/hmdb_metabolites.xml > /home/rstudio/hmdb_simple.xml

#LIPIDS
sed -z   's,> <EXACT_MASS>\n,SUMR_EM:,g' /usr/lipids.sdf | \
   sed -z 's,> <LM_ID>\n,SUMR_ID:,g' \
   > lipids_simple.tmp

#not all entries have NAME, some only have SYSTEMATIC_NAME
#sed -z 's,> <NAME>\n,SUMR_NAME:,g' | \
#   sed -z 's,> <NAME>\n,SUMR_NAME:,g' \
#   > lipids_simple.tmp

sed -i 's/,/-/g' lipids_simple.tmp
egrep "SUMR_EM.*|SUMR_ID.*|SUMR_NAME.*" lipids_simple.tmp > lipids_simple1.tmp
echo "IDENTITY, EXACT_MASS" > /home/rstudio/lipids_simple.csv
sed -z 's/SUMR_ID://g' lipids_simple1.tmp | \
sed -z 's/\nSUMR_EM:/, /g' >> /home/rstudio/lipids_simple.csv
  #sed -z 's/\nSUMR_NAME:/, /g' | \
rm /usr/hmdb_metabolites.xml
rm /usr/lipids.sdf
rm *.tmp
