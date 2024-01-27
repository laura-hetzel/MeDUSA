#!/usr/env bash


#HMDB
egrep "<.?hmdb( xml.*)?>|<.?metabolite>|^\s\s<accession>.*|^\s\s<name>.*|<monisotopic_molecular_weight>.*" \
   /usr/hmdb_metabolites.xml > /home/rstudio/hmdb_simple.xml

#LIPIDS
sed -z   's,> <EXACT_MASS>\n,SUMR_EM:,g' /usr/lipids.sdf | \
   sed -z 's,> <LM_ID>\nLMFA,SUMR_ID:LMFA,g' | \
   sed -z 's,> <NAME>\n,SUMR_NAME:,g' \
   > lipids_simple.tmp


egrep "SUMR_EM.*|SUMR_ID.*|SUMR_NAME.*" lipids_simple.tmp > lipids_simple1.tmp
sed -iz 's,SUMR_ID:,</mass>\n</lipid>\n<lipid>\n<id>,g' lipids_simple1.tmp
sed -iz 's,SUMR_NAME:,</id>\n<name>,g' lipids_simple1.tmp
sed -iz 's,SUMR_EM:,</name>\n<mass>,g' lipids_simple1.tmp
head -n 2 lipids_simple1.tmp >> lipids_simple1.tmp
tail -n +3 lipids_simple1.tmp > lipids_simple2.tmp && mv lipids_simple.tmp2 /home/rstudio/lipids_simple.xml

rm *.tmp