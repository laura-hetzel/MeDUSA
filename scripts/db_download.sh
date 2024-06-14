#!/usr/env bash
curl https://hmdb.ca/system/downloads/current/hmdb_metabolites.zip -o hmdb_metabolites.zip
echo "HMDB Downloaded"
unzip hmdb_metabolites.zip
mv hmdb_metabolites.xml /usr/hmdb_metabolites.xml

curl "https://lipidmaps.org/files/?file=LMSD&ext=sdf.zip" -o lipids.zip
echo "LIPIDS Downloaded"
unzip lipids.zip
mv structures.sdf /usr/lipids.sdf
