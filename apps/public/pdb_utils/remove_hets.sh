# unzip everything
find debugging -name '*.pdb.gz' | xargs gunzip
# remove all HETNAM lines
find debugging -name '*.pdb'  | xargs sed -i '' '/HETNAM/d'
# rezip everything
find debugging -name '*.pdb'  | xargs gzip