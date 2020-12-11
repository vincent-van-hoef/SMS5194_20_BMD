jb build ../report/
ls
cd ..
ls
cd report/
ghp-import -npf _build/html/
cd ..
jb build ../report/
jb build report/
ls
cd report/
ln -s ../scripts/analysis.Rmd 
cd ..
jb build report/
jb clean report/
jb build report/
ls
cd scripts/
cd ..
jb build report/
