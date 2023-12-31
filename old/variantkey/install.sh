git clone https://github.com/tecnickcom/variantkey.git
cd variantkey/r
make build
cd -
R -e 'install.packages("variantkey/r/variantkey", repos=NULL, type="source")'
