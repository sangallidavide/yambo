TARGET_DIR=/home/sangalli/Data/Lavoro/Codici/yambo/the_wole_project/branches/gpl/devel
TRUNK_DIR=/home/sangalli/Data/Lavoro/Codici/yambo/the_wole_project/trunk/devel
rm -r $TARGET_DIR/*
cp -r $TRUNK_DIR/* $TARGET_DIR/
cd $TARGET_DIR
 ./sbin/yamboo_new.pl -p="KERR SURF YPP_SURF ELPH YPP_ELPH FFTW FFTSG OPENMP"
 chmod +x delete.batch 
./delete.batch
svn status | grep -v libxc
