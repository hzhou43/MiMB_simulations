## Compile for linux binary
    #assuming the fftw3 library installed
    cd src
    make linux

## Compile for webassembly
### Compile fftw with emcc, assuming emsdk installed
    #setup emsdk
    source emsdkpath/emsdk_env.sh
    #download fftw to the parent directory of src and compile it with emcc
    cd ..
    wget http://www.fftw.org/fftw-3.3.10.tar.gz
    tar -zxf fftw-3.3.10.tar.gz
    cd fftw-3.3.10
    emconfigure ./configure
    emmake make
    ## comiple webassembly
    cd ../src
    make web

## Examples for running linux binary are in L8 directory to get excess chemical potential for L = 8
### Run brute-force insertion
    cd L8/bruteforce
    #generate script for wrapping run.sh with multiple CPUs
    bash genjobs.sh >jobs.sh
    #the real work
    bash jobs.sh
    #collect resuslt. The res.txt under each directory for temperature includes two columns. The first is the number density, the second is excess chemical potential in kT.
    bash col.sh
### Run FMAP insertion
    cd L8/fmap
    #generate script for wrapping run.sh with multiple CPUs
    bash genjobs.sh >jobs.sh
    #the real work
    bash jobs.sh
    #collect resuslt for res.txt
    bash col.sh
