#!/bin/sh

mkdir bin
mkdir bin/tmp
# rm -r source/tmp/*

printf "OSHUN version control: \n\n" > version.txt 
git log --pretty=format:'%h : %s' --graph >> version.txt
printf "\n\n (PICKSC made this)" >> version.txt

cd source
cp makefiles/makefile makefile
#Create temporary file with new line in place
cat main.cpp | sed -e "s/commit hash # /commit hash # `git rev-parse --short HEAD`/" > main2.cpp
#Copy the new file over the original file
mv main2.cpp main.cpp


make
#make debug
