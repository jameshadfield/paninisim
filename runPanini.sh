#!/bin/bash

#### given a folder (rtab/x) which contains a series of Rtab file, curl each one and save to json/x


echo "$1"

folname="genomes"


# curl -H "Content-type: text/plain; charset=UTF-8" -X POST --data-binary @./genomes/panGenomeSize.3.100.100.100.0.2.2.0.Rtab http://panini.wgsa.net/api/1.0/rtab


folderName=$1
if [[ -d json/$folderName ]]; then
	echo "json/$folderName already exists. Exiting!"
fi



for i in genomes/gamma/*Rtab; do j=json${i#genomes}; echo $i; echo $j; curl -H "Content-type: text/plain; charset=UTF-8" -X POST --data-binary @./${i} http://panini.wgsa.net/api/1.0/rtab > ./${j%Rtab}json; done