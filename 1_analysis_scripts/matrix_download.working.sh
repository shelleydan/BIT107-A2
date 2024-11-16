#!/usr/bin/bash

#Assessment GEOID is GSE217504

#Test GEOID is GSE217511

#I think spam kills it so heres another GSE217510

#Setting Variables
WORKDIR=$(pwd)
GEOID=GSE217504
GEOLINK=${GEOID::-3}nnn

# sed -i 's/$GEOLINK/$GEOLINK\nnn/'

echo ${GEOID}
echo ${GEOLINK}


#######################Downloading Matrix#########################

wget -O ${WORKDIR}/${GEOID}.matrix.txt.gz \
	https://ftp.ncbi.nlm.nih.gov/geo/series/${GEOLINK}/${GEOID}/matrix/${GEOID}_series_matrix.txt.gz

pid=$!
wait $pid

gunzip ${GEOID}.matrix.txt.gz
