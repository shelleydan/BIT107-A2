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


####################### Downloading Matrix #######################

#wget -O ${WORKDIR}/${GEOID}.matrix.txt.gz \
#	https://ftp.ncbi.nlm.nih.gov/geo/series/${GEOLINK}/${GEOID}/matrix/${GEOID}_series_matrix.txt.gz

#Assigned a jobid to the wget command to let the download finish before moving on
#pid=$!
#wait $pid

#The downloaded file is zipped and is useless without unzipping first.
#gunzip ${GEOID}.matrix.txt.gz 

#The above is all commented out to continue with the next part REMOVE THIS

####################### Extracting Info #######################

### Key Info to Extract ###
# - !Sample_description - sample name (e.g. SKA for this assignment) 
# - !Sample_source_name_ch1 - The Cell Line Used
# - !Sample_characteristics_ch1 - infected vs mock
# - !Sample_characteristics_ch1 - sample time

cat ${GEOID}.matrix.txt | grep \
				-e '!Sample_description'\
				-e '!Sample_source_name_ch1'\
				-e '!Sample_characteristics_ch1'\
				-e '!Sample_characteristics_ch1'\
				>> extracted_metadata.txt


