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
				> extracted_metadata.txt

#We need to take the extracted information and make a matrix with it:

awk -F '|' ' /^!Sample/ {
   # For lines starting with !Sample*, split them into fields using "|"
   for (i = 1; i <= NF; i++) { data[NR, i] = $i # Store each field in the "data" array (row,column)
   }
    max_fields = NF # Store the maximum number of fields encountered
}
END {
    # Now print the data transposed into columns
   for (i = 1; i <= max_fields; i++) 
	{ for (j = 1; j <= NR; j++) {
           # Print the field in the transposed format (columns)
		printf "%s", data[j, i]
		if (j < NR) { 
			printf "\t" # Tab after the field
		} #else {
   			#printf "%s", data[j, i]
		#}
	}
	}
} 
' extracted_metadata.txt > output.txt


