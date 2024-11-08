#!/usr/bin/bash

#Assessment GEOID is GSE217504

#Test GEOID is GSE217511

#I think spam kills it so heres another GSE217510

#Setting Variables
WORKDIR=$(pwd)

######################## USER INPUT GEO ID #######################
echo "Enter the geoprofile ID of the Project: (The GEO for this assignment is: GSE217504)" 
read GEOID #prompts the user for the ID

GEOLINK=${GEOID::-3}nnn #Customises the ID to work with the NCBI link

######################## FILE DOWNLOAD #######################

wget -O ${WORKDIR}/${GEOID}.matrix.txt.gz \
	https://ftp.ncbi.nlm.nih.gov/geo/series/${GEOLINK}/${GEOID}/matrix/${GEOID}_series_matrix.txt.gz

wget -O ${WORKDIR}/${GEOID}.counts.csv.gz \
	https://www.ncbi.nlm.nih.gov/geo/download/?acc=${GEOID}&format=file&file=${GEOID}%5Fhost%5Fcounts%5Fmatrix%2Ecsv%2Egz

#Assigned a jobid to the wget command to let the download finish before moving on
pid0=$!
wait $pid0
pid1=$!
wait $pid1

#The downloaded file is zipped and is useless without unzipping first.
gunzip ${GEOID}.matrix.txt.gz 
gunzip ${GEOID}.counts.csv.gz

#The above is all commented out to continue with the next part REMOVE THIS

####################### EXTRACTION #######################

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

# The data is extracted into extracted_metadata.txt

lineN=$(cat extracted_metadata.txt | wc -l)
echo ${lineN}

for i in $(seq 1 ${lineN}); do
	sed -n ${i}p extracted_metadata.txt > ${i}_matrix.txt
	sed -i 's/\t/\n/g' ${i}_matrix.txt
done

sed -i '1 s/^.*$/Cell Type/' 1_matrix.txt
sed -i '1 s/^.*$/Cell Type/' 2_matrix.txt
sed -i '1 s/^.*$/Cell Line/' 3_matrix.txt
sed -i '1 s/^.*$/Genotype/' 4_matrix.txt
sed -i '1 s/^.*$/Test Condition/' 5_matrix.txt
sed -i '1 s/^.*$/Sample Time/' 6_matrix.txt
sed -i '1 s/^.*$/ID/' 7_matrix.txt

paste *_matrix.txt > targets.txt
pid2=$!
wait $pid2

rm ${WORKDIR}/*_matrix.txt #Just some simple tidying

####################### CLEANING TARGETS #######################

sed -i 's/"//g' targets.txt
sed -i 's/cell line: //g' targets.txt
sed -i 's/cell type: //g' targets.txt
sed -i 's/genotype: //g' targets.txt
sed -i 's/treatment: //g' targets.txt
sed -i 's/time point_(in_hours): //g' targets.txt

echo "targets.txt is ready for use"
