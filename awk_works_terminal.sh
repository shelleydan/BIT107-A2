#!/bin/bash
awk -F '\t' ' BEGIN { col = 0 # Initialize column counter
}
# For lines starting with "!Sample", process the data
/^!Sample/ { col++ # Increment column counter for each new header headers[col] = $1 # Store 
    the header (first field) in the headers array
    # Store the values (all fields after the first) in the data array
    for (i = 2; i <= NF; i++) { data[col, i-1] = $i # Store data in the array
    }
}
END {
    # Print headers
    for (i = 1; i <= col; i++) { printf "%s", headers[i] if (i < col) { printf "\t" # Print 
            tab after the header
        }
    }
    print "" # Newline after the headers
    # Print data (transposed)
    max_rows = 0
    # Find the maximum number of rows
    for (i = 1; i <= col; i++) { for (j = 1; j <= NF-1; j++) { if (data[i, j] != "") { 
                max_rows = (j > max_rows) ? j : max_rows
            }
        }
    }
    # Print each row in transposed format
    for (i = 1; i <= max_rows; i++) { for (j = 1; j <= col; j++) { printf "%s", data[j, i] if 
            (j < col) {
                printf "\t" # Add tab between columns
            }
        }
        print "" # Newline after each row
    }
}
' extracted_metadata.txt > output.txt
