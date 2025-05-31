#!/bin/bash

# Define color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[0;33m'
NC='\033[0m' # No Color

# ==============================================
# OncoVI ClinVar RESOURCE MANAGER
# ==============================================
# The script manages ClinVar resource for OncoVI.
# First, the script check whether resources already
# exist, otherwise all the needed resources are 
# prepared for OncoVI analyisis.

echo ""
echo -e "${BLUE}==================================================================${NC}"
echo -e "${BLUE}===== OncoVI Resource Manager ====================================${NC}"
echo ""

# Create project folder structure
WORK_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

echo ""
echo -e "${BLUE}==================================================================${NC}"
echo -e "${BLUE}===== Working directory for OncoVI Project =======================${NC}"
echo ""
echo "${WORK_DIR}"
echo ""

# Create folder for resources 
RESOURCE_DIR=""$( dirname "${WORK_DIR}" )"/resources"

# Check if the resource folder already exists, otherwise create it
if [ ! -d "$RESOURCE_DIR" ]; then
    echo "Creating resources directory: $RESOURCE_DIR"
    mkdir -p "$RESOURCE_DIR"
else
    echo "Resources directory already exists"
fi

echo -e "${BLUE}==================================================================${NC}"
echo -e "${BLUE}===== Directory for OncoVI resources =============================${NC}"
echo ""
echo "${RESOURCE_DIR}"
echo ""

# Create folder for original resources
ORIGINAL_RESOURCE_DIR=""$( dirname "${WORK_DIR}" )"/original_resources"

# Check if the original resource folder already exists, otherwise create it
if [ ! -d "$ORIGINAL_RESOURCE_DIR" ]; then
    echo "Creating original resources directory: $ORIGINAL_RESOURCE_DIR"
    mkdir -p "$ORIGINAL_RESOURCE_DIR"
else
    echo "Original resources directory already exists"
fi

echo -e "${BLUE}==================================================================${NC}"
echo -e "${BLUE}===== Directory for OncoVI original resources ====================${NC}"
echo ""
echo "${ORIGINAL_RESOURCE_DIR}"
echo ""

# # Display help of the script
# print_help() {
    # echo ""
    # echo "The script prepares ClinVar resources for OncoVI"
    # echo ""
    # echo "How to use: $0 [COMMAND]"
    # echo ""
    # echo "Commands:"
    # echo "  cosmic_prepare Prepare resources from downloaded files"
    # echo "  cosmic_info    Show information about COSMIC files"
    # echo "  help           Display this help message"
    # echo ""
# }

# Check whether ClinVar underlying resources for OncoVI are available
check_ClinVar_resources() {
    echo -e "${BLUE}==================================================================${NC}"
    echo -e "${BLUE}===== Check ClinVar underlying resources for OncoVI ==============${NC}"
    echo ""
    
    # Check if ClinVar resource exists
    echo ""
    echo "Check if ClinVar resource exists:"

    REQUIRED_FILE=("clinvar_all_dictionary.txt")

    CLINVAR_RES="yes"  # Assume everything is present
    
    for file in "${REQUIRED_FILE[@]}"; do
        if [ ! -f "$RESOURCE_DIR/$file" ]; then
            echo ""
            echo "Required file is missing: $file"

            CLINVAR_RES="no"

            break  # Optional: exit loop early on first missing file
        else
            echo "Required file is found: $file"
            
        fi
    done

    # If ClinVar resource does not exist download it
    if [ "$CLINVAR_RES" = "no" ]; then
        
        # Check if original resources directory already exists
        if [ -d "$ORIGINAL_RESOURCE_DIR" ]; then
            echo ""
            cd "$ORIGINAL_RESOURCE_DIR"
            echo "Download ClinVar underlying resources:"
            echo ""

            wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
            
            echo "ClinVar resource successfully downloaded!"
        else 
            echo "Original resources directory does not exist"
        fi
    else
        echo "ClinVar resource is present"
        
    fi
	
	echo ""
}

# Prepare ClinVar underlying resources for OncoVI
prepare_ClinVar_resources() {
    echo -e "${BLUE}==================================================================${NC}"
    echo -e "${BLUE}===== Prepare ClinVar underlying resource for OncoVI =============${NC}"
    echo ""

    # Check if downloaded resource from ClinVar exists
    echo "Check if ClinVar downloaded resource exists:"
    echo ""
    
    DOWNLOADED_CLINVAR_FILE=("variant_summary.txt.gz")

    DOWNLOADED_CLINVAR_RES="yes"  # Assume everything is present
	
	# Check if file exists
    for file in "${DOWNLOADED_CLINVAR_FILE[@]}"; do
        if [ ! -f "$ORIGINAL_RESOURCE_DIR/$file" ]; then
            echo ""
            echo "Downloaded file from ClinVar is missing: $file"
            
	    DOWNLOADED_CLINVAR_RES="no"
            
	    break  # Optional: exit loop early on first missing file
        else
            echo "Downloaded file is found: $file"

        fi
    done

    # If downlaoded resource from ClinVar exists, process it 
    if [ "$DOWNLOADED_CLINVAR_RES" = "yes" ]; then

        # Prepare downloaded resource from ClinVar for OncoVI 
        echo ""
        echo "Preparing downloaded resource from ClinVar for OncoVI"
        echo ""
	echo "Select variants in genome reference build GRCh38:"
	echo ""

	zcat "$ORIGINAL_RESOURCE_DIR/$file" | awk -F '\t' '$17 == "GRCh38"' > variant_summary_GRCh38.txt

	echo "Reduce ClinVar resource to columns of interest:"
	echo ""
	awk -F'\t' -v OFS='|'  '{ print $5,$7,$19,$25,$26,$31,$32,$33,$34,$35,$37,$38,$40 }' variant_summary_GRCh38.txt > variant_summary_GRCh38_red.txt

	CLINVAR_GRCH38_VAR_RED="$ORIGINAL_RESOURCE_DIR/variant_summary_GRCh38_red.txt"

	echo "Create clinvar_all_dictionary.txt"
	echo ""

	# Run the python script to create the clinvar_all_dictionary.txt file
	python ${WORK_DIR}/create_clinvar_dict.py -i "$CLINVAR_GRCH38_VAR_RED" -r "$RESOURCE_DIR"
        
        # Remove intermediate files
        echo "Removing intermediate files from $ORIGINAL_RESOURCE_DIR"
        echo ""

        rm $ORIGINAL_RESOURCE_DIR/variant_summary_GRCh38.txt $CLINVAR_GRCH38_VAR_RED
    else
        echo "Downlaoded resource from ClinVar is not present"
       
    fi
}

# ===== MAIN PROGRAM =====
# Process command line arguments
print_help
check_ClinVar_resources
prepare_ClinVar_resources
