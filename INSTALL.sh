#!/usr/bin/env sh

# Check pipeline dependency
tool_desc_ver="1.5"

# Check Nextflow
if (command -v nextflow > /dev/null 2>&1); then {
	echo "Nextflow installed"
} else {
  		echo "Please install Nextflow by using:
curl -s https://get.nextflow.io | bash 
"
}
fi
# Check Singularity
if (command -v singularity > /dev/null 2>&1); then {
	echo "Singularity installed"
} else {
  		echo "Please ask your IT to install Singularity:
Go here https://singularity.lbl.gov/install-linux
"
}
fi

#Check Tool description for multiQC
if [ -e "conf_tools.txt" ]; then {
        echo "make_tool_desc_for_multiqc installed"
} else {
                echo "make_tool_desc_for_multiqc not installed. Installing it..."
                wget https://github.com/CRG-CNAG/make_tool_desc_for_multiqc/archive/v${tool_desc_ver}.tar.gz 
                tar -zvxf v${tool_desc_ver}.tar.gz
		mv make_tool_desc_for_multiqc-${tool_desc_ver}/conf_tools.txt .
                rm -fr make_tool_desc_for_multiqc-${tool_desc_ver} 
		rm -f v${tool_desc_ver}.tar.gz
}
fi



