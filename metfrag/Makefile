####
#
# The MIT License (MIT)
#
# Copyright 2020 Eric Bach <eric.bach@aalto.fi>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is furnished
# to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
####

# Path to the MetFrag jar-file
METFRAG_JAR:=~/.local/bin/MetFrag2.4.5.jar
METFRAG_NICE_VALUE:=19

# Input and output directories
CONF_DIR:=configs/
PEAK_DIR:=peaks/
CAND_DIR:=candidates/
RES_DIR:=results/
CLASS_DIR:=classyfire/
MOL_RT_INFO_DIR:=mol_rt_info/

# Add search paths for the recipes
VPATH:=$(CONF_DIR) $(PEAK_DIR) $(CAND_DIR)

INPUTS=$(subst $(CONF_DIR),,$(wildcard $(CONF_DIR)*/*.conf))  # config/AU_001/AU039212.conf -> AU_001/AU039212.conf
OUTPUTS_CANDIDATES=$(INPUTS:.conf=.csv.gz)                    # AU_001/AU039212.conf -> AU_001/AU039212.csv.gz
OUTPUTS_CLASSYFIRE:=$(addsuffix .classes,$(subst $(CONF_DIR),,$(wildcard $(CONF_DIR)*)))
                                                              # config/AU_001 -> AU_001.classes

# Multiple (merged) records cannot be processed by MetFrag (due to 'GC overhead limit exceeded'). We delete those
# records.
RECORDS_TO_DELETE:=BS_000/BS48628730

default: all

.PHONY: all scoring clean verify

all: scoring

scoring: remove_records verify $(addprefix $(RES_DIR), $(OUTPUTS_CANDIDATES))

verify:
	$(eval DATASETS := $(basename $(notdir $(wildcard $(MOL_RT_INFO_DIR)/*.csv))))
	$(eval NSPEC := $(foreach DS, $(DATASETS), $(words $(wildcard $(PEAK_DIR)/$(DS)/*.peaks))))
	$(eval NCAND := $(foreach DS, $(DATASETS), $(words $(wildcard $(CAND_DIR)/$(DS)/*.cands.gz))))
	$(eval NCONF := $(foreach DS, $(DATASETS), $(words $(wildcard $(CONF_DIR)/$(DS)/*.conf))))
	$(eval NINFO := $(foreach FILE, $(wildcard $(MOL_RT_INFO_DIR)/*.csv),\
		$(words $(wordlist 2, $(words $(file <$(FILE))), $(file <$(FILE))))))
	@echo "SPECTRA   :" $(NSPEC)
	@echo "CANDIDATES:" $(NCAND)
	@echo "CONFIGS   :" $(NCONF)
	@echo "INFO      :" $(NINFO)
	$(if $(findstring $(NSPEC), $(NCAND)),,$(error Number of spectra != Number of candidates))
	$(if $(findstring $(NSPEC), $(NCONF)),,$(error Number of spectra != Number of configurations))
	$(if $(findstring $(NSPEC), $(NINFO)),,$(error Number of spectra != Number of information))

scoring_progress:
	$(eval DATASETS := $(basename $(notdir $(wildcard $(MOL_RT_INFO_DIR)/*.csv))))
	$(eval NSPEC := $(foreach DS, $(DATASETS), $(words $(wildcard $(PEAK_DIR)/$(DS)/*.peaks))))
	$(eval NRES := $(foreach DS, $(DATASETS), $(words $(wildcard $(RES_DIR)/$(DS)/*.csv.gz))))
	@echo "DATASETS:" $(DATASETS)
	@echo "SPECTRA :" $(NSPEC)
	@echo "RESULTS :" $(NRES)

remove_records:
	# Delete (merged) records that cannot be processed by MetFrag
	$(foreach RECORD,$(notdir $(RECORDS_TO_DELETE)),find . -name $(notdir $(RECORD))* -delete;)
	# Remove those also from the mol_rt_info tables
	$(foreach RECORD,$(RECORDS_TO_DELETE),sed -i '/$(notdir $(RECORD))/d' $(MOL_RT_INFO_DIR)/$(subst /,,$(dir $(RECORD))).csv;)
	# Remove those also from the classyfire tables
	$(foreach RECORD,$(RECORDS_TO_DELETE),sed -i '/$(notdir $(RECORD))/d' $(CLASS_DIR)/$(subst /,,$(dir $(RECORD))).classes;)

clean:
	rm -f $(addsuffix .gz, $(addprefix $(RES_DIR), $(OUTPUTS_CANDIDATES)))
	rm -rf classyfire

$(RES_DIR)%.csv.gz: %.conf %.peaks %.cands.gz
	@echo "========"
	@echo "Scoring:" $<
	gunzip -fk $(word 3, $^)
	nice -$(METFRAG_NICE_VALUE) java -jar $(METFRAG_JAR) \
		ParameterFile=$(word 1, $^) \
		PeakListPath=$(word 2, $^) \
		LocalDatabasePath=$(subst .gz,,$(word 3, $^)) \
		ResultsPath=$(@D)
	rm -f $(subst .gz,, $(word 3, $^))
	@echo "---------"
	@echo "Simplify:" $(basename $@)
	python simplify_metfrag_output.py $(basename $@)
	gzip $(basename $@)

$(CLASS_DIR)%.classes: $(MOL_RT_INFO_DIR)/%.csv
	mkdir -p classyfire
	Rscript get_classyfire_classes.R $< $@

classyfire: remove_records verify $(addprefix $(CLASS_DIR), $(OUTPUTS_CLASSYFIRE))