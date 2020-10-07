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
RECORDS_TO_DELETE:=BML_000/BML05865 BS_000/BS483105 BS_000/BS576728 BS_000/BS629471 BS_000/BS679807 BS_000/BS790201 \
	BS_000/BS808533 EQ_000/EQ064253 EQ_000/EQ067784 EQ_000/EQ164898 EQ_000/EQ322405 EQ_000/EQ357642 EQ_000/EQ357642 \
	EQ_000/EQ445376 EQ_000/EQ636351 EQ_001/EQ032906 EQ_001/EQ055269 EQ_001/EQ110751 PR_000/PR127009 PR_000/PR312201 \
	PR_000/PR312648

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

remove_records:
	# Delete (merged) records that cannot be processed by MetFrag
	$(foreach RECORD,$(notdir $(RECORDS_TO_DELETE)),find . -name $(notdir $(RECORD))* -delete;)
	# Remove those also from the mol_rt_info tables
	$(foreach RECORD,$(RECORDS_TO_DELETE),sed -i '/$(notdir $(RECORD))/d' $(MOL_RT_INFO_DIR)/$(subst /,,$(dir $(RECORD))).csv;)

clean:
	rm -f $(addsuffix .gz, $(addprefix $(RES_DIR), $(OUTPUTS_CANDIDATES)))
	rm -rf classyfire

$(RES_DIR)%.csv.gz: %.conf %.peaks %.cands.gz
	@echo "========"
	@echo "Scoring:" $<
	gunzip -fk $(word 3, $^)
	java -jar $(METFRAG_JAR) \
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