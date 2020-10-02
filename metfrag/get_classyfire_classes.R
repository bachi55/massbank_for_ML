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

require(classyfireR)

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
    stop("Usage: R get_classyfire_classes INPUT_FILE OUTPUT_FILE")
}

# Read the molecule information
info <- read.csv(args[1])

pubchem_ids <- info$pubchem_id
inchikeys <- info$inchikey

# Set up output table
out <- data.frame(pubchem_id=pubchem_ids, inchikey=inchikeys,
                  kindom=NA, superclass=NA, class=NA, subclass=NA,
                  level_5=NA, level_6=NA, level_7=NA, level_8=NA)

# Get Classyfire results
for (idx in seq_along(inchikeys)) {
    inchikey <- inchikeys[idx]
    res <- get_classification(inchikey)
    if (is.null(res)) { next }
    out[idx, 3:10] <- unclass(classification(res))$Classification[1:8]
}

# Write out classifications
write.csv(out, file=args[2], row.names=FALSE, quote=TRUE)
