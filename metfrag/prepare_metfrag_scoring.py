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
import argparse
import os
import pandas as pd
import gzip

from massbank2db.db import MassbankDB
from massbank2db.spectrum import MBSpectrum


if __name__ == "__main__":
    argparser = argparse.ArgumentParser()
    argparser.add_argument("ion_mode", type=str, choices=["negative", "positive"],
                           help="Only datasets with the specified ionization mode are considered.")
    argparser.add_argument("--pubchemdb", type=str, help="Path to the PubChem SQLite database.",
                           default="/run/media/bach/EVO500GB/data/pubchem_24-06-2019/db/pubchem.sqlite")
    argparser.add_argument("--massbankdb", type=str, help="Path to the MassBank SQLite database.",
                           default="../massbank__v0_4_2.sqlite")
    args = argparser.parse_args()

    with MassbankDB(args.massbankdb) as mb_db:
        ds_df = mb_db.get_datasets_table()
        dss = ds_df[ds_df.ion_mode == args.ion_mode]["name"].values

        for idx, ds in enumerate(dss):
            print("DS: %s (%d/%d)" % (ds, idx + 1, len(dss)))

            os.makedirs(os.path.join("results", ds), exist_ok=True)

            mol_rt_df = []

            for mol, specs, cands in mb_db.iter_spectra(dataset=ds, return_candidates="mz", pc_dbfn=args.pubchemdb):
                # Spectra peaks are merged into a single spectrum.
                spec = MBSpectrum.merge_spectra(specs)

                # Output the merged spectrum in MetFrag format
                try:
                    metfrag_output = spec._to_metfrag_format(
                        cands=cands,
                        **{"MetFragScoreWeights": [1.0],
                           "MetFragScoreTypes": ["FragmenterScore"],
                           "LocalDatabasePath": os.path.join("..", "..", "candidates", ds),
                           "ResultsPath": os.path.join("..", "..", "results", ds),
                           "NumberThreads": 4,
                           "PeakListPath": os.path.join("..", "..", "peaks", ds),
                           "UseSmiles": True
                           }
                    )
                except ValueError as err:
                    print("Skip:", spec.get("original_accessions"))
                    print(err)

                # Collect information about the ground truth structure and retention time
                mol_rt_df.append([spec.get("accession"),
                                  ";".join(spec.get("original_accessions")),
                                  spec.get("pubchem_id"),
                                  spec.get("smiles_iso"),
                                  spec.get("retention_time"),
                                  spec.get("retention_time_unit"),
                                  spec.get("inchikey")])

                # Write out the configuration, peak list and candidates
                for k, v in metfrag_output.items():
                    if ".conf" in k:
                        odir = os.path.join("configs", ds)
                        ofn = os.path.join(odir, k)
                        opener = lambda fn: open(fn, "w")
                    elif ".peaks" in k:
                        odir = os.path.join("peaks", ds)
                        ofn = os.path.join(odir, k)
                        opener = lambda fn: open(fn, "w")
                    elif ".cands" in k:
                        odir = os.path.join("candidates", ds)
                        ofn = os.path.join(odir, k + ".gz")
                        opener = lambda fn: gzip.open(fn, "wt")
                    else:
                        ValueError("UPS")

                    os.makedirs(odir, exist_ok=True)
                    with opener(ofn) as ofile:
                        ofile.write(v)

            # Write out information about the ground truth structure and retention time
            os.makedirs(os.path.join("mol_rt_info"), exist_ok=True)
            pd.DataFrame(mol_rt_df,
                         columns=["accession", "original_accessions", "pubchem_id", "smiles_iso", "rt", "rt_unit",
                                  "inchikey"]) \
                .to_csv(os.path.join("mol_rt_info", ds + ".csv"), index=False)
