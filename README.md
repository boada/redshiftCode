# Redshift-Finder
Copy the python script somewhere in your PATH

e.g. make a $HOME/bin dir and then edit .bashrc to include:

export PATH=$HOME/bin:$PATH

and copy the script there.

Copy the templates directory somewhere, and edit the values of tempDir and 
tremontiFileName accordingly, so that the code can find them.


The code needs spectra in FITS table format. The script iraf_spec1dToFITSTable.py 
shows an example of how to do that.

Run with e.g.:

visualTemplateRedshift5.py spec1d*.fits results

assuming you're in a directory containing FITS table format spectra called 
spec1d*.fits

Plots and text file logging redshift measurements will be written under results/

An example spectrum in FITS table format is included to play with.

## Authors
- Matt Hilton (University of KwaZulu-Natal)
- Steven Boada (Texas A&M University)

## License
**Copyright 2015 the authors**.

Currently released under the **MIT license**.
See the file `LICENSE` for more details.

Any use of the content of this repository in academic research comes with a **request of acknowlegement**.
Contact the authors if it is unclear how to properly cite this work.

